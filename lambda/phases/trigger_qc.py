#!/usr/bin/env python3
"""
Trigger QC phase on EC2
"""

import json
import boto3
import os
from typing import Dict, Any

ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')
lambda_client = boto3.client('lambda')

ANALYSIS_AMI = os.environ['ANALYSIS_AMI_ID']
INSTANCE_TYPE = os.environ.get('QC_INSTANCE_TYPE', 'm5.xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Launch EC2 instance for QC analysis.
    """

    run_id = event['run_id']
    workflow_id = event['workflow_id']
    input_path = event['input_path']
    min_quality = event.get('config', {}).get('min_quality', 9)

    # Parse bucket and prefix from input path
    parts = input_path.replace('s3://', '').split('/', 1)
    bucket = parts[0]
    prefix = parts[1] if len(parts) > 1 else ''

    output_prefix = prefix.replace('basecalling', 'qc')

    # Launch EC2 instance
    instance_id = launch_analysis_instance(
        run_id, 'qc', bucket, prefix, output_prefix
    )

    # Wait for instance to be ready
    waiter = ec2.get_waiter('instance_running')
    waiter.wait(InstanceIds=[instance_id])

    # Execute QC scripts
    response = execute_qc_commands(
        instance_id, run_id, bucket, prefix, output_prefix, min_quality
    )

    command_id = response['Command']['CommandId']

    return {
        'phase': 'qc',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': command_id,
        'run_id': run_id,
        'workflow_id': workflow_id,
        'output_path': f's3://{bucket}/{output_prefix}'
    }

def launch_analysis_instance(run_id: str, phase: str, bucket: str,
                            input_prefix: str, output_prefix: str) -> str:
    """Launch EC2 instance for analysis."""

    user_data = f"""#!/bin/bash
export RUN_ID='{run_id}'
export PHASE='{phase}'
export S3_INPUT='s3://{bucket}/{input_prefix}'
export S3_OUTPUT='s3://{bucket}/{output_prefix}'

echo "Analysis instance started for $PHASE phase of run $RUN_ID" | logger

# Auto-terminate after 4 hours
echo "sudo shutdown -h +240" | at now + 4 hours
"""

    response = ec2.run_instances(
        ImageId=ANALYSIS_AMI,
        InstanceType=INSTANCE_TYPE,
        MinCount=1,
        MaxCount=1,
        SubnetId=SUBNET_ID,
        SecurityGroupIds=[SECURITY_GROUP],
        IamInstanceProfile={'Name': 'MinIONEC2Role'},
        UserData=user_data,
        BlockDeviceMappings=[
            {
                'DeviceName': '/dev/xvda',
                'Ebs': {
                    'VolumeSize': 200,
                    'VolumeType': 'gp3',
                    'DeleteOnTermination': True
                }
            }
        ],
        TagSpecifications=[
            {
                'ResourceType': 'instance',
                'Tags': [
                    {'Key': 'Name', 'Value': f'{phase}-{run_id}'},
                    {'Key': 'RunID', 'Value': run_id},
                    {'Key': 'Phase', 'Value': phase},
                    {'Key': 'AutoTerminate', 'Value': 'true'}
                ]
            }
        ]
    )

    return response['Instances'][0]['InstanceId']

def execute_qc_commands(instance_id: str, run_id: str, bucket: str,
                       input_prefix: str, output_prefix: str,
                       min_quality: int) -> Dict:
    """Execute QC analysis on EC2 instance."""

    command = f"""
#!/bin/bash
set -euo pipefail

# Set variables
export RUN_ID='{run_id}'
export S3_INPUT='s3://{bucket}/{input_prefix}'
export S3_OUTPUT='s3://{bucket}/{output_prefix}'

# Create working directory
mkdir -p /mnt/analysis/$RUN_ID/qc
cd /mnt/analysis/$RUN_ID/qc

# Download input files
aws s3 sync "$S3_INPUT" input/ --exclude "*.fast5" --exclude "*.pod5"

# Run NanoPlot for QC
/opt/minion/scripts/phase2_qc/nanoplot_qc.sh \\
    -i input/ \\
    -o qc_report/ \\
    -r "$RUN_ID"

# Check QC metrics
/opt/minion/scripts/phase2_qc/check_qc_metrics.py \\
    --input qc_report/NanoStats.txt \\
    --output qc_summary.json \\
    --min-quality {min_quality} \\
    --run-id "$RUN_ID"

# Check if QC passed
if grep -q '"qc_status": "PASS"' qc_summary.json; then
    echo "QC PASSED"
else
    echo "QC FAILED - check metrics"
    # Upload results anyway for troubleshooting
fi

# Upload results
aws s3 sync qc_report/ "$S3_OUTPUT" --exclude "*.fast5"
aws s3 cp qc_summary.json "$S3_OUTPUT"

# Clean up
rm -rf /mnt/analysis/$RUN_ID/qc
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={
            'commands': [command],
            'executionTimeout': ['14400']  # 4 hours
        }
    )

    return response