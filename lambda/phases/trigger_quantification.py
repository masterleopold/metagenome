#!/usr/bin/env python3
"""
Trigger quantification phase on EC2
"""

import json
import boto3
import os
from typing import Dict, Any

ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')

ANALYSIS_AMI = os.environ['ANALYSIS_AMI_ID']
INSTANCE_TYPE = os.environ.get('QUANT_INSTANCE_TYPE', 'm5.2xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Launch EC2 instance for pathogen quantification.
    """

    run_id = event['run_id']
    workflow_id = event['workflow_id']
    input_path = event['input_path']
    spike_in = event.get('config', {}).get('spike_in', 'PhiX174')
    plasma_volume = event.get('config', {}).get('plasma_volume', 10.0)

    # Parse S3 path
    parts = input_path.replace('s3://', '').split('/', 1)
    bucket = parts[0]
    prefix = parts[1] if len(parts) > 1 else ''

    output_prefix = prefix.replace('pathogen_detection', 'quantification')

    # Launch instance
    instance_id = launch_quant_instance(
        run_id, bucket, prefix, output_prefix
    )

    # Wait for instance
    waiter = ec2.get_waiter('instance_running')
    waiter.wait(InstanceIds=[instance_id])

    # Execute quantification
    response = execute_quantification(
        instance_id, run_id, bucket, prefix, output_prefix,
        spike_in, plasma_volume
    )

    return {
        'phase': 'quantification',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': response['Command']['CommandId'],
        'run_id': run_id,
        'workflow_id': workflow_id,
        'output_path': f's3://{bucket}/{output_prefix}'
    }

def launch_quant_instance(run_id: str, bucket: str,
                         input_prefix: str, output_prefix: str) -> str:
    """Launch EC2 instance for quantification."""

    user_data = f"""#!/bin/bash
export RUN_ID='{run_id}'
export S3_INPUT='s3://{bucket}/{input_prefix}'
export S3_OUTPUT='s3://{bucket}/{output_prefix}'

echo "Quantification instance started for run $RUN_ID" | logger

# Auto-terminate after 2 hours
echo "sudo shutdown -h +120" | at now + 2 hours
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
                    'VolumeSize': 100,
                    'VolumeType': 'gp3',
                    'DeleteOnTermination': True
                }
            }
        ],
        TagSpecifications=[
            {
                'ResourceType': 'instance',
                'Tags': [
                    {'Key': 'Name', 'Value': f'quantification-{run_id}'},
                    {'Key': 'RunID', 'Value': run_id},
                    {'Key': 'Phase', 'Value': 'quantification'}
                ]
            }
        ]
    )

    return response['Instances'][0]['InstanceId']

def execute_quantification(instance_id: str, run_id: str, bucket: str,
                          input_prefix: str, output_prefix: str,
                          spike_in: str, plasma_volume: float) -> Dict:
    """Execute quantification analysis on EC2 instance."""

    command = f"""
#!/bin/bash
set -euo pipefail

export RUN_ID='{run_id}'
WORK_DIR="/mnt/analysis/$RUN_ID/quantification"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Download pathogen detection results
aws s3 sync "s3://{bucket}/{input_prefix}" pathogen/ \\
    --exclude "*" --include "*.json" --include "*.txt"

# Kraken quantification
if [ -f pathogen/kraken2/report.txt ]; then
    /opt/minion/scripts/phase5_quantification/kraken_quantify.py \\
        --report pathogen/kraken2/report.txt \\
        --output kraken_quantification.json \\
        --run-id "$RUN_ID"
fi

# BLAST quantification
if [ -f pathogen/blast/results.tsv ]; then
    /opt/minion/scripts/phase5_quantification/blast_quantify.py \\
        --blast pathogen/blast/results.tsv \\
        --output blast_quantification.json \\
        --run-id "$RUN_ID"
fi

# Spike-in normalization
if [ -f pathogen/kraken2/report.txt ]; then
    /opt/minion/scripts/phase5_quantification/spike_in_normalization.py \\
        --report pathogen/kraken2/report.txt \\
        --spike-in "{spike_in}" \\
        --output normalized_quantification.json \\
        --run-id "$RUN_ID"
fi

# Calculate absolute copy numbers
/opt/minion/scripts/phase5_quantification/absolute_copy_number.py \\
    --kraken kraken_quantification.json \\
    --normalized normalized_quantification.json \\
    --plasma-volume {plasma_volume} \\
    --output absolute_quantification.json \\
    --run-id "$RUN_ID"

# Generate quantification report
cat > quantification_summary.json << EOF
{{
    "run_id": "$RUN_ID",
    "plasma_volume_ml": {plasma_volume},
    "spike_in": "{spike_in}",
    "timestamp": "$(date -Iseconds)",
    "files": {{
        "kraken": "kraken_quantification.json",
        "blast": "blast_quantification.json",
        "normalized": "normalized_quantification.json",
        "absolute": "absolute_quantification.json"
    }}
}}
EOF

# Upload results
aws s3 sync . "s3://{bucket}/{output_prefix}/" \\
    --exclude "pathogen/*"

# Clean up
rm -rf "$WORK_DIR"
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={
            'commands': [command],
            'executionTimeout': ['7200']  # 2 hours
        }
    )

    return response