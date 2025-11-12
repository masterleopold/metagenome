#!/usr/bin/env python3
"""
Trigger host removal phase on EC2
"""

import json
import boto3
import os
import re
import shlex
from typing import Dict, Any

ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')

ANALYSIS_AMI = os.environ['ANALYSIS_AMI_ID']
INSTANCE_TYPE = os.environ.get('HOST_REMOVAL_INSTANCE_TYPE', 'r5.2xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']
EFS_ID = os.environ['EFS_ID']


import re
import shlex

def validate_run_id(run_id: str) -> str:
    pattern = r'^[A-Z0-9][A-Z0-9_-]{0,63}$'
    if not re.match(pattern, run_id):
        raise ValueError(f"Invalid run_id format: {run_id}")
    return run_id

def validate_s3_path_component(component: str, name: str = "path") -> str:
    pattern = r'^[a-zA-Z0-9/_.-]+$'
    if not re.match(pattern, component) or '..' in component:
        raise ValueError(f"Invalid {name}: {component}")
    return component

def validate_bucket_name(bucket: str) -> str:
    pattern = r'^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$'
    if not re.match(pattern, bucket):
        raise ValueError(f"Invalid bucket: {bucket}")
    return bucket


def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Launch EC2 instance for host genome removal.
    """

    run_id = event['run_id']
    workflow_id = event['workflow_id']
    input_path = event['input_path']
    reference = event.get('config', {}).get('reference', 'sus_scrofa_11.1')

    # PMDA 4-virus support: check if RNA workflow needed
    target_viruses = event.get('config', {}).get('target_viruses', [])
    sample_type = event.get('config', {}).get('sample_type', 'dna')  # 'dna', 'rna', or 'dual'

    # Parse S3 path
    parts = input_path.replace('s3://', '').split('/', 1)
    bucket = parts[0]
    prefix = parts[1] if len(parts) > 1 else ''

    output_prefix = prefix.replace('qc', 'host_removal')

    # Launch instance with EFS mount for reference genome
    instance_id = launch_host_removal_instance(
        run_id, bucket, prefix, output_prefix
    )

    # Wait for instance
    waiter = ec2.get_waiter('instance_running')
    waiter.wait(InstanceIds=[instance_id])

    # Execute host removal (supports DNA, RNA, or DUAL workflows)
    response = execute_host_removal(
        instance_id, run_id, bucket, prefix, output_prefix, reference,
        sample_type, target_viruses
    )

    return {
        'phase': 'host_removal',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': response['Command']['CommandId'],
        'run_id': run_id,
        'workflow_id': workflow_id,
        'output_path': f's3://{bucket}/{output_prefix}'
    }

def launch_host_removal_instance(run_id: str, bucket: str,
                                input_prefix: str, output_prefix: str) -> str:
    """Launch EC2 instance with EFS mount."""

    user_data = f"""#!/bin/bash
# Mount EFS for reference genomes
apt-get update && apt-get install -y nfs-common
mkdir -p /mnt/efs
mount -t nfs4 {EFS_ID}.efs.{os.environ['AWS_REGION']}.amazonaws.com:/ /mnt/efs

export RUN_ID={shlex.quote(run_id)}
export S3_INPUT={shlex.quote(f's3://{bucket}/{input_prefix}')}
export S3_OUTPUT={shlex.quote(f's3://{bucket}/{output_prefix}')}

echo "Host removal instance started for run $RUN_ID" | logger

# Auto-terminate after 6 hours
echo "sudo shutdown -h +360" | at now + 6 hours
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
                    'VolumeSize': 300,
                    'VolumeType': 'gp3',
                    'DeleteOnTermination': True
                }
            }
        ],
        TagSpecifications=[
            {
                'ResourceType': 'instance',
                'Tags': [
                    {'Key': 'Name', 'Value': f'host-removal-{run_id}'},
                    {'Key': 'RunID', 'Value': run_id},
                    {'Key': 'Phase', 'Value': 'host_removal'}
                ]
            }
        ]
    )

    return response['Instances'][0]['InstanceId']

def execute_host_removal(instance_id: str, run_id: str, bucket: str,
                        input_prefix: str, output_prefix: str,
                        reference: str, sample_type: str = 'dna',
                        target_viruses: list = None) -> Dict:
    """Execute host removal on EC2 instance with DNA/RNA workflow support."""

    target_viruses = target_viruses or []
    target_viruses_str = ' '.join(target_viruses) if target_viruses else ''

    command = f"""
#!/bin/bash
set -euo pipefail

export RUN_ID={shlex.quote(run_id)}
WORK_DIR="/mnt/analysis/$RUN_ID/host_removal"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Download FASTQ files
aws s3 sync "s3://{bucket}/{input_prefix}" input/ --exclude "*" --include "*.fastq*"

# Determine input FASTQ
INPUT_FASTQ=$(find input/ -name "*.fastq*" | head -1)

# Use host removal orchestrator (supports DNA, RNA, DUAL workflows)
python3 /opt/minion/scripts/phase3_host_removal/host_removal_orchestrator.py \\
    -i "$INPUT_FASTQ" \\
    -o filtered/ \\
    -r "$RUN_ID" \\
    --sample-type {sample_type} \\
    --host-genome /mnt/efs/host_genome/{reference}.fa \\
    --rrna-db /mnt/efs/databases/host/pig_rrna/minimap2/pig_rrna.mmi \\
    -t 16

# Calculate depletion metrics (combined DNA+RNA if dual)
/opt/minion/scripts/phase3_host_removal/calculate_depletion.py \\
    --before input/ \\
    --after filtered/ \\
    --output depletion_stats.json \\
    --run-id "$RUN_ID"

# Upload results
aws s3 sync filtered/ "s3://{bucket}/{output_prefix}/" --exclude "*.sam" --exclude "*.bam"
aws s3 cp depletion_stats.json "s3://{bucket}/{output_prefix}/"

# Upload orchestrator results if available
if [ -f filtered/host_removal_orchestrator_results.json ]; then
    aws s3 cp filtered/host_removal_orchestrator_results.json "s3://{bucket}/{output_prefix}/"
fi

# Clean up
rm -rf "$WORK_DIR"
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={
            'commands': [command],
            'executionTimeout': ['21600']  # 6 hours
        }
    )

    return response