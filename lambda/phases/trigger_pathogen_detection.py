#!/usr/bin/env python3
"""
Trigger pathogen detection phase on EC2
"""

import json
import boto3
import os
from typing import Dict, Any, List

ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')
sns = boto3.client('sns')

ANALYSIS_AMI = os.environ['ANALYSIS_AMI_ID']
INSTANCE_TYPE = os.environ.get('PATHOGEN_INSTANCE_TYPE', 'r5.4xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']
EFS_ID = os.environ['EFS_ID']
SNS_TOPIC = os.environ['SNS_TOPIC_ARN']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Launch EC2 instance for pathogen detection.
    """

    run_id = event['run_id']
    workflow_id = event['workflow_id']
    input_path = event['input_path']
    databases = event.get('config', {}).get('databases', ['kraken2', 'rvdb', 'pmda'])

    # Parse S3 path
    parts = input_path.replace('s3://', '').split('/', 1)
    bucket = parts[0]
    prefix = parts[1] if len(parts) > 1 else ''

    output_prefix = prefix.replace('host_removal', 'pathogen_detection')

    # Launch high-memory instance for pathogen detection
    instance_id = launch_pathogen_instance(
        run_id, bucket, prefix, output_prefix
    )

    # Wait for instance
    waiter = ec2.get_waiter('instance_running')
    waiter.wait(InstanceIds=[instance_id])

    # Execute pathogen detection
    response = execute_pathogen_detection(
        instance_id, run_id, bucket, prefix, output_prefix, databases
    )

    return {
        'phase': 'pathogen_detection',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': response['Command']['CommandId'],
        'run_id': run_id,
        'workflow_id': workflow_id,
        'output_path': f's3://{bucket}/{output_prefix}'
    }

def launch_pathogen_instance(run_id: str, bucket: str,
                            input_prefix: str, output_prefix: str) -> str:
    """Launch high-memory EC2 instance."""

    user_data = f"""#!/bin/bash
# Mount EFS for databases
apt-get update && apt-get install -y nfs-common
mkdir -p /mnt/efs
mount -t nfs4 {EFS_ID}.efs.{os.environ['AWS_REGION']}.amazonaws.com:/ /mnt/efs

export RUN_ID='{run_id}'
export S3_INPUT='s3://{bucket}/{input_prefix}'
export S3_OUTPUT='s3://{bucket}/{output_prefix}'

# Optimize for Kraken2 (uses lots of memory)
echo 3 > /proc/sys/vm/drop_caches
echo 1 > /proc/sys/vm/compact_memory

echo "Pathogen detection instance started for run $RUN_ID" | logger

# Auto-terminate after 8 hours
echo "sudo shutdown -h +480" | at now + 8 hours
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
                    'VolumeSize': 500,
                    'VolumeType': 'gp3',
                    'Iops': 10000,
                    'DeleteOnTermination': True
                }
            }
        ],
        TagSpecifications=[
            {
                'ResourceType': 'instance',
                'Tags': [
                    {'Key': 'Name', 'Value': f'pathogen-{run_id}'},
                    {'Key': 'RunID', 'Value': run_id},
                    {'Key': 'Phase', 'Value': 'pathogen_detection'}
                ]
            }
        ]
    )

    return response['Instances'][0]['InstanceId']

def execute_pathogen_detection(instance_id: str, run_id: str, bucket: str,
                              input_prefix: str, output_prefix: str,
                              databases: List[str]) -> Dict:
    """Execute pathogen detection on EC2 instance."""

    # Build database commands
    db_commands = []
    if 'kraken2' in databases:
        db_commands.append("""
# Run Kraken2
/opt/minion/scripts/phase4_pathogen/kraken2_search.sh \\
    -i filtered/ \\
    -o kraken2/ \\
    -d /mnt/efs/databases/kraken2/standard \\
    -r "$RUN_ID"

# Extract PMDA pathogens
/opt/minion/scripts/phase4_pathogen/extract_pmda_pathogens.py \\
    --kraken kraken2/report.txt \\
    --output kraken2/pmda_pathogens.json \\
    --run-id "$RUN_ID"
""")

    if 'rvdb' in databases:
        db_commands.append("""
# Run BLAST against RVDB
/opt/minion/scripts/phase4_pathogen/blast_search.sh \\
    -i filtered/ \\
    -o blast/ \\
    -d /mnt/efs/databases/rvdb/rvdb.fasta \\
    -r "$RUN_ID"
""")

    if 'pmda' in databases:
        db_commands.append("""
# Run targeted PMDA pathogen search
/opt/minion/scripts/phase4_pathogen/pmda_targeted_search.py \\
    --input filtered/ \\
    --output pmda/ \\
    --database /mnt/efs/databases/pmda/pmda_91.fasta \\
    --run-id "$RUN_ID"
""")

    # Always run PERV analysis (critical for xenotransplantation)
    db_commands.append("""
# PERV-specific analysis
/opt/minion/scripts/phase4_pathogen/perv_analysis.sh \\
    -i filtered/ \\
    -o perv/ \\
    -r "$RUN_ID"

# Check for PERV detection
if grep -q '"perv_detected": true' perv/perv_summary.json; then
    echo "CRITICAL: PERV detected in sample!"
    aws sns publish \\
        --topic-arn {SNS_TOPIC} \\
        --subject "CRITICAL: PERV Detection - Run $RUN_ID" \\
        --message "PERV sequences detected in run $RUN_ID. Immediate review required."
fi
""")

    command = f"""
#!/bin/bash
set -euo pipefail

export RUN_ID='{run_id}'
WORK_DIR="/mnt/analysis/$RUN_ID/pathogen"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Download filtered reads
aws s3 sync "s3://{bucket}/{input_prefix}" filtered/ --exclude "*" --include "*.fastq*"

# Get read statistics
echo "Total reads after host removal: $(cat filtered/*.fastq | wc -l)/4 | bc" > stats.txt

{''.join(db_commands)}

# Aggregate results
/opt/minion/scripts/phase4_pathogen/aggregate_results.py \\
    --kraken kraken2/ \\
    --blast blast/ \\
    --perv perv/ \\
    --output pathogen_summary.json \\
    --run-id "$RUN_ID"

# Upload all results
aws s3 sync . "s3://{bucket}/{output_prefix}/" \\
    --exclude "filtered/*" \\
    --exclude "*.fastq*"

# Clean up
rm -rf "$WORK_DIR"
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={
            'commands': [command.format(SNS_TOPIC=SNS_TOPIC)],
            'executionTimeout': ['28800']  # 8 hours
        }
    )

    return response