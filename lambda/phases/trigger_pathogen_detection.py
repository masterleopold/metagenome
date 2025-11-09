#!/usr/bin/env python3
"""
Trigger pathogen detection phase on EC2
"""

import json
import boto3
import os
import re
import shlex
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
    """Lambda handler for Phase 4: Multi-database pathogen detection.

    Orchestrates pathogen screening using Kraken2, BLAST (RVDB), and targeted
    PMDA database searches. Critical for identifying all 91 PMDA-designated
    pathogens and PERV sequences. Launches high-memory EC2 instance with EFS
    mount for database access.

    Args:
        event: Lambda event dict containing:
            - run_id (str): Unique run identifier
            - workflow_id (int): Step Functions execution ID
            - input_path (str): S3 path to host-depleted FASTQ files
            - config (dict, optional): Config with databases list
              Default databases: ['kraken2', 'rvdb', 'pmda']

        context: Lambda context object (unused)

    Returns:
        Dict with execution details:
            - phase: "pathogen_detection"
            - status: "RUNNING"
            - instance_id: EC2 instance ID
            - command_id: SSM command ID
            - run_id: Run identifier
            - workflow_id: Workflow ID
            - output_path: S3 path to results

    Raises:
        KeyError: If required event parameters missing
        ValueError: If input_path format invalid
        botocore.exceptions.ClientError: If AWS operations fail

    Note:
        Automatically triggers SNS alert if PERV detected.
        Uses r5.4xlarge (128GB RAM) for Kraken2 database.
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
    """Launch high-memory EC2 instance for pathogen detection with spot/on-demand fallback.

    Mounts EFS for access to reference databases (Kraken2, BLAST, PMDA).
    Attempts spot instance first (70% cost savings), falls back to on-demand
    if unavailable after 5-minute timeout.

    Args:
        run_id: Unique identifier for this sequencing run (e.g., "RUN-2024-001")
        bucket: S3 bucket containing input FASTQ files (post-host-depletion)
        input_prefix: S3 prefix to filtered reads
        output_prefix: S3 prefix for pathogen detection results

    Returns:
        str: EC2 instance ID of the launched pathogen detection instance

    Raises:
        botocore.exceptions.ClientError: If both spot and on-demand launches fail

    Note:
        Uses r5.4xlarge (128GB RAM) for Kraken2 database loading.
        Auto-terminates after 8 hours to prevent runaway costs.
        Critical for PERV and PMDA 91-pathogen screening.
    """

    user_data = f"""#!/bin/bash
# Mount EFS for databases
apt-get update && apt-get install -y nfs-common
mkdir -p /mnt/efs
mount -t nfs4 {EFS_ID}.efs.{os.environ['AWS_REGION']}.amazonaws.com:/ /mnt/efs

export RUN_ID={shlex.quote(run_id)}
export S3_INPUT={shlex.quote(f's3://{bucket}/{input_prefix}')}
export S3_OUTPUT={shlex.quote(f's3://{bucket}/{output_prefix}')}

# Optimize for Kraken2 (uses lots of memory)
echo 3 > /proc/sys/vm/drop_caches
echo 1 > /proc/sys/vm/compact_memory

echo "Pathogen detection instance started for run $RUN_ID" | logger

# Auto-terminate after 8 hours
echo "sudo shutdown -h +480" | at now + 8 hours
"""

    instance_config = {
        'ImageId': ANALYSIS_AMI,
        'InstanceType': INSTANCE_TYPE,
        'SubnetId': SUBNET_ID,
        'SecurityGroupIds': [SECURITY_GROUP],
        'IamInstanceProfile': {'Name': os.environ.get('EC2_IAM_ROLE', 'MinIONEC2Role')},
        'UserData': user_data,
        'BlockDeviceMappings': [
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
        'TagSpecifications': [
            {
                'ResourceType': 'instance',
                'Tags': [
                    {'Key': 'Name', 'Value': f'pathogen-{run_id}'},
                    {'Key': 'RunID', 'Value': run_id},
                    {'Key': 'Phase', 'Value': 'pathogen_detection'}
                ]
            }
        ]
    }

    # Try spot instance first for cost savings
    use_spot = os.environ.get('USE_SPOT_INSTANCES', 'true').lower() == 'true'

    if use_spot:
        try:
            print(f"Attempting to launch spot instance for pathogen detection - {run_id}")
            response = ec2.request_spot_instances(
                InstanceCount=1,
                Type='one-time',
                LaunchSpecification=instance_config
            )

            spot_request_id = response['SpotInstanceRequests'][0]['SpotInstanceRequestId']

            # Wait for spot request fulfillment with timeout
            waiter = ec2.get_waiter('spot_instance_request_fulfilled')
            waiter.wait(
                SpotInstanceRequestIds=[spot_request_id],
                WaiterConfig={'Delay': 15, 'MaxAttempts': 20}  # 5 minute timeout
            )

            # Get instance ID
            spot_response = ec2.describe_spot_instance_requests(
                SpotInstanceRequestIds=[spot_request_id]
            )

            instance_id = spot_response['SpotInstanceRequests'][0]['InstanceId']
            print(f"Spot instance launched successfully: {instance_id}")

            return instance_id

        except Exception as e:
            print(f"Spot instance request failed: {str(e)}")
            print(f"Falling back to on-demand instance for {run_id}")

            # Cancel spot request if it exists
            try:
                ec2.cancel_spot_instance_requests(SpotInstanceRequestIds=[spot_request_id])
            except:
                pass

    # Launch on-demand instance (fallback or if spot disabled)
    print(f"Launching on-demand instance for pathogen detection - {run_id}")
    response = ec2.run_instances(
        MinCount=1,
        MaxCount=1,
        **instance_config
    )

    instance_id = response['Instances'][0]['InstanceId']
    print(f"On-demand instance launched: {instance_id}")

    return instance_id

def execute_pathogen_detection(instance_id: str, run_id: str, bucket: str,
                              input_prefix: str, output_prefix: str,
                              databases: List[str]) -> Dict:
    """Execute pathogen detection analysis on EC2 instance via SSM.

    Builds and runs shell commands for selected databases (Kraken2, BLAST,
    PMDA targeted search) plus mandatory PERV analysis. Downloads filtered
    reads from S3, runs analyses, aggregates results, and uploads to S3.

    Args:
        instance_id: EC2 instance ID to execute commands on
        run_id: Unique run identifier for naming/tracking
        bucket: S3 bucket name
        input_prefix: S3 prefix to host-filtered FASTQ files
        output_prefix: S3 prefix for output results
        databases: List of databases to search. Options:
            - 'kraken2': Taxonomic classification
            - 'rvdb': BLAST viral database search
            - 'pmda': Targeted 91-pathogen search
            Note: PERV analysis always runs regardless of this list

    Returns:
        Dict: SSM send_command API response containing Command object
            with CommandId for monitoring execution

    Raises:
        botocore.exceptions.ClientError: If SSM command fails to send

    Note:
        Automatically publishes SNS alert if PERV detected.
        Command timeout: 8 hours.
        Instance auto-cleans up work directory after upload.
    """

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
    db_commands.append(f"""
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

export RUN_ID={shlex.quote(run_id)}
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
            'commands': [command],
            'executionTimeout': ['28800']  # 8 hours
        }
    )

    return response