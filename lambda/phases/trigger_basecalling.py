#!/usr/bin/env python3
"""
Trigger basecalling phase on EC2
"""

import json
import boto3
import os
from typing import Dict, Any

ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')
sns = boto3.client('sns')

SNS_TOPIC = os.environ['SNS_TOPIC_ARN']
BASECALLING_AMI = os.environ['BASECALLING_AMI_ID']
INSTANCE_TYPE = os.environ.get('BASECALLING_INSTANCE_TYPE', 'g4dn.xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Launch EC2 instance for basecalling and trigger analysis.
    """

    run_id = event['run_id']
    workflow_id = event['workflow_id']
    bucket = event['bucket']
    input_prefix = event['input_prefix']
    output_prefix = event['output_prefix']
    skip_duplex = event.get('config', {}).get('skip_duplex', False)

    # Launch EC2 instance
    instance_id = launch_basecalling_instance(
        run_id, bucket, input_prefix, output_prefix, skip_duplex
    )

    # Wait for instance to be ready
    waiter = ec2.get_waiter('instance_running')
    waiter.wait(InstanceIds=[instance_id])

    # Execute basecalling script via SSM
    response = execute_basecalling_command(
        instance_id, run_id, bucket, input_prefix, output_prefix, skip_duplex
    )

    command_id = response['Command']['CommandId']

    # Store instance and command info
    result = {
        'phase': 'basecalling',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': command_id,
        'run_id': run_id,
        'workflow_id': workflow_id
    }

    # Send notification
    send_notification(
        f'Basecalling Started - Run {run_id}',
        f'Basecalling phase initiated on instance {instance_id}\n'
        f'Input: s3://{bucket}/{input_prefix}\n'
        f'Duplex mode: {not skip_duplex}'
    )

    return result

def launch_basecalling_instance(run_id: str, bucket: str, input_prefix: str,
                               output_prefix: str, skip_duplex: bool) -> str:
    """Launch GPU instance for basecalling."""

    user_data = f"""#!/bin/bash
# Set environment
export RUN_ID='{run_id}'
export S3_INPUT='s3://{bucket}/{input_prefix}'
export S3_OUTPUT='s3://{bucket}/{output_prefix}'
export SKIP_DUPLEX='{str(skip_duplex).lower()}'

# Log startup
echo "Basecalling instance started for run $RUN_ID" | logger

# Auto-terminate after 12 hours
echo "sudo shutdown -h +720" | at now + 12 hours
"""

    # Use spot instance for cost savings
    response = ec2.request_spot_instances(
        InstanceCount=1,
        Type='one-time',
        LaunchSpecification={
            'ImageId': BASECALLING_AMI,
            'InstanceType': INSTANCE_TYPE,
            'SubnetId': SUBNET_ID,
            'SecurityGroupIds': [SECURITY_GROUP],
            'IamInstanceProfile': {
                'Name': 'MinIONEC2Role'
            },
            'UserData': user_data,
            'BlockDeviceMappings': [
                {
                    'DeviceName': '/dev/xvda',
                    'Ebs': {
                        'VolumeSize': 500,
                        'VolumeType': 'gp3',
                        'DeleteOnTermination': True
                    }
                }
            ],
            'TagSpecifications': [
                {
                    'ResourceType': 'instance',
                    'Tags': [
                        {'Key': 'Name', 'Value': f'basecalling-{run_id}'},
                        {'Key': 'RunID', 'Value': run_id},
                        {'Key': 'Phase', 'Value': 'basecalling'},
                        {'Key': 'AutoTerminate', 'Value': 'true'}
                    ]
                }
            ]
        }
    )

    # Get instance ID from spot request
    spot_request_id = response['SpotInstanceRequests'][0]['SpotInstanceRequestId']

    # Wait for spot request fulfillment
    waiter = ec2.get_waiter('spot_instance_request_fulfilled')
    waiter.wait(SpotInstanceRequestIds=[spot_request_id])

    # Get instance ID
    spot_response = ec2.describe_spot_instance_requests(
        SpotInstanceRequestIds=[spot_request_id]
    )

    instance_id = spot_response['SpotInstanceRequests'][0]['InstanceId']

    return instance_id

def execute_basecalling_command(instance_id: str, run_id: str, bucket: str,
                               input_prefix: str, output_prefix: str,
                               skip_duplex: bool) -> Dict:
    """Execute basecalling script on EC2 instance."""

    command = f"""
#!/bin/bash
set -euo pipefail

# Set variables
export RUN_ID='{run_id}'
export S3_INPUT='s3://{bucket}/{input_prefix}'
export S3_OUTPUT='s3://{bucket}/{output_prefix}'
export SKIP_DUPLEX='{str(skip_duplex).lower()}'

# Create working directory
mkdir -p /mnt/analysis/$RUN_ID
cd /mnt/analysis/$RUN_ID

# Run basecalling script
/opt/minion/scripts/phase1_basecalling/basecall_duplex.sh \\
    -i "$S3_INPUT" \\
    -o "$S3_OUTPUT" \\
    -r "$RUN_ID" \\
    {'-s' if skip_duplex else ''}

# Upload results
aws s3 sync /mnt/analysis/$RUN_ID/output/ "$S3_OUTPUT" \\
    --exclude "*.fast5" \\
    --exclude "*.pod5"

# Clean up local files
rm -rf /mnt/analysis/$RUN_ID

# Signal completion
aws sns publish \\
    --topic-arn {SNS_TOPIC} \\
    --subject "Basecalling Complete - Run $RUN_ID" \\
    --message "Basecalling phase completed for run $RUN_ID"
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={
            'commands': [command],
            'executionTimeout': ['43200']  # 12 hours
        }
    )

    return response

def send_notification(subject: str, message: str):
    """Send SNS notification."""

    sns.publish(
        TopicArn=SNS_TOPIC,
        Subject=subject,
        Message=message
    )