#!/usr/bin/env python3
"""
IMPROVED Phase 4 Pathogen Detection Lambda Handler (v2).

Demonstrates:
- AWS Lambda Powertools structured logging
- Pydantic type safety
- Repository pattern for database access
- Proper error handling with audit logs
- CloudWatch Metrics integration

This is an example of refactored code following expert recommendations.
"""

import json
import os
from typing import Any, Dict

import boto3
from aws_lambda_powertools import Logger, Tracer, Metrics
from aws_lambda_powertools.logging import correlation_paths
from aws_lambda_powertools.metrics import MetricUnit
from botocore.exceptions import ClientError

# Import new type-safe models
from lib.models.workflow import WorkflowStatus, WorkflowPhase
from lib.models.database import WorkflowRecord
from lib.repositories.rds_repository import RDSWorkflowRepository
from lib.logging.logger import AuditLogger

# Initialize Lambda Powertools (singleton pattern)
logger = Logger(service="pathogen-detection")
tracer = Tracer(service="pathogen-detection")
metrics = Metrics(namespace="MinION/Pipeline", service="pathogen-detection")
audit_logger = AuditLogger(service="pathogen-detection")

# AWS clients
ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')
sns = boto3.client('sns')

# Environment variables
ANALYSIS_AMI = os.environ['ANALYSIS_AMI_ID']
INSTANCE_TYPE = os.environ.get('PATHOGEN_INSTANCE_TYPE', 'r5.4xlarge')
SUBNET_ID = os.environ['SUBNET_ID']
SECURITY_GROUP = os.environ['SECURITY_GROUP_ID']
EFS_ID = os.environ['EFS_ID']
SNS_TOPIC = os.environ['SNS_TOPIC_ARN']

# Database repository (initialized lazily)
_workflow_repo = None


def get_workflow_repository() -> RDSWorkflowRepository:
    """Get or create workflow repository instance."""
    global _workflow_repo
    if _workflow_repo is None:
        _workflow_repo = RDSWorkflowRepository(
            cluster_arn=os.environ['RDS_CLUSTER_ARN'],
            secret_arn=os.environ['RDS_SECRET_ARN'],
            database=os.environ.get('RDS_DATABASE', 'minion_metadata'),
        )
    return _workflow_repo


@logger.inject_lambda_context(correlation_id_path=correlation_paths.EVENT_BRIDGE)
@tracer.capture_lambda_handler
@metrics.log_metrics(capture_cold_start_metric=True)
def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Lambda handler for Phase 4: Multi-database pathogen detection.

    IMPROVED VERSION with:
    - Structured logging (CloudWatch Logs Insights compatible)
    - Automatic correlation IDs for request tracing
    - X-Ray distributed tracing
    - CloudWatch Metrics for monitoring
    - Audit logging for PMDA compliance

    Args:
        event: Lambda event dict containing:
            - run_id (str): Unique run identifier
            - workflow_id (int): Step Functions execution ID
            - input_path (str): S3 path to host-depleted FASTQ files
            - operator_email (str): Operator email for audit trail
            - config (dict, optional): Configuration

        context: Lambda context object

    Returns:
        Dict with execution details

    Raises:
        ValueError: If required parameters missing or invalid
        RuntimeError: If AWS operations fail
    """

    # Extract and validate inputs
    try:
        run_id = event['run_id']
        workflow_id = event.get('workflow_id', 'unknown')
        input_path = event['input_path']
        operator_email = event.get('operator_email', 'unknown@example.com')
        databases = event.get('config', {}).get('databases', ['kraken2', 'rvdb', 'pmda'])

        logger.info(
            "Phase 4 pathogen detection initiated",
            extra={
                "run_id": run_id,
                "workflow_id": workflow_id,
                "input_path": input_path,
                "databases": databases,
                "instance_type": INSTANCE_TYPE,
            }
        )

    except KeyError as e:
        logger.error(f"Missing required parameter: {e}")
        metrics.add_metric(name="ValidationErrors", unit=MetricUnit.Count, value=1)
        raise ValueError(f"Missing required parameter: {e}")

    # Parse S3 path
    parts = input_path.replace('s3://', '').split('/', 1)
    bucket = parts[0]
    prefix = parts[1] if len(parts) > 1 else ''
    output_prefix = prefix.replace('host_removal', 'pathogen_detection')

    # Update workflow status in database
    try:
        repo = get_workflow_repository()
        repo.update_status(
            run_id=run_id,
            status=WorkflowStatus.PATHOGEN_DETECTION,
            error_message=None
        )

        audit_logger.log_database_operation(
            operation="UPDATE",
            table="workflow_executions",
            run_id=run_id,
            success=True,
            details={"status": "pathogen_detection"}
        )

    except Exception as e:
        logger.warning(
            "Failed to update workflow status in database",
            extra={"error": str(e), "run_id": run_id}
        )
        # Don't fail the entire workflow for DB errors
        metrics.add_metric(name="DatabaseErrors", unit=MetricUnit.Count, value=1)

    # Launch EC2 instance
    try:
        instance_id = launch_pathogen_instance(
            run_id=run_id,
            bucket=bucket,
            input_prefix=prefix,
            output_prefix=output_prefix
        )

        logger.info(
            "Pathogen detection instance launched",
            extra={
                "run_id": run_id,
                "instance_id": instance_id,
                "instance_type": INSTANCE_TYPE,
            }
        )

        metrics.add_metric(name="InstancesLaunched", unit=MetricUnit.Count, value=1)

    except ClientError as e:
        logger.error(
            "Failed to launch EC2 instance",
            extra={
                "error": str(e),
                "error_code": e.response['Error']['Code'],
                "run_id": run_id,
            }
        )

        # Log error to audit trail
        audit_logger.log_error(
            run_id=run_id,
            phase="pathogen_detection",
            error_message=str(e),
            error_type=type(e).__name__,
            operator_email=operator_email,
        )

        metrics.add_metric(name="LaunchErrors", unit=MetricUnit.Count, value=1)
        raise RuntimeError(f"Failed to launch instance: {e}")

    # Wait for instance to be running
    try:
        with tracer.capture_method(name="wait_for_instance"):
            waiter = ec2.get_waiter('instance_running')
            waiter.wait(InstanceIds=[instance_id])

        logger.info(
            "Instance is running",
            extra={"run_id": run_id, "instance_id": instance_id}
        )

    except Exception as e:
        logger.error(
            "Instance failed to start",
            extra={"error": str(e), "instance_id": instance_id}
        )
        raise

    # Execute pathogen detection via SSM
    try:
        command_id = execute_pathogen_detection(
            instance_id=instance_id,
            run_id=run_id,
            bucket=bucket,
            input_prefix=prefix,
            output_prefix=output_prefix,
            databases=databases
        )

        logger.info(
            "Pathogen detection command sent",
            extra={
                "run_id": run_id,
                "command_id": command_id,
                "databases": databases,
            }
        )

    except Exception as e:
        logger.error(
            "Failed to send SSM command",
            extra={"error": str(e), "run_id": run_id}
        )
        raise

    # Prepare response
    response = {
        'phase': 'pathogen_detection',
        'status': 'RUNNING',
        'instance_id': instance_id,
        'command_id': command_id,
        'run_id': run_id,
        'workflow_id': workflow_id,
        'output_path': f's3://{bucket}/{output_prefix}'
    }

    logger.info(
        "Phase 4 pathogen detection started successfully",
        extra=response
    )

    return response


@tracer.capture_method
def launch_pathogen_instance(
    run_id: str,
    bucket: str,
    input_prefix: str,
    output_prefix: str
) -> str:
    """
    Launch high-memory EC2 instance for pathogen detection.

    Uses spot instances by default for cost savings (70% reduction).
    Falls back to on-demand if spot unavailable.

    Args:
        run_id: Unique run identifier
        bucket: S3 bucket name
        input_prefix: S3 input prefix
        output_prefix: S3 output prefix

    Returns:
        Instance ID

    Raises:
        ClientError: If launch fails
    """
    import shlex

    user_data = f"""#!/bin/bash
# Mount EFS for databases
apt-get update && apt-get install -y nfs-common
mkdir -p /mnt/efs
mount -t nfs4 {EFS_ID}.efs.{os.environ['AWS_REGION']}.amazonaws.com:/ /mnt/efs

export RUN_ID={shlex.quote(run_id)}
export S3_INPUT={shlex.quote(f's3://{bucket}/{input_prefix}')}
export S3_OUTPUT={shlex.quote(f's3://{bucket}/{output_prefix}')}

# Optimize for Kraken2 (memory-intensive)
echo 3 > /proc/sys/vm/drop_caches
echo 1 > /proc/sys/vm/compact_memory

echo "Pathogen detection instance started for run $RUN_ID" | logger

# Auto-terminate after 8 hours to prevent runaway costs
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

    # Try spot instance first
    use_spot = os.environ.get('USE_SPOT_INSTANCES', 'true').lower() == 'true'

    if use_spot:
        try:
            logger.debug("Attempting spot instance launch")

            instance_config['InstanceMarketOptions'] = {
                'MarketType': 'spot',
                'SpotOptions': {
                    'MaxPrice': os.environ.get('SPOT_MAX_PRICE', '0.50'),
                    'SpotInstanceType': 'one-time'
                }
            }

            response = ec2.run_instances(MinCount=1, MaxCount=1, **instance_config)
            instance_id = response['Instances'][0]['InstanceId']

            logger.info("Spot instance launched", extra={"instance_id": instance_id})
            metrics.add_metric(name="SpotInstancesLaunched", unit=MetricUnit.Count, value=1)

            return instance_id

        except ClientError as e:
            logger.warning(
                "Spot instance unavailable, falling back to on-demand",
                extra={"error": str(e)}
            )
            # Remove spot config for on-demand fallback
            instance_config.pop('InstanceMarketOptions', None)

    # On-demand instance (or fallback from spot)
    response = ec2.run_instances(MinCount=1, MaxCount=1, **instance_config)
    instance_id = response['Instances'][0]['InstanceId']

    logger.info("On-demand instance launched", extra={"instance_id": instance_id})
    metrics.add_metric(name="OnDemandInstancesLaunched", unit=MetricUnit.Count, value=1)

    return instance_id


@tracer.capture_method
def execute_pathogen_detection(
    instance_id: str,
    run_id: str,
    bucket: str,
    input_prefix: str,
    output_prefix: str,
    databases: list[str]
) -> str:
    """
    Execute pathogen detection script via SSM.

    Args:
        instance_id: EC2 instance ID
        run_id: Run identifier
        bucket: S3 bucket
        input_prefix: Input prefix
        output_prefix: Output prefix
        databases: List of databases to use

    Returns:
        SSM command ID
    """
    import shlex

    # Build safe command
    db_args = ' '.join([f'--db {shlex.quote(db)}' for db in databases])

    command = f"""#!/bin/bash
set -euo pipefail

# Download input
aws s3 sync s3://{shlex.quote(bucket)}/{shlex.quote(input_prefix)} /workspace/input/

# Run pathogen detection (all 91 PMDA pathogens + PERV)
/opt/scripts/phase4_pathogen/detect_pmda_all_91_pathogens.py \\
    --input /workspace/input \\
    --output /workspace/output \\
    --run-id {shlex.quote(run_id)} \\
    {db_args}

# Upload results
aws s3 sync /workspace/output/ s3://{shlex.quote(bucket)}/{shlex.quote(output_prefix)}/

# Terminate instance when done
sudo shutdown -h now
"""

    response = ssm.send_command(
        InstanceIds=[instance_id],
        DocumentName='AWS-RunShellScript',
        Parameters={'commands': [command]},
        CloudWatchOutputConfig={
            'CloudWatchLogGroupName': f'/minion/pathogen-detection/{run_id}',
            'CloudWatchOutputEnabled': True
        }
    )

    return response['Command']['CommandId']
