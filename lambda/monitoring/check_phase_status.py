#!/usr/bin/env python3
"""
Check phase execution status
"""

import json
import boto3
import os
from typing import Dict, Any

ssm = boto3.client('ssm')
ec2 = boto3.client('ec2')
rds = boto3.client('rds-data')

CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Check the status of a running phase.

    Called by Step Functions to monitor phase execution.
    """

    phase = event.get('phase')
    instance_id = event.get('instance_id')
    command_id = event.get('command_id')
    run_id = event.get('run_id')
    workflow_id = event.get('workflow_id')

    # Check SSM command status
    command_status = check_command_status(instance_id, command_id)

    # Check instance status
    instance_status = check_instance_status(instance_id)

    # Determine phase status
    if command_status == 'Success':
        phase_status = 'COMPLETED'
        update_phase_status(workflow_id, phase, 'COMPLETED')
    elif command_status == 'Failed':
        phase_status = 'FAILED'
        update_phase_status(workflow_id, phase, 'FAILED')
        # Get error details
        error_details = get_command_error(instance_id, command_id)
        log_phase_error(workflow_id, phase, error_details)
    elif command_status in ['InProgress', 'Pending']:
        phase_status = 'RUNNING'
    else:
        phase_status = 'UNKNOWN'

    # Check if instance is still running
    if instance_status != 'running' and phase_status == 'RUNNING':
        # Instance terminated while command was running
        phase_status = 'FAILED'
        error_details = {'message': f'Instance {instance_id} terminated unexpectedly'}
        log_phase_error(workflow_id, phase, error_details)

    result = {
        'phase': phase,
        'status': phase_status,
        'instance_id': instance_id,
        'instance_status': instance_status,
        'command_id': command_id,
        'command_status': command_status,
        'run_id': run_id,
        'workflow_id': workflow_id
    }

    # Add error details if failed
    if phase_status == 'FAILED':
        result['error'] = error_details

    return result

def check_command_status(instance_id: str, command_id: str) -> str:
    """Check SSM command execution status."""

    try:
        response = ssm.get_command_invocation(
            CommandId=command_id,
            InstanceId=instance_id
        )
        return response['Status']
    except Exception as e:
        print(f'Error checking command status: {e}')
        return 'Unknown'

def check_instance_status(instance_id: str) -> str:
    """Check EC2 instance status."""

    try:
        response = ec2.describe_instances(InstanceIds=[instance_id])
        if response['Reservations']:
            instance = response['Reservations'][0]['Instances'][0]
            return instance['State']['Name']
    except Exception as e:
        print(f'Error checking instance status: {e}')

    return 'unknown'

def get_command_error(instance_id: str, command_id: str) -> Dict:
    """Get command error details."""

    try:
        response = ssm.get_command_invocation(
            CommandId=command_id,
            InstanceId=instance_id
        )

        error = {
            'message': response.get('StatusDetails', 'Command failed'),
            'stdout': response.get('StandardOutputContent', '')[-1000:],  # Last 1000 chars
            'stderr': response.get('StandardErrorContent', '')[-1000:]
        }

        return error
    except Exception as e:
        return {'message': str(e)}

def update_phase_status(workflow_id: str, phase: str, status: str):
    """Update phase execution status in database."""

    sql = """
        UPDATE phase_executions
        SET status = :status, updated_at = NOW()
        WHERE workflow_id = :workflow_id AND phase_name = :phase_name
    """

    rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': workflow_id}},
            {'name': 'phase_name', 'value': {'stringValue': phase}},
            {'name': 'status', 'value': {'stringValue': status}}
        ]
    )

def log_phase_error(workflow_id: str, phase: str, error: Dict):
    """Log phase execution error to database."""

    sql = """
        INSERT INTO phase_errors
        (workflow_id, phase_name, error_message, error_details, created_at)
        VALUES (:workflow_id, :phase_name, :error_message, :error_details, NOW())
    """

    rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': workflow_id}},
            {'name': 'phase_name', 'value': {'stringValue': phase}},
            {'name': 'error_message', 'value': {'stringValue': error.get('message', 'Unknown error')}},
            {'name': 'error_details', 'value': {'stringValue': json.dumps(error)}}
        ]
    )