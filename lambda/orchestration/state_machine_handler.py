#!/usr/bin/env python3
"""
Step Functions state machine handler
"""

import json
import boto3
import os
from typing import Dict, Any

ec2 = boto3.client('ec2')
s3 = boto3.client('s3')
rds = boto3.client('rds-data')

CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Handle state transitions in Step Functions.

    Manages phase transitions and error handling.
    """

    phase = event.get('phase')
    action = event.get('action')
    run_id = event.get('run_id')
    workflow_id = event.get('workflow_id')

    if action == 'check_phase_completion':
        return check_phase_completion(phase, run_id, workflow_id)

    elif action == 'handle_phase_error':
        return handle_phase_error(phase, run_id, workflow_id, event.get('error'))

    elif action == 'prepare_next_phase':
        return prepare_next_phase(phase, run_id, workflow_id, event)

    elif action == 'finalize_workflow':
        return finalize_workflow(run_id, workflow_id, event)

    else:
        raise ValueError(f'Unknown action: {action}')

def check_phase_completion(phase: str, run_id: str, workflow_id: str) -> Dict[str, Any]:
    """Check if a phase has completed successfully."""

    # Check output files in S3
    bucket = f'minion-analysis-{os.environ["ENVIRONMENT"]}'
    prefix = f'runs/{run_id}/analysis/{phase}/'

    response = s3.list_objects_v2(
        Bucket=bucket,
        Prefix=prefix,
        MaxKeys=1
    )

    if 'Contents' not in response:
        return {
            'completed': False,
            'phase': phase,
            'message': f'No output files found for {phase}'
        }

    # Update phase status in database
    update_phase_status(workflow_id, phase, 'COMPLETED')

    return {
        'completed': True,
        'phase': phase,
        'output_path': f's3://{bucket}/{prefix}'
    }

def handle_phase_error(phase: str, run_id: str, workflow_id: str, error: Dict) -> Dict[str, Any]:
    """Handle phase execution errors."""

    # Log error to database
    log_phase_error(workflow_id, phase, error)

    # Determine if error is recoverable
    recoverable = is_recoverable_error(error)

    if recoverable:
        # Retry logic
        return {
            'action': 'RETRY',
            'phase': phase,
            'retry_count': error.get('retry_count', 0) + 1,
            'max_retries': 3
        }
    else:
        # Fail workflow
        update_workflow_status(workflow_id, 'FAILED')
        return {
            'action': 'FAIL',
            'phase': phase,
            'error': error
        }

def prepare_next_phase(current_phase: str, run_id: str, workflow_id: str, event: Dict) -> Dict[str, Any]:
    """Prepare configuration for next phase."""

    phase_order = [
        'basecalling', 'qc', 'host_removal',
        'pathogen_detection', 'quantification', 'reporting'
    ]

    current_index = phase_order.index(current_phase)

    if current_index >= len(phase_order) - 1:
        # All phases complete
        return {
            'next_phase': None,
            'complete': True
        }

    next_phase = phase_order[current_index + 1]

    # Prepare phase-specific configuration
    phase_config = event['phases'].get(next_phase, {})

    # Get outputs from previous phase
    bucket = f'minion-analysis-{os.environ["ENVIRONMENT"]}'
    prev_output = f's3://{bucket}/runs/{run_id}/analysis/{current_phase}/'

    return {
        'next_phase': next_phase,
        'complete': False,
        'config': phase_config,
        'input_path': prev_output
    }

def finalize_workflow(run_id: str, workflow_id: str, event: Dict) -> Dict[str, Any]:
    """Finalize workflow execution."""

    # Update workflow status
    update_workflow_status(workflow_id, 'COMPLETED')

    # Collect all outputs
    outputs = collect_workflow_outputs(run_id, workflow_id)

    # Generate final summary
    summary = generate_workflow_summary(run_id, workflow_id, outputs)

    # Store summary
    bucket = f'minion-analysis-{os.environ["ENVIRONMENT"]}'
    key = f'runs/{run_id}/analysis/workflow_summary.json'

    s3.put_object(
        Bucket=bucket,
        Key=key,
        Body=json.dumps(summary, indent=2),
        ContentType='application/json'
    )

    return {
        'status': 'COMPLETED',
        'run_id': run_id,
        'workflow_id': workflow_id,
        'summary_path': f's3://{bucket}/{key}',
        'outputs': outputs
    }

def update_phase_status(workflow_id: str, phase: str, status: str):
    """Update phase execution status in database."""

    sql = """
        INSERT INTO phase_executions (workflow_id, phase_name, status)
        VALUES (:workflow_id, :phase_name, :status)
        ON CONFLICT (workflow_id, phase_name)
        DO UPDATE SET status = :status, updated_at = NOW()
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

def update_workflow_status(workflow_id: str, status: str):
    """Update workflow execution status."""

    sql = """
        UPDATE workflow_executions
        SET status = :status, updated_at = NOW()
        WHERE id = :workflow_id
    """

    rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': workflow_id}},
            {'name': 'status', 'value': {'stringValue': status}}
        ]
    )

def log_phase_error(workflow_id: str, phase: str, error: Dict):
    """Log phase execution error."""

    sql = """
        INSERT INTO phase_errors
        (workflow_id, phase_name, error_message, error_details)
        VALUES (:workflow_id, :phase_name, :error_message, :error_details)
    """

    rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': workflow_id}},
            {'name': 'phase_name', 'value': {'stringValue': phase}},
            {'name': 'error_message', 'value': {'stringValue': str(error.get('message', 'Unknown error'))}},
            {'name': 'error_details', 'value': {'stringValue': json.dumps(error)}}
        ]
    )

def is_recoverable_error(error: Dict) -> bool:
    """Determine if error is recoverable."""

    recoverable_errors = [
        'EC2 capacity error',
        'S3 timeout',
        'Network error'
    ]

    error_message = str(error.get('message', '')).lower()

    for recoverable in recoverable_errors:
        if recoverable.lower() in error_message:
            return True

    return False

def collect_workflow_outputs(run_id: str, workflow_id: str) -> Dict[str, str]:
    """Collect all workflow output paths."""

    bucket = f'minion-analysis-{os.environ["ENVIRONMENT"]}'
    phases = ['basecalling', 'qc', 'host_removal', 'pathogen_detection', 'quantification', 'reporting']

    outputs = {}
    for phase in phases:
        outputs[phase] = f's3://{bucket}/runs/{run_id}/analysis/{phase}/'

    return outputs

def generate_workflow_summary(run_id: str, workflow_id: str, outputs: Dict) -> Dict:
    """Generate workflow execution summary."""

    # Get workflow metrics from database
    sql = """
        SELECT
            w.status, w.created_at, w.updated_at,
            COUNT(DISTINCT p.phase_name) as phases_completed,
            COUNT(DISTINCT e.id) as errors_encountered
        FROM workflow_executions w
        LEFT JOIN phase_executions p ON w.id = p.workflow_id AND p.status = 'COMPLETED'
        LEFT JOIN phase_errors e ON w.id = e.workflow_id
        WHERE w.id = :workflow_id
        GROUP BY w.id, w.status, w.created_at, w.updated_at
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': workflow_id}}
        ]
    )

    if response['records']:
        record = response['records'][0]
        summary = {
            'run_id': run_id,
            'workflow_id': workflow_id,
            'status': record[0]['stringValue'],
            'created_at': record[1]['stringValue'],
            'updated_at': record[2]['stringValue'],
            'phases_completed': record[3]['longValue'],
            'errors_encountered': record[4]['longValue'],
            'outputs': outputs
        }
    else:
        summary = {
            'run_id': run_id,
            'workflow_id': workflow_id,
            'outputs': outputs
        }

    return summary