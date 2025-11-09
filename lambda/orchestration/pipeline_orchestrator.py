#!/usr/bin/env python3
"""
Main pipeline orchestrator Lambda function
"""

import json
import boto3
import os
from datetime import datetime
from typing import Dict, Any

s3 = boto3.client('s3')
ec2 = boto3.client('ec2')
stepfunctions = boto3.client('stepfunctions')
sns = boto3.client('sns')
rds = boto3.client('rds-data')

STATE_MACHINE_ARN = os.environ['STATE_MACHINE_ARN']
CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']
SNS_TOPIC = os.environ['SNS_TOPIC_ARN']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Main orchestrator for MinION analysis pipeline.

    Triggered by S3 upload of MinION data.
    """

    # Extract S3 event details
    bucket = event['Records'][0]['s3']['bucket']['name']
    key = event['Records'][0]['s3']['object']['key']

    # Parse run information from S3 key
    # Expected format: runs/{run_id}/fast5/{file}.fast5
    parts = key.split('/')
    if len(parts) < 3:
        return {
            'statusCode': 400,
            'body': 'Invalid S3 key format'
        }

    run_id = parts[1]

    # Create workflow execution record
    workflow_id = create_workflow_execution(run_id, bucket, key)

    # Prepare Step Functions input
    execution_input = {
        'run_id': run_id,
        'workflow_id': workflow_id,
        'bucket': bucket,
        'input_prefix': f'runs/{run_id}/fast5/',
        'output_prefix': f'runs/{run_id}/analysis/',
        'phases': {
            'basecalling': {'enabled': True, 'skip_duplex': False},
            'qc': {'enabled': True, 'min_quality': 9},
            'host_removal': {'enabled': True, 'reference': 'sus_scrofa_11.1'},
            'pathogen_detection': {'enabled': True, 'databases': ['kraken2', 'rvdb', 'pmda']},
            'quantification': {'enabled': True, 'spike_in': 'PhiX174'},
            'reporting': {'enabled': True, 'formats': ['pdf', 'json']}
        }
    }

    # Start Step Functions execution
    response = stepfunctions.start_execution(
        stateMachineArn=STATE_MACHINE_ARN,
        name=f'{run_id}-{datetime.now().strftime("%Y%m%d%H%M%S")}',
        input=json.dumps(execution_input)
    )

    # Update workflow status
    update_workflow_status(workflow_id, 'RUNNING', response['executionArn'])

    # Send notification
    send_notification(
        f'Pipeline Started for Run {run_id}',
        f'Analysis pipeline has been initiated for run {run_id}.\n'
        f'Workflow ID: {workflow_id}\n'
        f'Input: s3://{bucket}/{key}'
    )

    return {
        'statusCode': 200,
        'body': json.dumps({
            'workflow_id': workflow_id,
            'execution_arn': response['executionArn'],
            'status': 'STARTED'
        })
    }

def create_workflow_execution(run_id: str, bucket: str, key: str) -> str:
    """Create workflow execution record in database."""

    sql = """
        INSERT INTO workflow_executions
        (run_id, status, input_path, created_at)
        VALUES (:run_id, :status, :input_path, :created_at)
        RETURNING id
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'run_id', 'value': {'stringValue': run_id}},
            {'name': 'status', 'value': {'stringValue': 'INITIATED'}},
            {'name': 'input_path', 'value': {'stringValue': f's3://{bucket}/{key}'}},
            {'name': 'created_at', 'value': {'stringValue': datetime.now().isoformat()}}
        ]
    )

    return response['records'][0][0]['longValue']

def update_workflow_status(workflow_id: str, status: str, execution_arn: str):
    """Update workflow execution status."""

    sql = """
        UPDATE workflow_executions
        SET status = :status,
            execution_arn = :execution_arn,
            updated_at = :updated_at
        WHERE id = :workflow_id
    """

    rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': int(workflow_id)}},
            {'name': 'status', 'value': {'stringValue': status}},
            {'name': 'execution_arn', 'value': {'stringValue': execution_arn}},
            {'name': 'updated_at', 'value': {'stringValue': datetime.now().isoformat()}}
        ]
    )

def send_notification(subject: str, message: str):
    """Send SNS notification."""

    sns.publish(
        TopicArn=SNS_TOPIC,
        Subject=subject,
        Message=message
    )