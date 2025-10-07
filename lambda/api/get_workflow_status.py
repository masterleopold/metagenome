#!/usr/bin/env python3
"""
API endpoint to get workflow status
"""

import json
import boto3
import os
from typing import Dict, Any

rds = boto3.client('rds-data')
s3 = boto3.client('s3')

CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']
BUCKET = f'minion-analysis-{os.environ["ENVIRONMENT"]}'

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Get workflow execution status.

    API Gateway endpoint: GET /workflows/{workflow_id}
    """

    # Get workflow ID from path parameters
    workflow_id = event.get('pathParameters', {}).get('workflow_id')

    if not workflow_id:
        return {
            'statusCode': 400,
            'headers': {'Content-Type': 'application/json'},
            'body': json.dumps({'error': 'Missing workflow_id'})
        }

    try:
        # Get workflow status from database
        workflow = get_workflow_status(workflow_id)

        if not workflow:
            return {
                'statusCode': 404,
                'headers': {'Content-Type': 'application/json'},
                'body': json.dumps({'error': 'Workflow not found'})
            }

        # Get phase statuses
        phases = get_phase_statuses(workflow_id)

        # Get latest metrics
        metrics = get_workflow_metrics(workflow_id)

        # Get pathogen detection results if available
        pathogens = get_pathogen_results(workflow_id)

        # Build response
        response_data = {
            'workflow_id': workflow_id,
            'run_id': workflow['run_id'],
            'status': workflow['status'],
            'created_at': workflow['created_at'],
            'updated_at': workflow['updated_at'],
            'phases': phases,
            'metrics': metrics,
            'pathogens': pathogens
        }

        # Add report URLs if completed
        if workflow['status'] == 'COMPLETED':
            reports = get_report_urls(workflow['run_id'])
            response_data['reports'] = reports

        return {
            'statusCode': 200,
            'headers': {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*'
            },
            'body': json.dumps(response_data, default=str)
        }

    except Exception as e:
        print(f'Error getting workflow status: {e}')
        return {
            'statusCode': 500,
            'headers': {'Content-Type': 'application/json'},
            'body': json.dumps({'error': 'Internal server error'})
        }

def get_workflow_status(workflow_id: str) -> Dict:
    """Get workflow status from database."""

    sql = """
        SELECT run_id, status, created_at, updated_at, execution_arn
        FROM workflow_executions
        WHERE id = :workflow_id
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': int(workflow_id)}}
        ]
    )

    if response['records']:
        record = response['records'][0]
        return {
            'run_id': record[0]['stringValue'],
            'status': record[1]['stringValue'],
            'created_at': record[2]['stringValue'],
            'updated_at': record[3]['stringValue'],
            'execution_arn': record[4].get('stringValue', '')
        }

    return None

def get_phase_statuses(workflow_id: str) -> list:
    """Get phase execution statuses."""

    sql = """
        SELECT phase_name, status, created_at, updated_at
        FROM phase_executions
        WHERE workflow_id = :workflow_id
        ORDER BY created_at
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': int(workflow_id)}}
        ]
    )

    phases = []
    for record in response['records']:
        phases.append({
            'name': record[0]['stringValue'],
            'status': record[1]['stringValue'],
            'started_at': record[2]['stringValue'],
            'completed_at': record[3].get('stringValue', '')
        })

    return phases

def get_workflow_metrics(workflow_id: str) -> Dict:
    """Get aggregated workflow metrics."""

    sql = """
        SELECT phase_name, metrics
        FROM phase_metrics
        WHERE workflow_id = :workflow_id
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': int(workflow_id)}}
        ]
    )

    metrics = {}
    for record in response['records']:
        phase = record[0]['stringValue']
        phase_metrics = json.loads(record[1]['stringValue'])
        metrics[phase] = phase_metrics

    return metrics

def get_pathogen_results(workflow_id: str) -> list:
    """Get pathogen detection results."""

    sql = """
        SELECT pathogen_code, pathogen_name, read_count, confidence_score
        FROM pathogen_detections
        WHERE workflow_id = :workflow_id
        ORDER BY read_count DESC
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': int(workflow_id)}}
        ]
    )

    pathogens = []
    for record in response['records']:
        pathogens.append({
            'code': record[0]['stringValue'],
            'name': record[1]['stringValue'],
            'read_count': record[2]['longValue'],
            'confidence': record[3]['doubleValue']
        })

    return pathogens

def get_report_urls(run_id: str) -> Dict:
    """Generate presigned URLs for reports."""

    reports = {}
    report_files = [
        ('pdf', f'{run_id}_report.pdf'),
        ('html', f'{run_id}_report.html'),
        ('json', f'{run_id}_pmda_report.json'),
        ('checklist', 'pmda_checklist.json')
    ]

    for report_type, filename in report_files:
        key = f'runs/{run_id}/analysis/reports/{filename}'

        # Check if file exists
        try:
            s3.head_object(Bucket=BUCKET, Key=key)

            # Generate presigned URL (valid for 24 hours)
            url = s3.generate_presigned_url(
                'get_object',
                Params={'Bucket': BUCKET, 'Key': key},
                ExpiresIn=86400
            )
            reports[report_type] = url

        except:
            pass  # File doesn't exist

    return reports