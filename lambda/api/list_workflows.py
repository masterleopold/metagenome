#!/usr/bin/env python3
"""
API endpoint to list workflow executions
"""

import json
import boto3
import os
from typing import Dict, Any
from datetime import datetime, timedelta

rds = boto3.client('rds-data')

CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    List workflow executions with filtering.

    API Gateway endpoint: GET /workflows
    Query parameters:
    - status: Filter by status (RUNNING, COMPLETED, FAILED)
    - run_id: Filter by run ID
    - limit: Maximum number of results (default 20, max 100)
    - offset: Pagination offset
    - days: Get workflows from last N days
    """

    # Parse query parameters
    params = event.get('queryStringParameters', {}) or {}
    status_filter = params.get('status')
    run_id_filter = params.get('run_id')
    limit = min(int(params.get('limit', 20)), 100)
    offset = int(params.get('offset', 0))
    days = int(params.get('days', 7))

    try:
        # Build SQL query
        conditions = []
        sql_params = []

        # Add time filter
        cutoff_date = datetime.utcnow() - timedelta(days=days)
        conditions.append("created_at >= :cutoff_date")
        sql_params.append({
            'name': 'cutoff_date',
            'value': {'stringValue': cutoff_date.isoformat()}
        })

        # Add status filter
        if status_filter:
            conditions.append("status = :status")
            sql_params.append({
                'name': 'status',
                'value': {'stringValue': status_filter}
            })

        # Add run_id filter
        if run_id_filter:
            conditions.append("run_id LIKE :run_id")
            sql_params.append({
                'name': 'run_id',
                'value': {'stringValue': f'%{run_id_filter}%'}
            })

        where_clause = " AND ".join(conditions) if conditions else "1=1"

        # Get total count
        count_sql = f"""
            SELECT COUNT(*) FROM workflow_executions
            WHERE {where_clause}
        """

        count_response = rds.execute_statement(
            resourceArn=CLUSTER_ARN,
            secretArn=SECRET_ARN,
            database=DATABASE,
            sql=count_sql,
            parameters=sql_params
        )

        total_count = count_response['records'][0][0]['longValue'] if count_response['records'] else 0

        # Get workflows
        list_sql = f"""
            SELECT
                w.id, w.run_id, w.status, w.created_at, w.updated_at,
                COUNT(DISTINCT p.phase_name) as phases_completed,
                COUNT(DISTINCT e.id) as error_count,
                MAX(pd.pathogen_code IS NOT NULL) as has_pathogens
            FROM workflow_executions w
            LEFT JOIN phase_executions p ON w.id = p.workflow_id AND p.status = 'COMPLETED'
            LEFT JOIN phase_errors e ON w.id = e.workflow_id
            LEFT JOIN pathogen_detections pd ON w.id = pd.workflow_id
            WHERE {where_clause}
            GROUP BY w.id, w.run_id, w.status, w.created_at, w.updated_at
            ORDER BY w.created_at DESC
            LIMIT :limit OFFSET :offset
        """

        sql_params.extend([
            {'name': 'limit', 'value': {'longValue': limit}},
            {'name': 'offset', 'value': {'longValue': offset}}
        ])

        response = rds.execute_statement(
            resourceArn=CLUSTER_ARN,
            secretArn=SECRET_ARN,
            database=DATABASE,
            sql=list_sql,
            parameters=sql_params
        )

        # Format workflows
        workflows = []
        for record in response.get('records', []):
            workflow = {
                'workflow_id': record[0]['longValue'],
                'run_id': record[1]['stringValue'],
                'status': record[2]['stringValue'],
                'created_at': record[3]['stringValue'],
                'updated_at': record[4].get('stringValue', ''),
                'phases_completed': record[5]['longValue'],
                'error_count': record[6]['longValue'],
                'has_pathogens': bool(record[7].get('booleanValue', False))
            }

            # Calculate duration
            if workflow['updated_at']:
                start = datetime.fromisoformat(workflow['created_at'])
                end = datetime.fromisoformat(workflow['updated_at'])
                duration = (end - start).total_seconds()
                workflow['duration_seconds'] = int(duration)

            workflows.append(workflow)

        # Build response
        response_data = {
            'workflows': workflows,
            'pagination': {
                'total': total_count,
                'limit': limit,
                'offset': offset,
                'has_more': (offset + limit) < total_count
            },
            'filters': {
                'status': status_filter,
                'run_id': run_id_filter,
                'days': days
            }
        }

        return {
            'statusCode': 200,
            'headers': {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*'
            },
            'body': json.dumps(response_data, default=str)
        }

    except Exception as e:
        print(f'Error listing workflows: {e}')
        return {
            'statusCode': 500,
            'headers': {'Content-Type': 'application/json'},
            'body': json.dumps({'error': 'Internal server error'})
        }