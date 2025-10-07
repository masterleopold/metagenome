#!/usr/bin/env python3
"""
API endpoint to get pathogen detection results
"""

import json
import boto3
import os
from typing import Dict, Any, List

rds = boto3.client('rds-data')
s3 = boto3.client('s3')

CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']
BUCKET = f'minion-analysis-{os.environ["ENVIRONMENT"]}'

# PMDA 91 pathogen codes
PMDA_PATHOGENS = set([
    'PERV-A', 'PERV-B', 'PERV-C', 'HEV', 'JEV', 'PRRSV', 'PCV2', 'PRV',
    'FMDV', 'ASFV', 'CSFV', 'SIV', 'PPV', 'EMCV', 'RV', 'PEDV', 'TGEV',
    'SA', 'SP', 'SS', 'EC', 'SE', 'CT', 'CP', 'LA', 'BA', 'BP', 'FT',
    'YP', 'MT', 'MB', 'MA', 'CJ', 'HP', 'BB', 'LP', 'TP', 'CS', 'PM',
    'AP', 'HPS', 'EP', 'BR', 'CN', 'CA', 'AF', 'PJ', 'TG', 'TC', 'TS',
    'EC-P', 'CP-P', 'GD', 'SS-P', 'PRION'
])

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Get pathogen detection results for a run.

    API Gateway endpoint: GET /pathogens/{run_id}
    Query parameters:
    - pmda_only: Show only PMDA 91 pathogens
    - min_reads: Minimum read count threshold
    - include_quantification: Include absolute quantification data
    """

    # Get run ID from path
    run_id = event.get('pathParameters', {}).get('run_id')

    if not run_id:
        return {
            'statusCode': 400,
            'headers': {'Content-Type': 'application/json'},
            'body': json.dumps({'error': 'Missing run_id'})
        }

    # Parse query parameters
    params = event.get('queryStringParameters', {}) or {}
    pmda_only = params.get('pmda_only', 'false').lower() == 'true'
    min_reads = int(params.get('min_reads', 0))
    include_quant = params.get('include_quantification', 'false').lower() == 'true'

    try:
        # Get workflow ID for run
        workflow_id = get_workflow_id(run_id)

        if not workflow_id:
            return {
                'statusCode': 404,
                'headers': {'Content-Type': 'application/json'},
                'body': json.dumps({'error': 'Run not found'})
            }

        # Get pathogen detections
        pathogens = get_pathogen_detections(workflow_id, pmda_only, min_reads)

        # Get quantification data if requested
        if include_quant:
            quantification = get_quantification_data(run_id)
            # Merge quantification with detection data
            for pathogen in pathogens:
                code = pathogen['code']
                if code in quantification:
                    pathogen['quantification'] = quantification[code]

        # Get PERV analysis results
        perv_results = get_perv_analysis(workflow_id)

        # Build PMDA compliance summary
        pmda_summary = build_pmda_summary(pathogens)

        response_data = {
            'run_id': run_id,
            'workflow_id': workflow_id,
            'total_pathogens': len(pathogens),
            'pmda_pathogens_detected': pmda_summary['detected'],
            'pmda_compliance': pmda_summary,
            'perv_analysis': perv_results,
            'pathogens': pathogens
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
        print(f'Error getting pathogen results: {e}')
        return {
            'statusCode': 500,
            'headers': {'Content-Type': 'application/json'},
            'body': json.dumps({'error': 'Internal server error'})
        }

def get_workflow_id(run_id: str) -> int:
    """Get workflow ID for a run."""

    sql = """
        SELECT id FROM workflow_executions
        WHERE run_id = :run_id
        ORDER BY created_at DESC
        LIMIT 1
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'run_id', 'value': {'stringValue': run_id}}
        ]
    )

    if response['records']:
        return response['records'][0][0]['longValue']

    return None

def get_pathogen_detections(workflow_id: int, pmda_only: bool, min_reads: int) -> List[Dict]:
    """Get pathogen detection results."""

    conditions = ["workflow_id = :workflow_id"]
    params = [{'name': 'workflow_id', 'value': {'longValue': workflow_id}}]

    if pmda_only:
        pmda_list = "'" + "','".join(PMDA_PATHOGENS) + "'"
        conditions.append(f"pathogen_code IN ({pmda_list})")

    if min_reads > 0:
        conditions.append("read_count >= :min_reads")
        params.append({'name': 'min_reads', 'value': {'longValue': min_reads}})

    where_clause = " AND ".join(conditions)

    sql = f"""
        SELECT
            pathogen_code, pathogen_name, read_count,
            confidence_score, detection_method
        FROM pathogen_detections
        WHERE {where_clause}
        ORDER BY read_count DESC
    """

    response = rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=params
    )

    pathogens = []
    for record in response.get('records', []):
        pathogen = {
            'code': record[0]['stringValue'],
            'name': record[1]['stringValue'],
            'read_count': record[2]['longValue'],
            'confidence': record[3]['doubleValue'],
            'method': record[4]['stringValue'],
            'is_pmda': record[0]['stringValue'] in PMDA_PATHOGENS
        }

        # Classify risk level
        if pathogen['code'].startswith('PERV'):
            pathogen['risk_level'] = 'CRITICAL'
        elif pathogen['code'] in ['ASFV', 'CSFV', 'FMDV']:
            pathogen['risk_level'] = 'HIGH'
        elif pathogen['is_pmda']:
            pathogen['risk_level'] = 'MEDIUM'
        else:
            pathogen['risk_level'] = 'LOW'

        pathogens.append(pathogen)

    return pathogens

def get_quantification_data(run_id: str) -> Dict:
    """Get absolute quantification data from S3."""

    quantification = {}

    try:
        # Read absolute quantification file
        key = f'runs/{run_id}/analysis/quantification/absolute_quantification.json'
        response = s3.get_object(Bucket=BUCKET, Key=key)
        data = json.loads(response['Body'].read())

        for code, quant in data.get('pathogens', {}).items():
            quantification[code] = {
                'copies_per_ml': quant.get('copies_per_ml', 0),
                'log10_copies': quant.get('log10_copies_per_ml', 0),
                'ci_lower': quant.get('ci_lower', 0),
                'ci_upper': quant.get('ci_upper', 0)
            }

    except Exception as e:
        print(f'Error reading quantification data: {e}')

    return quantification

def get_perv_analysis(workflow_id: int) -> Dict:
    """Get PERV-specific analysis results."""

    sql = """
        SELECT
            perv_subtype, integration_sites, copy_number,
            recombination_detected, phylogenetic_cluster
        FROM perv_analysis
        WHERE workflow_id = :workflow_id
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
        return {
            'detected': True,
            'subtype': record[0]['stringValue'],
            'integration_sites': record[1]['longValue'],
            'copy_number': record[2]['doubleValue'],
            'recombination': record[3].get('booleanValue', False),
            'cluster': record[4].get('stringValue', '')
        }

    return {'detected': False}

def build_pmda_summary(pathogens: List[Dict]) -> Dict:
    """Build PMDA compliance summary."""

    detected_pmda = set()
    for pathogen in pathogens:
        if pathogen['is_pmda']:
            detected_pmda.add(pathogen['code'])

    return {
        'total': len(PMDA_PATHOGENS),
        'detected': len(detected_pmda),
        'not_detected': len(PMDA_PATHOGENS - detected_pmda),
        'detection_rate': len(detected_pmda) / len(PMDA_PATHOGENS),
        'detected_list': list(detected_pmda),
        'critical_detected': any(p.startswith('PERV') for p in detected_pmda)
    }