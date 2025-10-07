#!/usr/bin/env python3
"""
Collect and store analysis metrics
"""

import json
import boto3
import os
from typing import Dict, Any
from datetime import datetime
from decimal import Decimal

s3 = boto3.client('s3')
cloudwatch = boto3.client('cloudwatch')
rds = boto3.client('rds-data')

CLUSTER_ARN = os.environ['CLUSTER_ARN']
SECRET_ARN = os.environ['SECRET_ARN']
DATABASE = os.environ['DATABASE']
BUCKET = f'minion-analysis-{os.environ["ENVIRONMENT"]}'

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """
    Collect metrics from completed analysis phases.

    Triggered after each phase completion.
    """

    phase = event.get('phase')
    run_id = event.get('run_id')
    workflow_id = event.get('workflow_id')
    output_path = event.get('output_path', '')

    # Parse S3 path
    if output_path.startswith('s3://'):
        parts = output_path.replace('s3://', '').split('/', 1)
        bucket = parts[0]
        prefix = parts[1] if len(parts) > 1 else ''
    else:
        bucket = BUCKET
        prefix = f'runs/{run_id}/analysis/{phase}'

    # Collect phase-specific metrics
    metrics = {}

    if phase == 'basecalling':
        metrics = collect_basecalling_metrics(bucket, prefix)
    elif phase == 'qc':
        metrics = collect_qc_metrics(bucket, prefix)
    elif phase == 'host_removal':
        metrics = collect_host_removal_metrics(bucket, prefix)
    elif phase == 'pathogen_detection':
        metrics = collect_pathogen_metrics(bucket, prefix)
    elif phase == 'quantification':
        metrics = collect_quantification_metrics(bucket, prefix)
    elif phase == 'reporting':
        metrics = collect_reporting_metrics(bucket, prefix)

    # Store metrics in database
    store_metrics(workflow_id, phase, metrics)

    # Send metrics to CloudWatch
    send_cloudwatch_metrics(phase, run_id, metrics)

    return {
        'phase': phase,
        'run_id': run_id,
        'workflow_id': workflow_id,
        'metrics': metrics
    }

def collect_basecalling_metrics(bucket: str, prefix: str) -> Dict:
    """Collect basecalling phase metrics."""

    metrics = {
        'total_reads': 0,
        'mean_quality': 0,
        'mean_length': 0,
        'total_bases': 0,
        'duplex_ratio': 0
    }

    # Read summary file if exists
    try:
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/summary.json'
        )
        summary = json.loads(response['Body'].read())

        metrics['total_reads'] = summary.get('total_reads', 0)
        metrics['mean_quality'] = summary.get('mean_quality', 0)
        metrics['mean_length'] = summary.get('mean_length', 0)
        metrics['total_bases'] = summary.get('total_bases', 0)
        metrics['duplex_ratio'] = summary.get('duplex_ratio', 0)

    except Exception as e:
        print(f'Error reading basecalling metrics: {e}')

    return metrics

def collect_qc_metrics(bucket: str, prefix: str) -> Dict:
    """Collect QC phase metrics."""

    metrics = {
        'qc_status': 'UNKNOWN',
        'mean_quality_score': 0,
        'reads_passed': 0,
        'reads_failed': 0,
        'n50': 0
    }

    try:
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/qc_summary.json'
        )
        summary = json.loads(response['Body'].read())

        metrics['qc_status'] = summary.get('qc_status', 'UNKNOWN')
        metrics['mean_quality_score'] = summary.get('mean_quality', 0)
        metrics['reads_passed'] = summary.get('reads_passed', 0)
        metrics['reads_failed'] = summary.get('reads_failed', 0)
        metrics['n50'] = summary.get('n50', 0)

    except Exception as e:
        print(f'Error reading QC metrics: {e}')

    return metrics

def collect_host_removal_metrics(bucket: str, prefix: str) -> Dict:
    """Collect host removal metrics."""

    metrics = {
        'reads_before': 0,
        'reads_after': 0,
        'depletion_rate': 0,
        'host_contamination': 0
    }

    try:
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/depletion_stats.json'
        )
        stats = json.loads(response['Body'].read())

        metrics['reads_before'] = stats.get('reads_before', 0)
        metrics['reads_after'] = stats.get('reads_after', 0)
        metrics['depletion_rate'] = stats.get('depletion_rate', 0)
        metrics['host_contamination'] = stats.get('host_percentage', 0)

    except Exception as e:
        print(f'Error reading host removal metrics: {e}')

    return metrics

def collect_pathogen_metrics(bucket: str, prefix: str) -> Dict:
    """Collect pathogen detection metrics."""

    metrics = {
        'pathogens_detected': 0,
        'pmda_pathogens': 0,
        'perv_detected': False,
        'critical_pathogens': []
    }

    try:
        # Read pathogen summary
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/pathogen_summary.json'
        )
        summary = json.loads(response['Body'].read())

        metrics['pathogens_detected'] = summary.get('total_pathogens', 0)
        metrics['pmda_pathogens'] = summary.get('pmda_pathogens', 0)
        metrics['perv_detected'] = summary.get('perv_detected', False)
        metrics['critical_pathogens'] = summary.get('critical_pathogens', [])

    except Exception as e:
        print(f'Error reading pathogen metrics: {e}')

    return metrics

def collect_quantification_metrics(bucket: str, prefix: str) -> Dict:
    """Collect quantification metrics."""

    metrics = {
        'spike_in_recovery': 0,
        'normalization_factor': 0,
        'highest_copy_number': 0
    }

    try:
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/absolute_quantification.json'
        )
        quant = json.loads(response['Body'].read())

        # Find highest copy number
        max_copies = 0
        for pathogen in quant.get('pathogens', {}).values():
            copies = pathogen.get('copies_per_ml', 0)
            if copies > max_copies:
                max_copies = copies

        metrics['highest_copy_number'] = max_copies

        # Get normalization info
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/normalized_quantification.json'
        )
        norm = json.loads(response['Body'].read())

        metrics['spike_in_recovery'] = norm.get('recovery_rate', 0)
        metrics['normalization_factor'] = norm.get('normalization_factor', 0)

    except Exception as e:
        print(f'Error reading quantification metrics: {e}')

    return metrics

def collect_reporting_metrics(bucket: str, prefix: str) -> Dict:
    """Collect reporting metrics."""

    metrics = {
        'reports_generated': 0,
        'pmda_compliant': False,
        'formats': []
    }

    try:
        response = s3.get_object(
            Bucket=bucket,
            Key=f'{prefix}/analysis_complete.json'
        )
        summary = json.loads(response['Body'].read())

        metrics['reports_generated'] = len(summary.get('reports_generated', []))
        metrics['pmda_compliant'] = summary.get('pmda_compliance', False)
        metrics['formats'] = summary.get('reports_generated', [])

    except Exception as e:
        print(f'Error reading reporting metrics: {e}')

    return metrics

def store_metrics(workflow_id: str, phase: str, metrics: Dict):
    """Store metrics in database."""

    sql = """
        INSERT INTO phase_metrics
        (workflow_id, phase_name, metrics, created_at)
        VALUES (:workflow_id, :phase_name, :metrics, NOW())
    """

    rds.execute_statement(
        resourceArn=CLUSTER_ARN,
        secretArn=SECRET_ARN,
        database=DATABASE,
        sql=sql,
        parameters=[
            {'name': 'workflow_id', 'value': {'longValue': workflow_id}},
            {'name': 'phase_name', 'value': {'stringValue': phase}},
            {'name': 'metrics', 'value': {'stringValue': json.dumps(metrics, default=str)}}
        ]
    )

def send_cloudwatch_metrics(phase: str, run_id: str, metrics: Dict):
    """Send metrics to CloudWatch."""

    namespace = 'MinION/Analysis'

    # Convert metrics to CloudWatch format
    metric_data = []

    for key, value in metrics.items():
        if isinstance(value, (int, float)):
            metric_data.append({
                'MetricName': f'{phase}_{key}',
                'Value': float(value),
                'Unit': 'None',
                'Dimensions': [
                    {'Name': 'Phase', 'Value': phase},
                    {'Name': 'RunID', 'Value': run_id}
                ]
            })

    if metric_data:
        cloudwatch.put_metric_data(
            Namespace=namespace,
            MetricData=metric_data[:20]  # CloudWatch limit is 20 metrics per call
        )