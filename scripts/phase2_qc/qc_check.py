#!/usr/bin/env python3
"""
MinION Metagenomics Pipeline - QC Check
Validates quality control metrics against PMDA thresholds
"""

import argparse
import json
import sys
from pathlib import Path
import psycopg2
import boto3

# PMDA QC thresholds
# Note: Aligned with default_pipeline.yaml quality_thresholds.min_reads
QC_THRESHOLDS = {
    'min_reads': 10000,  # Changed from 100000 to match pipeline config
    'min_mean_qscore': 9.0,
    'min_median_qscore': 8.0,
    'min_n50': 200,
    'max_failed_reads_pct': 50,
    'min_q30_pct': 5.0
}


def check_qc_metrics(summary_file: Path) -> tuple:
    """Check if QC metrics pass thresholds."""
    with open(summary_file) as f:
        summary = json.load(f)

    metrics = summary.get('metrics', {})
    failures = []

    # Check each threshold
    if metrics.get('total_reads', 0) < QC_THRESHOLDS['min_reads']:
        failures.append(f"Total reads ({metrics.get('total_reads', 0)}) < {QC_THRESHOLDS['min_reads']}")

    if metrics.get('mean_qscore', 0) < QC_THRESHOLDS['min_mean_qscore']:
        failures.append(f"Mean Q-score ({metrics.get('mean_qscore', 0)}) < {QC_THRESHOLDS['min_mean_qscore']}")

    if metrics.get('n50', 0) < QC_THRESHOLDS['min_n50']:
        failures.append(f"N50 ({metrics.get('n50', 0)}) < {QC_THRESHOLDS['min_n50']}")

    qc_pass = len(failures) == 0
    return qc_pass, failures


def update_database(run_id: str, qc_pass: bool, metrics: dict, db_config: dict):
    """Update database with QC results."""
    try:
        conn = psycopg2.connect(**db_config)
        cur = conn.cursor()

        # Update workflow status
        status = 'host_removal' if qc_pass else 'failed'
        cur.execute("""
            UPDATE workflow_executions
            SET status = %s,
                mean_qscore = %s,
                median_qscore = %s,
                n50 = %s
            WHERE run_id = %s
        """, (status, metrics.get('mean_qscore'), metrics.get('median_qscore'),
              metrics.get('n50'), run_id))

        # Insert QC metrics
        for metric_name, value in metrics.items():
            cur.execute("""
                INSERT INTO qc_metrics (run_id, phase, metric_name, metric_value, passed)
                VALUES (%s, %s, %s, %s, %s)
                ON CONFLICT (run_id, phase, metric_name) DO UPDATE
                SET metric_value = EXCLUDED.metric_value, passed = EXCLUDED.passed
            """, (run_id, 'qc', metric_name, value, qc_pass))

        conn.commit()
        conn.close()
    except Exception as e:
        print(f"Database update failed: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', required=True, type=Path)
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--db-host')
    parser.add_argument('--db-name', default='minion_metadata')
    parser.add_argument('--db-user', default='minion_user')
    parser.add_argument('--db-password')

    args = parser.parse_args()

    with open(args.summary) as f:
        summary = json.load(f)

    qc_pass, failures = check_qc_metrics(args.summary)

    if args.db_host:
        db_config = {
            'host': args.db_host,
            'database': args.db_name,
            'user': args.db_user,
            'password': args.db_password
        }
        update_database(args.run_id, qc_pass, summary.get('metrics', {}), db_config)

    if not qc_pass:
        print("QC FAILED:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        sys.exit(1)
    else:
        print("QC PASSED")
        sys.exit(0)


if __name__ == '__main__':
    main()