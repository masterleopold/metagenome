#!/usr/bin/env python3
"""
Check QC metrics from NanoPlot output and determine if run passes quality thresholds.
"""

import argparse
import json
import sys
from pathlib import Path
import re

def parse_nanostats(nanostats_file: Path) -> dict:
    """Parse NanoStats.txt file into metrics dictionary."""

    metrics = {}

    with open(nanostats_file) as f:
        for line in f:
            line = line.strip()
            if ':' not in line:
                continue

            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip()

            # Remove commas from numbers
            value = value.replace(',', '')

            # Try to convert to numeric
            try:
                if '.' in value:
                    metrics[key] = float(value)
                else:
                    metrics[key] = int(value)
            except ValueError:
                metrics[key] = value

    return metrics

def check_qc_thresholds(metrics: dict, min_quality: float) -> dict:
    """Check if metrics pass QC thresholds."""

    # Default thresholds (PMDA requirements)
    thresholds = {
        'min_reads': 10000,
        'min_mean_quality': min_quality,
        'min_median_quality': min_quality - 1,
        'min_n50': 200,
        'min_mean_length': 400
    }

    checks = {}

    # Check number of reads
    total_reads = metrics.get('Number of reads', 0)
    checks['reads_count'] = {
        'value': total_reads,
        'threshold': thresholds['min_reads'],
        'pass': total_reads >= thresholds['min_reads']
    }

    # Check mean quality
    mean_quality = metrics.get('Mean read quality', 0)
    checks['mean_quality'] = {
        'value': round(mean_quality, 2),
        'threshold': thresholds['min_mean_quality'],
        'pass': mean_quality >= thresholds['min_mean_quality']
    }

    # Check median quality
    median_quality = metrics.get('Median read quality', 0)
    checks['median_quality'] = {
        'value': round(median_quality, 2),
        'threshold': thresholds['min_median_quality'],
        'pass': median_quality >= thresholds['min_median_quality']
    }

    # Check N50
    n50 = metrics.get('Read length N50', 0)
    checks['n50'] = {
        'value': n50,
        'threshold': thresholds['min_n50'],
        'pass': n50 >= thresholds['min_n50']
    }

    # Check mean read length
    mean_length = metrics.get('Mean read length', 0)
    checks['mean_length'] = {
        'value': round(mean_length, 2),
        'threshold': thresholds['min_mean_length'],
        'pass': mean_length >= thresholds['min_mean_length']
    }

    # Overall pass/fail
    all_passed = all(check['pass'] for check in checks.values())

    return {
        'checks': checks,
        'qc_status': 'PASS' if all_passed else 'FAIL'
    }

def main():
    parser = argparse.ArgumentParser(
        description='Check QC metrics from NanoPlot output'
    )
    parser.add_argument('--input', required=True, type=Path,
                       help='NanoStats.txt file from NanoPlot')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output JSON file path')
    parser.add_argument('--min-quality', type=float, default=9.0,
                       help='Minimum mean quality score (default: 9.0)')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')

    args = parser.parse_args()

    # Validate input file exists
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Checking QC metrics for run {args.run_id}...")
    print(f"Input: {args.input}")
    print(f"Minimum quality: {args.min_quality}")

    # Parse metrics
    metrics = parse_nanostats(args.input)

    # Check thresholds
    qc_results = check_qc_thresholds(metrics, args.min_quality)

    # Build output
    output_data = {
        'run_id': args.run_id,
        'qc_status': qc_results['qc_status'],
        'checks': qc_results['checks'],
        'raw_metrics': metrics
    }

    # Write results
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)

    # Print summary
    print(f"\nQC Summary:")
    print(f"  Status: {qc_results['qc_status']}")
    print(f"\nChecks:")
    for check_name, check_data in qc_results['checks'].items():
        status_symbol = '✓' if check_data['pass'] else '✗'
        print(f"  {status_symbol} {check_name}: {check_data['value']} (threshold: {check_data['threshold']})")

    print(f"\nResults written to: {args.output}")

    # Exit with appropriate code
    if qc_results['qc_status'] == 'FAIL':
        print("\n[WARNING] QC checks failed")
        sys.exit(1)
    else:
        print("\nAll QC checks passed")
        sys.exit(0)

if __name__ == '__main__':
    main()
