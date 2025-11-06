#!/usr/bin/env python3
"""
Calculate host depletion metrics by comparing read counts before and after host removal.
"""

import argparse
import json
import sys
from pathlib import Path
import glob

def count_fastq_reads(directory: Path) -> int:
    """Count total reads in all FASTQ files in a directory."""
    total_reads = 0

    # Support both .fastq and .fastq.gz
    patterns = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']

    for pattern in patterns:
        for fastq_file in directory.glob(pattern):
            # Count lines and divide by 4 (FASTQ format)
            try:
                if str(fastq_file).endswith('.gz'):
                    import gzip
                    with gzip.open(fastq_file, 'rt') as f:
                        line_count = sum(1 for _ in f)
                else:
                    with open(fastq_file) as f:
                        line_count = sum(1 for _ in f)

                total_reads += line_count // 4
            except Exception as e:
                print(f"Warning: Failed to count reads in {fastq_file}: {e}", file=sys.stderr)

    return total_reads

def calculate_depletion_metrics(before_dir: Path, after_dir: Path, run_id: str) -> dict:
    """Calculate host depletion statistics."""

    reads_before = count_fastq_reads(before_dir)
    reads_after = count_fastq_reads(after_dir)

    if reads_before == 0:
        print("ERROR: No reads found in input directory", file=sys.stderr)
        sys.exit(1)

    host_reads = reads_before - reads_after
    depletion_rate = (host_reads / reads_before * 100) if reads_before > 0 else 0
    retention_rate = (reads_after / reads_before * 100) if reads_before > 0 else 0

    # PMDA requirement: >90% host depletion efficiency
    pmda_compliant = depletion_rate >= 90.0

    metrics = {
        'run_id': run_id,
        'reads_before_depletion': reads_before,
        'reads_after_depletion': reads_after,
        'host_reads_removed': host_reads,
        'depletion_rate_percent': round(depletion_rate, 2),
        'retention_rate_percent': round(retention_rate, 2),
        'pmda_compliant': pmda_compliant,
        'pmda_threshold_percent': 90.0
    }

    return metrics

def main():
    parser = argparse.ArgumentParser(
        description='Calculate host depletion metrics for MinION pipeline'
    )
    parser.add_argument('--before', required=True, type=Path,
                       help='Directory containing FASTQ files before host removal')
    parser.add_argument('--after', required=True, type=Path,
                       help='Directory containing FASTQ files after host removal')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output JSON file path')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')

    args = parser.parse_args()

    # Validate directories exist
    if not args.before.exists():
        print(f"ERROR: Input directory not found: {args.before}", file=sys.stderr)
        sys.exit(1)

    if not args.after.exists():
        print(f"ERROR: Output directory not found: {args.after}", file=sys.stderr)
        sys.exit(1)

    # Calculate metrics
    print(f"Calculating depletion metrics for run {args.run_id}...")
    print(f"Before: {args.before}")
    print(f"After: {args.after}")

    metrics = calculate_depletion_metrics(args.before, args.after, args.run_id)

    # Write results
    with open(args.output, 'w') as f:
        json.dump(metrics, f, indent=2)

    # Print summary
    print(f"\nDepletion Summary:")
    print(f"  Reads before: {metrics['reads_before_depletion']:,}")
    print(f"  Reads after:  {metrics['reads_after_depletion']:,}")
    print(f"  Host removed: {metrics['host_reads_removed']:,}")
    print(f"  Depletion rate: {metrics['depletion_rate_percent']}%")
    print(f"  PMDA compliant: {'YES' if metrics['pmda_compliant'] else 'NO'}")

    if not metrics['pmda_compliant']:
        print(f"\n[WARNING] Depletion rate below PMDA threshold of {metrics['pmda_threshold_percent']}%")
        sys.exit(2)  # Warning exit code

    print(f"\nResults written to: {args.output}")

if __name__ == '__main__':
    main()
