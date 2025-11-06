#!/usr/bin/env python3
"""
Calculate host depletion rate from alignment statistics
"""

import argparse
import json
import sys
from pathlib import Path
import pysam

def calculate_depletion(bam_file: Path) -> dict:
    """Calculate host depletion statistics."""

    bam = pysam.AlignmentFile(str(bam_file), "rb")

    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0

    for read in bam:
        total_reads += 1
        if read.is_unmapped:
            unmapped_reads += 1
        else:
            mapped_reads += 1

    bam.close()

    depletion_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0

    return {
        'total_reads': total_reads,
        'host_reads': mapped_reads,
        'non_host_reads': unmapped_reads,
        'host_depletion_rate': round(depletion_rate, 2),
        'non_host_percentage': round(100 - depletion_rate, 2)
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-bam', required=True, type=Path)
    parser.add_argument('-o', '--output', required=True, type=Path)
    parser.add_argument('-r', '--run-id', required=True)

    args = parser.parse_args()

    stats = calculate_depletion(args.input_bam)
    stats['run_id'] = args.run_id

    with open(args.output, 'w') as f:
        json.dump(stats, f, indent=2)

    print(f"Host depletion rate: {stats['host_depletion_rate']}%")
    print(f"Non-host reads: {stats['non_host_reads']:,} ({stats['non_host_percentage']}%)")

if __name__ == '__main__':
    main()