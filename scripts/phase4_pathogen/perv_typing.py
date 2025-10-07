#!/usr/bin/env python3
"""
PERV typing from aligned reads
Identifies PERV-A, PERV-B, and PERV-C subtypes
"""

import argparse
import json
from pathlib import Path
from collections import defaultdict
import pysam
import numpy as np

# PERV subtype specific regions
PERV_MARKERS = {
    'PERV-A': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC']
    },
    'PERV-B': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAACCACCGTAGC', 'TGGAAACCTGGAAAACC']
    },
    'PERV-C': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCATAGG', 'TGGAGACCTGGAAGAAC']
    }
}


def identify_perv_type(bam_file: Path) -> dict:
    """Identify PERV subtypes from aligned reads."""

    bam = pysam.AlignmentFile(str(bam_file), "rb")
    perv_reads = defaultdict(list)
    perv_coverage = defaultdict(lambda: defaultdict(int))

    for read in bam:
        if read.is_unmapped:
            continue

        reference = read.reference_name
        if 'PERV' not in reference:
            continue

        # Count reads per PERV type
        for perv_type in ['PERV-A', 'PERV-B', 'PERV-C']:
            if perv_type in reference:
                perv_reads[perv_type].append(read.query_name)

                # Calculate coverage
                for pos in range(read.reference_start, read.reference_end):
                    perv_coverage[perv_type][pos] += 1

        # Check for specific motifs
        seq = read.query_sequence
        if seq:
            for perv_type, markers in PERV_MARKERS.items():
                for motif in markers['specific_motifs']:
                    if motif in seq:
                        perv_reads[f"{perv_type}_specific"].append(read.query_name)

    bam.close()

    # Calculate statistics
    results = {
        'perv_a_detected': len(perv_reads['PERV-A']) > 0,
        'perv_b_detected': len(perv_reads['PERV-B']) > 0,
        'perv_c_detected': len(perv_reads['PERV-C']) > 0,
        'perv_a_reads': len(set(perv_reads['PERV-A'])),
        'perv_b_reads': len(set(perv_reads['PERV-B'])),
        'perv_c_reads': len(set(perv_reads['PERV-C'])),
    }

    # Calculate coverage breadth
    for perv_type in ['PERV-A', 'PERV-B', 'PERV-C']:
        if perv_coverage[perv_type]:
            coverage_array = list(perv_coverage[perv_type].values())
            results[f"{perv_type.lower().replace('-', '_')}_coverage_depth"] = np.mean(coverage_array)
            results[f"{perv_type.lower().replace('-', '_')}_coverage_breadth"] = (
                len([x for x in coverage_array if x > 0]) / 8000 * 100  # Assuming ~8kb PERV genome
            )

    # Determine predominant type
    read_counts = {
        'PERV-A': results['perv_a_reads'],
        'PERV-B': results['perv_b_reads'],
        'PERV-C': results['perv_c_reads']
    }

    if max(read_counts.values()) > 0:
        results['predominant_type'] = max(read_counts, key=read_counts.get)
    else:
        results['predominant_type'] = 'None'

    # Check for co-infection
    detected_types = [k for k, v in read_counts.items() if v > 10]
    if len(detected_types) > 1:
        results['co_infection'] = True
        results['co_infection_types'] = detected_types
    else:
        results['co_infection'] = False

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    results = identify_perv_type(args.bam)
    results['run_id'] = args.run_id

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"PERV-A: {'Detected' if results['perv_a_detected'] else 'Not detected'} ({results['perv_a_reads']} reads)")
    print(f"PERV-B: {'Detected' if results['perv_b_detected'] else 'Not detected'} ({results['perv_b_reads']} reads)")
    print(f"PERV-C: {'Detected' if results['perv_c_detected'] else 'Not detected'} ({results['perv_c_reads']} reads)")

    if results.get('co_infection'):
        print(f"WARNING: Co-infection detected with {', '.join(results['co_infection_types'])}")


if __name__ == '__main__':
    main()