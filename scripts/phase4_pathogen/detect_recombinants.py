#!/usr/bin/env python3
"""
Detect PERV recombinants
Identifies potential recombination events between PERV subtypes
"""

import argparse
import json
from pathlib import Path
import pysam
from collections import defaultdict

def detect_recombination(bam_file: Path) -> dict:
    """Detect potential PERV recombination events."""

    # Validate BAM file exists and is readable
    if not bam_file.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    try:
        bam = pysam.AlignmentFile(str(bam_file), "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file {bam_file}: {e}")

    chimeric_reads = []
    split_alignments = defaultdict(list)

    for read in bam:
        if read.is_unmapped:
            continue

        # Check for supplementary alignments (chimeric reads)
        if read.is_supplementary:
            split_alignments[read.query_name].append({
                'reference': read.reference_name,
                'start': read.reference_start,
                'end': read.reference_end
            })

        # Check if read maps to multiple PERV types
        if read.has_tag('SA'):  # Supplementary alignment tag
            sa_tag = read.get_tag('SA')
            chimeric_reads.append({
                'read_name': read.query_name,
                'primary_ref': read.reference_name,
                'supplementary': sa_tag
            })

    bam.close()

    # Analyze split alignments for recombination
    recombinants = []
    for read_name, alignments in split_alignments.items():
        if len(alignments) > 1:
            refs = [a['reference'] for a in alignments]
            perv_types = set()
            for ref in refs:
                if 'PERV-A' in ref:
                    perv_types.add('PERV-A')
                elif 'PERV-B' in ref:
                    perv_types.add('PERV-B')
                elif 'PERV-C' in ref:
                    perv_types.add('PERV-C')

            if len(perv_types) > 1:
                recombinants.append({
                    'read': read_name,
                    'types': list(perv_types),
                    'recombination_type': 'x'.join(sorted(perv_types))
                })

    results = {
        'recombinant_detected': len(recombinants) > 0,
        'recombinant_count': len(recombinants),
        'chimeric_read_count': len(chimeric_reads)
    }

    if recombinants:
        # Count recombination types
        recomb_types = defaultdict(int)
        for r in recombinants:
            recomb_types[r['recombination_type']] += 1

        results['recombination_types'] = dict(recomb_types)
        results['predominant_recombination'] = max(recomb_types, key=recomb_types.get)

        # Estimate breakpoints (simplified)
        results['estimated_breakpoints'] = []
        for r in recombinants[:10]:  # Analyze first 10 recombinants
            alignments = split_alignments[r['read']]
            if len(alignments) >= 2:
                breakpoint = (alignments[0]['end'] + alignments[1]['start']) // 2
                results['estimated_breakpoints'].append(breakpoint)

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    results = detect_recombination(args.bam)
    results['run_id'] = args.run_id

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    if results['recombinant_detected']:
        print(f"WARNING: PERV recombinants detected! Count: {results['recombinant_count']}")
        print(f"Recombination types: {results.get('recombination_types', {})}")
    else:
        print("No PERV recombinants detected")

if __name__ == '__main__':
    main()