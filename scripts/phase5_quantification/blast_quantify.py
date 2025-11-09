#!/usr/bin/env python3
"""
Quantify pathogens from BLAST results
"""

import argparse
import json
from pathlib import Path
from collections import defaultdict

def parse_blast_results(blast_file: Path) -> dict:
    """Parse BLAST results and quantify pathogens."""

    pathogen_hits = defaultdict(int)
    pathogen_scores = defaultdict(list)

    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue

            subject_title = parts[12]
            evalue = float(parts[10])
            bitscore = float(parts[11])

            # Check for PMDA pathogens
            if 'endogenous retrovirus' in subject_title.lower():
                pathogen_hits['PERV'] += 1
                pathogen_scores['PERV'].append(bitscore)
            elif 'hepatitis e' in subject_title.lower():
                pathogen_hits['HEV'] += 1
                pathogen_scores['HEV'].append(bitscore)
            # Add other pathogens...

    results = {
        'pathogens': {}
    }

    for pathogen, count in pathogen_hits.items():
        scores = pathogen_scores[pathogen]
        avg_score = sum(scores) / len(scores) if len(scores) > 0 else 0.0
        results['pathogens'][pathogen] = {
            'hit_count': count,
            'avg_bitscore': round(avg_score, 2)
        }

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    results = parse_blast_results(args.blast)
    results['run_id'] = args.run_id
    results['method'] = 'blast'

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == '__main__':
    main()