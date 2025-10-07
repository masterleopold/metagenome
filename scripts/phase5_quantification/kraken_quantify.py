#!/usr/bin/env python3
"""
Quantify pathogens from Kraken2 report
"""

import argparse
import json
from pathlib import Path
import pandas as pd

def parse_kraken_report(report_file: Path) -> dict:
    """Parse Kraken2 report and extract quantification data."""

    # Read Kraken2 report
    df = pd.read_csv(report_file, sep='\t', header=None,
                     names=['percentage', 'reads_clade', 'reads_taxon',
                            'rank', 'ncbi_id', 'name'])

    # Calculate total reads
    total_reads = df[df['rank'] == 'U']['reads_clade'].sum()  # Unclassified
    total_reads += df[df['rank'] == 'R']['reads_clade'].sum()  # Root

    results = {
        'total_reads': int(total_reads),
        'pathogens': {}
    }

    # PMDA pathogen list (abbreviated)
    pmda_pathogens = {
        'Porcine endogenous retrovirus': 'PERV',
        'Hepatitis E virus': 'HEV',
        'Japanese encephalitis virus': 'JEV',
        'Streptococcus suis': 'SS',
        'Escherichia coli': 'EC'
        # Add all 91 pathogens
    }

    # Extract pathogen quantities
    for _, row in df.iterrows():
        name = row['name'].strip()

        for pathogen_name, code in pmda_pathogens.items():
            if pathogen_name.lower() in name.lower():
                reads = int(row['reads_clade'])
                rpm = (reads / total_reads * 1e6) if total_reads > 0 else 0

                results['pathogens'][code] = {
                    'name': pathogen_name,
                    'reads': reads,
                    'rpm': round(rpm, 2),
                    'percentage': round(row['percentage'], 4),
                    'rank': row['rank'],
                    'ncbi_id': row['ncbi_id']
                }

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    results = parse_kraken_report(args.report)
    results['run_id'] = args.run_id
    results['method'] = 'kraken2'

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"Total reads: {results['total_reads']:,}")
    print(f"Pathogens detected: {len(results['pathogens'])}")

    if results['pathogens']:
        print("\nTop pathogens by RPM:")
        sorted_pathogens = sorted(results['pathogens'].items(),
                                  key=lambda x: x[1]['rpm'], reverse=True)
        for code, data in sorted_pathogens[:5]:
            print(f"  {code}: {data['rpm']} RPM ({data['reads']} reads)")

if __name__ == '__main__':
    main()