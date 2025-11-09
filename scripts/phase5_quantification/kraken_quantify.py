#!/usr/bin/env python3
"""
Quantify pathogens from Kraken2 report

Processes all 91 PMDA-designated pathogens for xenotransplantation screening.
"""

import argparse
import json
from pathlib import Path
import pandas as pd
from typing import Dict, Tuple


def load_pmda_pathogens(config_file: Path) -> Tuple[Dict[str, str], Dict[str, Dict]]:
    """
    Load PMDA pathogen database from config file.

    Args:
        config_file: Path to pmda_pathogens.json

    Returns:
        Tuple of (name_to_code mapping, code_to_info mapping)
    """
    with open(config_file) as f:
        pmda_data = json.load(f)

    name_to_code = {}
    code_to_info = {}

    for category, data in pmda_data['categories'].items():
        for pathogen in data.get('pathogens', []):
            name = pathogen['name']
            code = pathogen['code']

            # Store both full name and common variations
            name_to_code[name.lower()] = code

            # Also add genus-level matching for bacteria/fungi
            if category in ['bacteria', 'fungi', 'parasites']:
                # Extract genus (first word)
                genus = name.split()[0].lower()
                if genus not in name_to_code:  # Don't overwrite if already exists
                    name_to_code[genus] = code

            code_to_info[code] = {
                'name': name,
                'category': category,
                'risk_level': pathogen['risk_level']
            }

    return name_to_code, code_to_info


def parse_kraken_report(report_file: Path, config_file: Path) -> dict:
    """
    Parse Kraken2 report and extract quantification data for all PMDA pathogens.

    Args:
        report_file: Path to Kraken2 report file
        config_file: Path to PMDA pathogen config

    Returns:
        Dictionary with quantification results for all detected PMDA pathogens
    """
    # Load PMDA pathogen database (all 91 pathogens)
    name_to_code, code_to_info = load_pmda_pathogens(config_file)

    # Read Kraken2 report
    df = pd.read_csv(report_file, sep='\t', header=None,
                     names=['percentage', 'reads_clade', 'reads_taxon',
                            'rank', 'ncbi_id', 'name'])

    # Calculate total reads
    total_reads = df[df['rank'] == 'U']['reads_clade'].sum()  # Unclassified
    total_reads += df[df['rank'] == 'R']['reads_clade'].sum()  # Root

    results = {
        'total_reads': int(total_reads),
        'pathogens': {},
        'pmda_pathogens_detected': 0,
        'pmda_config_loaded': len(code_to_info)  # Should be 91
    }

    # Extract pathogen quantities using fuzzy matching
    for _, row in df.iterrows():
        taxon_name = row['name'].strip().lower()

        # Try to match against PMDA pathogen list
        matched_code = None
        matched_name = None

        # Exact match first
        if taxon_name in name_to_code:
            matched_code = name_to_code[taxon_name]
            matched_name = code_to_info[matched_code]['name']
        else:
            # Fuzzy match - check if any PMDA pathogen name is contained in taxon
            for pmda_name, code in name_to_code.items():
                if pmda_name in taxon_name or taxon_name in pmda_name:
                    matched_code = code
                    matched_name = code_to_info[code]['name']
                    break

        if matched_code:
            reads = int(row['reads_clade'])
            rpm = (reads / total_reads * 1e6) if total_reads > 0 else 0

            # Only add if we have more reads than existing entry (prefer species over genus)
            if matched_code not in results['pathogens'] or reads > results['pathogens'][matched_code]['reads']:
                results['pathogens'][matched_code] = {
                    'name': matched_name,
                    'code': matched_code,
                    'reads': reads,
                    'rpm': round(rpm, 2),
                    'percentage': round(row['percentage'], 4),
                    'rank': row['rank'],
                    'ncbi_id': row['ncbi_id'],
                    'kraken_taxon': row['name'].strip(),
                    'category': code_to_info[matched_code]['category'],
                    'risk_level': code_to_info[matched_code]['risk_level']
                }

    results['pmda_pathogens_detected'] = len(results['pathogens'])

    return results

def main():
    parser = argparse.ArgumentParser(
        description='Quantify all 91 PMDA pathogens from Kraken2 report'
    )
    parser.add_argument('--report', required=True, type=Path,
                       help='Kraken2 report file')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output JSON file')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('--config', type=Path,
                       help='PMDA pathogen config file (default: auto-detect)')

    args = parser.parse_args()

    # Auto-detect config file if not provided
    if args.config is None:
        # Try common locations
        possible_paths = [
            Path('/opt/minion/templates/config/pmda_pathogens.json'),  # AMI location
            Path(__file__).parent.parent.parent / 'templates' / 'config' / 'pmda_pathogens.json',  # Repo location
            Path.cwd() / 'templates' / 'config' / 'pmda_pathogens.json'  # Working directory
        ]
        for path in possible_paths:
            if path.exists():
                args.config = path
                break

        if args.config is None:
            print("ERROR: Could not find pmda_pathogens.json config file")
            print("Please specify with --config parameter")
            return 1

    # Parse Kraken2 report with PMDA pathogen matching
    results = parse_kraken_report(args.report, args.config)
    results['run_id'] = args.run_id
    results['method'] = 'kraken2'
    results['config_file'] = str(args.config)

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"Loaded PMDA pathogen database: {results['pmda_config_loaded']} pathogens")
    print(f"Total reads in Kraken2 report: {results['total_reads']:,}")
    print(f"PMDA pathogens detected: {results['pmda_pathogens_detected']}")

    if results['pathogens']:
        print("\nTop PMDA pathogens by RPM:")
        sorted_pathogens = sorted(results['pathogens'].items(),
                                  key=lambda x: x[1]['rpm'], reverse=True)
        for code, data in sorted_pathogens[:10]:
            risk_marker = "[CRITICAL]" if data['risk_level'] == 'CRITICAL' else ""
            print(f"  {code} ({data['category']}): {data['rpm']} RPM ({data['reads']} reads) {risk_marker}")

        # Alert if critical pathogens detected
        critical_detected = [code for code, data in results['pathogens'].items()
                           if data['risk_level'] == 'CRITICAL']
        if critical_detected:
            print(f"\n[WARNING] {len(critical_detected)} CRITICAL pathogen(s) detected: {', '.join(critical_detected)}")

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())