#!/usr/bin/env python3
"""
Extract PMDA-designated 91 pathogens from Kraken2 report
Filters Kraken2 output for regulatory compliance tracking
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Set
import sys


def load_pmda_pathogens(config_file: Path) -> Dict[str, Dict]:
    """Load PMDA pathogen database."""
    with open(config_file) as f:
        pmda_data = json.load(f)

    # Build lookup dictionary: taxon_name -> pathogen_info
    pathogen_lookup = {}
    pathogen_codes = {}

    for category, data in pmda_data['categories'].items():
        for pathogen in data.get('pathogens', []):
            name = pathogen['name'].lower()
            code = pathogen['code']
            pathogen_lookup[name] = {
                'code': code,
                'name': pathogen['name'],
                'category': category,
                'risk_level': pathogen['risk_level']
            }
            pathogen_codes[code] = pathogen

    return pathogen_lookup, pathogen_codes


def parse_kraken2_report(report_file: Path) -> List[Dict]:
    """
    Parse Kraken2 report format.

    Format:
    Percentage  Clade_reads  Taxon_reads  Rank  TaxID  Name
    """
    results = []

    with open(report_file) as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            try:
                percentage = float(parts[0])
                clade_reads = int(parts[1])
                taxon_reads = int(parts[2])
                rank = parts[3].strip()
                taxid = parts[4].strip()
                name = parts[5].strip()

                results.append({
                    'percentage': percentage,
                    'clade_reads': clade_reads,
                    'taxon_reads': taxon_reads,
                    'rank': rank,
                    'taxid': taxid,
                    'name': name
                })
            except (ValueError, IndexError) as e:
                # Skip malformed lines
                continue

    return results


def match_pmda_pathogens(kraken_results: List[Dict],
                         pathogen_lookup: Dict[str, Dict],
                         pathogen_codes: Dict[str, Dict],
                         min_reads: int = 10) -> Dict:
    """
    Match Kraken2 results to PMDA pathogen list.

    Uses fuzzy matching for genus/species names.
    """
    detected_pathogens = {}
    total_classified_reads = 0

    for result in kraken_results:
        name_lower = result['name'].lower().strip()
        reads = result['taxon_reads']

        if reads < min_reads:
            continue

        # Track total classified reads
        if result['rank'] == 'U':  # Unclassified
            continue
        total_classified_reads += reads

        # Direct name match
        if name_lower in pathogen_lookup:
            pathogen = pathogen_lookup[name_lower]
            code = pathogen['code']

            if code not in detected_pathogens:
                detected_pathogens[code] = {
                    'code': code,
                    'name': pathogen['name'],
                    'category': pathogen['category'],
                    'risk_level': pathogen['risk_level'],
                    'reads': 0,
                    'percentage': 0.0,
                    'rank': result['rank'],
                    'taxid': result['taxid'],
                    'detection_confidence': 'high'
                }

            # Accumulate reads if multiple entries
            detected_pathogens[code]['reads'] += reads
            detected_pathogens[code]['percentage'] += result['percentage']

        # Partial name matching for genus/species
        else:
            for pathogen_name, pathogen in pathogen_lookup.items():
                # Check if pathogen name is substring of Kraken result
                # or vice versa (handles partial matches)
                if (pathogen_name in name_lower or
                    name_lower in pathogen_name or
                    any(word in name_lower.split() for word in pathogen_name.split() if len(word) > 4)):

                    code = pathogen['code']

                    if code not in detected_pathogens:
                        detected_pathogens[code] = {
                            'code': code,
                            'name': pathogen['name'],
                            'category': pathogen['category'],
                            'risk_level': pathogen['risk_level'],
                            'reads': 0,
                            'percentage': 0.0,
                            'rank': result['rank'],
                            'taxid': result['taxid'],
                            'detection_confidence': 'medium',
                            'matched_name': result['name']
                        }

                    detected_pathogens[code]['reads'] += reads
                    detected_pathogens[code]['percentage'] += result['percentage']
                    break

    # Calculate detection statistics
    detected_critical = [p for p in detected_pathogens.values()
                        if p['risk_level'] == 'CRITICAL']
    detected_high = [p for p in detected_pathogens.values()
                    if p['risk_level'] == 'HIGH']

    summary = {
        'total_pmda_pathogens_detected': len(detected_pathogens),
        'critical_pathogens_detected': len(detected_critical),
        'high_risk_pathogens_detected': len(detected_high),
        'pathogens': detected_pathogens,
        'total_classified_reads': total_classified_reads
    }

    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Extract PMDA-designated pathogens from Kraken2 report'
    )
    parser.add_argument('--report', '--kraken', dest='report', required=True, type=Path,
                       help='Kraken2 report file')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output JSON file')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('--config', type=Path,
                       default='/home/user/metagenome/templates/config/pmda_pathogens.json',
                       help='PMDA pathogen configuration file')
    parser.add_argument('--min-reads', type=int, default=10,
                       help='Minimum reads for detection (default: 10)')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Validate inputs
    if not args.report.exists():
        print(f"ERROR: Kraken2 report not found: {args.report}", file=sys.stderr)
        sys.exit(1)

    # Try multiple locations for config file
    config_file = args.config
    if not config_file.exists():
        # Try relative to script location
        script_dir = Path(__file__).parent.parent.parent
        config_file = script_dir / 'templates' / 'config' / 'pmda_pathogens.json'

    if not config_file.exists():
        print(f"ERROR: PMDA pathogen config not found: {config_file}", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Loading PMDA pathogen database from {config_file}")

    # Load PMDA pathogen database
    pathogen_lookup, pathogen_codes = load_pmda_pathogens(config_file)

    if args.verbose:
        print(f"Loaded {len(pathogen_lookup)} PMDA pathogens")

    # Parse Kraken2 report
    if args.verbose:
        print(f"Parsing Kraken2 report: {args.report}")

    kraken_results = parse_kraken2_report(args.report)

    if args.verbose:
        print(f"Parsed {len(kraken_results)} taxa from Kraken2 report")

    # Match to PMDA pathogens
    results = match_pmda_pathogens(
        kraken_results, pathogen_lookup, pathogen_codes, args.min_reads
    )

    # Add metadata
    results['run_id'] = args.run_id
    results['kraken2_report'] = str(args.report)
    results['min_reads_threshold'] = args.min_reads
    results['pmda_version'] = '2024.1'

    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"PMDA Pathogen Detection Summary:")
    print(f"  Total PMDA pathogens detected: {results['total_pmda_pathogens_detected']}")
    print(f"  Critical pathogens detected: {results['critical_pathogens_detected']}")
    print(f"  High-risk pathogens detected: {results['high_risk_pathogens_detected']}")

    if results['critical_pathogens_detected'] > 0:
        print("\nCRITICAL PATHOGENS DETECTED:")
        for code, pathogen in results['pathogens'].items():
            if pathogen['risk_level'] == 'CRITICAL':
                print(f"  - {code}: {pathogen['name']} ({pathogen['reads']} reads)")

    print(f"\nResults written to: {args.output}")

    # Exit with warning if critical pathogens detected
    if results['critical_pathogens_detected'] > 0:
        sys.exit(2)  # Warning exit code

    sys.exit(0)


if __name__ == '__main__':
    main()
