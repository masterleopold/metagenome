#!/usr/bin/env python3
"""
Aggregate pathogen detection results from multiple sources
Combines Kraken2, BLAST, and PERV-specific analysis
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Set
from collections import defaultdict


def load_json_safe(file_path: Path) -> Dict:
    """Safely load JSON file, return empty dict if not found."""
    if not file_path.exists():
        return {}

    try:
        with open(file_path) as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"WARNING: Failed to parse {file_path}: {e}", file=sys.stderr)
        return {}


def aggregate_kraken_results(kraken_dir: Path) -> Dict:
    """Aggregate Kraken2 results."""
    results = {
        'method': 'kraken2',
        'pathogens': {},
        'pmda_pathogens': {}
    }

    # Load PMDA pathogen extraction results
    pmda_file = kraken_dir / 'pmda_pathogens.json'
    if pmda_file.exists():
        pmda_data = load_json_safe(pmda_file)
        results['pmda_pathogens'] = pmda_data.get('pathogens', {})
        results['pmda_summary'] = {
            'total_detected': pmda_data.get('total_pmda_pathogens_detected', 0),
            'critical_detected': pmda_data.get('critical_pathogens_detected', 0),
            'high_risk_detected': pmda_data.get('high_risk_pathogens_detected', 0)
        }

    # Load Bracken species-level results
    bracken_file = kraken_dir / 'bracken_S.txt'
    if bracken_file.exists():
        try:
            with open(bracken_file) as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        name = parts[0]
                        taxid = parts[1]
                        reads = int(parts[5])
                        results['pathogens'][name] = {
                            'taxid': taxid,
                            'reads': reads,
                            'method': 'kraken2_bracken'
                        }
        except Exception as e:
            print(f"WARNING: Failed to parse Bracken results: {e}", file=sys.stderr)

    return results


def aggregate_blast_results(blast_dir: Path) -> Dict:
    """Aggregate BLAST results."""
    results = {
        'method': 'blast',
        'viral_hits': [],
        'top_hits': []
    }

    # Look for BLAST output files
    blast_files = list(blast_dir.glob('*.blast.txt')) + list(blast_dir.glob('blast_results.txt'))

    for blast_file in blast_files:
        try:
            with open(blast_file) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue

                    # BLAST outfmt 6 format
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        hit = {
                            'query': parts[0],
                            'subject': parts[1],
                            'identity': float(parts[2]),
                            'length': int(parts[3]),
                            'evalue': float(parts[10]),
                            'bitscore': float(parts[11])
                        }

                        if hit['identity'] >= 90 and hit['evalue'] < 1e-5:
                            results['viral_hits'].append(hit)

        except Exception as e:
            print(f"WARNING: Failed to parse BLAST file {blast_file}: {e}", file=sys.stderr)

    # Get top hits by bitscore
    if results['viral_hits']:
        sorted_hits = sorted(results['viral_hits'],
                           key=lambda x: x['bitscore'],
                           reverse=True)
        results['top_hits'] = sorted_hits[:50]

    return results


def aggregate_perv_results(perv_dir: Path) -> Dict:
    """Aggregate PERV-specific analysis results."""
    results = {
        'method': 'perv_specific',
        'perv_detected': False,
        'perv_types': {},
        'recombinants': {},
        'phylogenetics': {}
    }

    # Load PERV typing results
    typing_file = perv_dir / 'perv_types.json'
    if typing_file.exists():
        typing_data = load_json_safe(typing_file)
        results['perv_types'] = typing_data
        results['perv_detected'] = (
            typing_data.get('perv_a_detected', False) or
            typing_data.get('perv_b_detected', False) or
            typing_data.get('perv_c_detected', False)
        )

    # Load recombinant detection
    recombinant_file = perv_dir / 'perv_recombinants.json'
    if recombinant_file.exists():
        results['recombinants'] = load_json_safe(recombinant_file)

    # Load phylogenetic analysis
    phylo_file = perv_dir / 'perv_phylogenetics.json'
    if phylo_file.exists():
        results['phylogenetics'] = load_json_safe(phylo_file)

    # Load summary
    summary_file = perv_dir / 'perv_summary.json'
    if summary_file.exists():
        summary = load_json_safe(summary_file)
        results['perv_detected'] = summary.get('perv_detected', False)

    return results


def merge_detections(kraken_results: Dict, blast_results: Dict,
                    perv_results: Dict) -> Dict:
    """
    Merge all detection results with conflict resolution.

    Priority: PERV-specific > Kraken2 PMDA > BLAST
    """
    merged = {
        'detection_methods_used': [],
        'all_pathogens': {},
        'pmda_pathogens': {},
        'perv_analysis': {},
        'critical_findings': [],
        'summary': {}
    }

    # Track which methods were used
    if kraken_results.get('pathogens'):
        merged['detection_methods_used'].append('kraken2')
    if blast_results.get('viral_hits'):
        merged['detection_methods_used'].append('blast')
    if perv_results.get('perv_detected'):
        merged['detection_methods_used'].append('perv_specific')

    # Add PMDA pathogens from Kraken2
    if kraken_results.get('pmda_pathogens'):
        merged['pmda_pathogens'] = kraken_results['pmda_pathogens']

        # Add to critical findings if critical pathogens detected
        for code, pathogen in kraken_results['pmda_pathogens'].items():
            if pathogen['risk_level'] == 'CRITICAL':
                merged['critical_findings'].append({
                    'type': 'critical_pathogen',
                    'pathogen_code': code,
                    'pathogen_name': pathogen['name'],
                    'reads': pathogen['reads'],
                    'method': 'kraken2'
                })

    # Add PERV results (highest priority)
    if perv_results:
        merged['perv_analysis'] = perv_results

        if perv_results.get('perv_detected'):
            merged['critical_findings'].append({
                'type': 'perv_detection',
                'severity': 'CRITICAL',
                'message': 'PERV sequences detected - xenotransplantation safety concern',
                'details': perv_results.get('perv_types', {})
            })

    # Add all Kraken2 pathogens
    for name, data in kraken_results.get('pathogens', {}).items():
        merged['all_pathogens'][name] = {
            'name': name,
            'reads': data['reads'],
            'detection_method': 'kraken2',
            'confidence': 'high'
        }

    # Add BLAST viral hits
    if blast_results.get('top_hits'):
        merged['blast_viral_hits'] = blast_results['top_hits'][:20]

    # Generate summary statistics
    merged['summary'] = {
        'total_pathogens_detected': len(merged['all_pathogens']),
        'pmda_pathogens_detected': len(merged['pmda_pathogens']),
        'critical_findings_count': len(merged['critical_findings']),
        'perv_detected': perv_results.get('perv_detected', False),
        'methods_used': merged['detection_methods_used']
    }

    # Add PMDA summary if available
    if kraken_results.get('pmda_summary'):
        merged['summary']['pmda_summary'] = kraken_results['pmda_summary']

    return merged


def generate_alert_recommendations(merged_results: Dict) -> List[Dict]:
    """Generate alert recommendations based on findings."""
    alerts = []

    # PERV detection (highest priority)
    if merged_results['summary'].get('perv_detected'):
        alerts.append({
            'priority': 'CRITICAL',
            'alert_type': 'PERV_DETECTION',
            'action': 'IMMEDIATE',
            'message': 'PERV sequences detected - quarantine donor immediately',
            'details': merged_results.get('perv_analysis', {})
        })

    # Critical PMDA pathogens
    critical_count = merged_results['summary']['critical_findings_count']
    if critical_count > 0:
        for finding in merged_results['critical_findings']:
            if finding['type'] == 'critical_pathogen':
                alerts.append({
                    'priority': 'CRITICAL',
                    'alert_type': 'CRITICAL_PATHOGEN',
                    'action': 'IMMEDIATE',
                    'pathogen': finding['pathogen_name'],
                    'pathogen_code': finding['pathogen_code'],
                    'reads': finding['reads']
                })

    # High pathogen count
    if merged_results['summary']['total_pathogens_detected'] > 50:
        alerts.append({
            'priority': 'WARNING',
            'alert_type': 'HIGH_PATHOGEN_COUNT',
            'action': 'REVIEW',
            'message': f"High pathogen count detected: {merged_results['summary']['total_pathogens_detected']}"
        })

    return alerts


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate pathogen detection results from multiple sources'
    )
    parser.add_argument('--kraken', required=True, type=Path,
                       help='Kraken2 results directory')
    parser.add_argument('--blast', required=True, type=Path,
                       help='BLAST results directory')
    parser.add_argument('--perv', required=True, type=Path,
                       help='PERV analysis directory')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output JSON file')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    if args.verbose:
        print(f"Aggregating results for run: {args.run_id}")
        print(f"  Kraken2 dir: {args.kraken}")
        print(f"  BLAST dir: {args.blast}")
        print(f"  PERV dir: {args.perv}")

    # Aggregate results from each source
    kraken_results = aggregate_kraken_results(args.kraken)
    blast_results = aggregate_blast_results(args.blast)
    perv_results = aggregate_perv_results(args.perv)

    if args.verbose:
        print(f"\nKraken2 results:")
        print(f"  PMDA pathogens: {len(kraken_results.get('pmda_pathogens', {}))}")
        print(f"  Total pathogens: {len(kraken_results.get('pathogens', {}))}")

        print(f"\nBLAST results:")
        print(f"  Viral hits: {len(blast_results.get('viral_hits', []))}")

        print(f"\nPERV results:")
        print(f"  PERV detected: {perv_results.get('perv_detected', False)}")

    # Merge all results
    merged_results = merge_detections(kraken_results, blast_results, perv_results)

    # Generate alerts
    alerts = generate_alert_recommendations(merged_results)
    merged_results['alerts'] = alerts

    # Add metadata
    merged_results['run_id'] = args.run_id
    merged_results['aggregation_version'] = '1.0.0'

    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(merged_results, f, indent=2)

    # Print summary
    print("\n" + "="*70)
    print(f"PATHOGEN DETECTION SUMMARY - Run {args.run_id}")
    print("="*70)
    print(f"Total pathogens detected: {merged_results['summary']['total_pathogens_detected']}")
    print(f"PMDA pathogens detected: {merged_results['summary']['pmda_pathogens_detected']}")
    print(f"Detection methods used: {', '.join(merged_results['summary']['methods_used'])}")
    print(f"Critical findings: {merged_results['summary']['critical_findings_count']}")

    if merged_results['summary'].get('perv_detected'):
        print("\n[CRITICAL] PERV DETECTED")

    if alerts:
        print(f"\n[WARNING] {len(alerts)} ALERT(S) GENERATED:")
        for alert in alerts:
            print(f"  [{alert['priority']}] {alert['alert_type']}: {alert.get('message', alert.get('pathogen', ''))}")

    print(f"\nResults written to: {args.output}")
    print("="*70)

    # Exit code based on findings
    if merged_results['summary'].get('perv_detected'):
        sys.exit(2)  # PERV detected
    elif merged_results['summary']['critical_findings_count'] > 0:
        sys.exit(1)  # Critical pathogen detected
    else:
        sys.exit(0)  # Success


if __name__ == '__main__':
    main()
