#!/usr/bin/env python3
"""
Normalize pathogen quantities using spike-in control
"""

import argparse
import json
from pathlib import Path

def calculate_spike_in_recovery(report_file: Path, spike_in: str) -> float:
    """Calculate spike-in recovery rate."""

    spike_in_reads = 0
    total_reads = 0

    with open(report_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            reads = int(parts[1])
            name = parts[5].strip()

            if spike_in.lower() in name.lower():
                spike_in_reads += reads

            if parts[3] == 'U' or parts[3] == 'R':  # Unclassified or Root
                total_reads += reads

    # Expected spike-in percentage (e.g., 1%)
    expected_percentage = 1.0
    observed_percentage = (spike_in_reads / total_reads * 100) if total_reads > 0 else 0

    recovery_rate = observed_percentage / expected_percentage if expected_percentage > 0 else 1.0

    return recovery_rate

def normalize_quantities(quantification: dict, recovery_rate: float) -> dict:
    """Normalize pathogen quantities based on spike-in recovery."""

    normalized = {
        'recovery_rate': round(recovery_rate, 4),
        'normalization_factor': round(1.0 / recovery_rate, 4) if recovery_rate > 0 else 1.0,
        'pathogens': {}
    }

    for pathogen, data in quantification.get('pathogens', {}).items():
        normalized['pathogens'][pathogen] = {
            'reads_normalized': round(data['reads'] / recovery_rate) if recovery_rate > 0 else data['reads'],
            'rpm_normalized': round(data['rpm'] / recovery_rate, 2) if recovery_rate > 0 else data['rpm']
        }

    return normalized

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', required=True, type=Path)
    parser.add_argument('--spike-in', default='PhiX174')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    recovery_rate = calculate_spike_in_recovery(args.report, args.spike_in)

    # Load quantification data
    quant_file = args.report.parent / 'kraken_quantification.json'
    if quant_file.exists():
        with open(quant_file) as f:
            quantification = json.load(f)
    else:
        quantification = {'pathogens': {}}

    results = normalize_quantities(quantification, recovery_rate)
    results['run_id'] = args.run_id
    results['spike_in'] = args.spike_in

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Spike-in recovery rate: {recovery_rate:.2%}")
    print(f"Normalization factor: {results['normalization_factor']:.2f}x")

if __name__ == '__main__':
    main()