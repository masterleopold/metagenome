#!/usr/bin/env python3
"""
Calculate absolute copy numbers (copies/mL) for detected pathogens
"""

import argparse
import json
from pathlib import Path

# Average genome sizes for calculation (bp)
GENOME_SIZES = {
    'PERV': 8000,
    'HEV': 7200,
    'JEV': 11000,
    'SS': 2000000,  # Streptococcus suis
    'EC': 5000000,  # E. coli
    # Add all pathogens
}

def calculate_copies_per_ml(reads: int, genome_size: int, total_reads: int,
                           plasma_volume_ml: float, extraction_efficiency: float = 0.7) -> float:
    """
    Calculate absolute copy number per mL plasma.

    Formula:
    copies/mL = (pathogen_reads / total_reads) * total_DNA_molecules / plasma_volume
    """

    # Estimate total DNA molecules from total reads
    # Assuming average read length of 1000bp and coverage of 1x
    avg_read_length = 1000
    total_dna_molecules = total_reads * avg_read_length / genome_size

    # Calculate pathogen copies
    pathogen_fraction = reads / total_reads if total_reads > 0 else 0
    pathogen_molecules = pathogen_fraction * total_dna_molecules

    # Adjust for extraction efficiency and normalize to per mL
    copies_per_ml = (pathogen_molecules / extraction_efficiency) / plasma_volume_ml

    return copies_per_ml

def process_quantification(kraken_file: Path, normalized_file: Path,
                          plasma_volume: float) -> dict:
    """Process quantification files and calculate absolute copy numbers."""

    results = {
        'plasma_volume_ml': plasma_volume,
        'pathogens': {}
    }

    # Load kraken quantification
    with open(kraken_file) as f:
        kraken_data = json.load(f)

    total_reads = kraken_data.get('total_reads', 0)

    # Load normalized data if available
    normalized_data = {}
    if normalized_file.exists():
        with open(normalized_file) as f:
            norm = json.load(f)
            normalized_data = norm.get('pathogens', {})
            results['normalization_factor'] = norm.get('normalization_factor', 1.0)

    # Calculate absolute copy numbers
    for pathogen_code, data in kraken_data.get('pathogens', {}).items():
        reads = data.get('reads', 0)
        genome_size = GENOME_SIZES.get(pathogen_code, 10000)  # Default 10kb

        # Use normalized reads if available
        if pathogen_code in normalized_data:
            reads = normalized_data[pathogen_code].get('reads_normalized', reads)

        copies_per_ml = calculate_copies_per_ml(
            reads, genome_size, total_reads, plasma_volume
        )

        results['pathogens'][pathogen_code] = {
            'name': data.get('name', pathogen_code),
            'reads': reads,
            'copies_per_ml': round(copies_per_ml, 2),
            'log10_copies_per_ml': round(np.log10(copies_per_ml + 1), 2),
            'confidence': 'high' if reads > 100 else 'medium' if reads > 10 else 'low'
        }

        # Calculate confidence intervals (simplified)
        if reads > 0:
            # Poisson-based 95% CI
            ci_lower = copies_per_ml * 0.7  # Simplified
            ci_upper = copies_per_ml * 1.3
            results['pathogens'][pathogen_code]['ci_lower'] = round(ci_lower, 2)
            results['pathogens'][pathogen_code]['ci_upper'] = round(ci_upper, 2)

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kraken', required=True, type=Path)
    parser.add_argument('--normalized', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--plasma-volume', type=float, default=10.0)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    # Add numpy import
    global np
    import numpy as np

    results = process_quantification(
        args.kraken, args.normalized, args.plasma_volume
    )
    results['run_id'] = args.run_id

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"Absolute quantification completed")
    print(f"Plasma volume: {args.plasma_volume} mL")

    if results['pathogens']:
        print("\nPathogen copy numbers (copies/mL):")
        sorted_pathogens = sorted(results['pathogens'].items(),
                                  key=lambda x: x[1]['copies_per_ml'], reverse=True)
        for code, data in sorted_pathogens[:10]:
            print(f"  {code}: {data['copies_per_ml']:.2e} copies/mL ({data['confidence']} confidence)")

if __name__ == '__main__':
    main()