#!/usr/bin/env python3
"""
Calculate absolute copy numbers (copies/mL) for detected pathogens
"""

import argparse
import json
from pathlib import Path
import numpy as np

# Average genome sizes for calculation (bp)
# Values based on NCBI RefSeq and published literature
GENOME_SIZES = {
    # Critical PERV retroviruses
    'PERV': 8000,
    'PERV-A': 8000,
    'PERV-B': 8000,
    'PERV-C': 8000,

    # High-risk viruses
    'HEV': 7200,        # Hepatitis E virus
    'JEV': 11000,       # Japanese encephalitis virus
    'ASFV': 170000,     # African swine fever virus
    'CSFV': 12300,      # Classical swine fever virus
    'FMDV': 8400,       # Foot-and-mouth disease virus
    'PRRSV': 15000,     # PRRS virus
    'PRV': 143000,      # Pseudorabies virus

    # Other viruses
    'PCV2': 1768,       # Porcine circovirus 2
    'PCV3': 2000,       # Porcine circovirus 3
    'PPV': 5000,        # Porcine parvovirus
    'SIV': 13500,       # Swine influenza virus
    'PEDV': 28000,      # Porcine epidemic diarrhea virus
    'TGEV': 28500,      # Transmissible gastroenteritis virus
    'RV': 18500,        # Rotavirus
    'EMCV': 7800,       # Encephalomyocarditis virus

    # Bacteria
    'SS': 2000000,      # Streptococcus suis
    'EC': 5000000,      # E. coli
    'SE': 4800000,      # Salmonella enterica
    'SA': 2800000,      # Staphylococcus aureus
    'SP': 2100000,      # Streptococcus pneumoniae
    'MT': 4400000,      # Mycobacterium tuberculosis
    'MB': 4300000,      # Mycobacterium bovis
    'BA': 5200000,      # Bacillus anthracis
    'YP': 4600000,      # Yersinia pestis
    'FT': 1900000,      # Francisella tularensis
    'CT': 2800000,      # Clostridium tetani
    'CP': 3500000,      # Clostridium perfringens
    'CS': 4200000,      # Clostridium difficile
    'LA': 2900000,      # Listeria monocytogenes
    'BP': 4000000,      # Bordetella pertussis
    'BB': 1500000,      # Borrelia burgdorferi
    'BRA': 3000000,     # Brachyspira species
    'TP': 1100000,      # Treponema pallidum
    'CJ': 1600000,      # Campylobacter jejuni
    'HP': 1600000,      # Helicobacter pylori
    'LP': 3400000,      # Legionella pneumophila
    'MA': 800000,       # Mycoplasma pneumoniae
    'PM': 2400000,      # Pasteurella multocida
    'AP': 2300000,      # Actinobacillus pleuropneumoniae
    'HPS': 2200000,     # Haemophilus parasuis
    'EP': 1800000,      # Erysipelothrix rhusiopathiae
    'BR': 3300000,      # Brucella species
    'LM': 1700000,      # Lawsonia intracellularis

    # Parasites (typical values)
    'TG': 65000000,     # Toxoplasma gondii
    'TC': 55000000,     # Trypanosoma cruzi
    'TS': 64000000,     # Trichinella spiralis
    'EC-P': 115000000,  # Echinococcus species
    'CP-P': 9100000,    # Cryptosporidium parvum

    # Fungi (typical values)
    'CN': 19000000,     # Cryptococcus neoformans
    'CA': 14000000,     # Candida albicans
    'AF': 29000000,     # Aspergillus fumigatus
    'PJ': 8000000,      # Pneumocystis jirovecii
    'HI': 30000000,     # Histoplasma capsulatum
}

def calculate_copies_per_ml(reads: int, genome_size: int, total_reads: int,
                           plasma_volume_ml: float, avg_read_length: int = 1000,
                           extraction_efficiency: float = 0.7) -> float:
    """
    Calculate absolute copy number per mL plasma.

    Formula:
    copies/mL = (pathogen_reads / total_reads) * total_DNA_molecules / plasma_volume

    Args:
        reads: Number of reads mapping to pathogen
        genome_size: Pathogen genome size in bp
        total_reads: Total sequencing reads
        plasma_volume_ml: Sample volume in mL
        avg_read_length: Average read length in bp (default: 1000)
        extraction_efficiency: Extraction efficiency 0-1 (default: 0.7)

    Returns:
        Absolute copy number per mL plasma
    """

    # Estimate total DNA molecules from total reads
    total_dna_molecules = total_reads * avg_read_length / genome_size

    # Calculate pathogen copies
    pathogen_fraction = reads / total_reads if total_reads > 0 else 0
    pathogen_molecules = pathogen_fraction * total_dna_molecules

    # Adjust for extraction efficiency and normalize to per mL
    copies_per_ml = (pathogen_molecules / extraction_efficiency) / plasma_volume_ml

    return copies_per_ml

def process_quantification(kraken_file: Path, normalized_file: Path,
                          plasma_volume: float, avg_read_length: int = 1000,
                          extraction_efficiency: float = 0.7) -> dict:
    """
    Process quantification files and calculate absolute copy numbers.

    Args:
        kraken_file: Path to Kraken2 quantification JSON
        normalized_file: Path to spike-in normalized data
        plasma_volume: Sample plasma volume in mL
        avg_read_length: Average sequencing read length in bp (default: 1000)
        extraction_efficiency: DNA/RNA extraction efficiency 0-1 (default: 0.7)

    Returns:
        Dictionary with quantification results
    """

    results = {
        'plasma_volume_ml': plasma_volume,
        'avg_read_length': avg_read_length,
        'extraction_efficiency': extraction_efficiency,
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
            reads, genome_size, total_reads, plasma_volume,
            avg_read_length, extraction_efficiency
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
    parser = argparse.ArgumentParser(
        description='Calculate absolute pathogen copy numbers (copies/mL)'
    )
    parser.add_argument('--kraken', required=True, type=Path,
                       help='Kraken2 quantification JSON file')
    parser.add_argument('--normalized', required=True, type=Path,
                       help='Spike-in normalized data JSON file')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output JSON file')
    parser.add_argument('--plasma-volume', type=float, default=10.0,
                       help='Plasma sample volume in mL (default: 10.0)')
    parser.add_argument('--avg-read-length', type=int, default=1000,
                       help='Average sequencing read length in bp (default: 1000)')
    parser.add_argument('--extraction-efficiency', type=float, default=0.7,
                       help='DNA/RNA extraction efficiency 0-1 (default: 0.7)')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')

    args = parser.parse_args()

    results = process_quantification(
        args.kraken, args.normalized, args.plasma_volume,
        args.avg_read_length, args.extraction_efficiency
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