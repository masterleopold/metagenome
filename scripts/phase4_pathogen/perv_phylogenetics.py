#!/usr/bin/env python3
"""
PERV phylogenetic analysis
Determines phylogenetic relationship of detected PERV sequences
"""

import argparse
import json
from pathlib import Path
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess

def run_phylogenetic_analysis(consensus_file: Path, reference_file: Path) -> dict:
    """Perform phylogenetic analysis of PERV sequences."""

    # Combine consensus and reference sequences
    sequences = []
    for record in SeqIO.parse(consensus_file, "fasta"):
        sequences.append(record)
    for record in SeqIO.parse(reference_file, "fasta"):
        sequences.append(record)

    # Write combined file
    combined_file = consensus_file.parent / "combined_sequences.fasta"
    SeqIO.write(sequences, combined_file, "fasta")

    # Run MAFFT alignment
    aligned_file = consensus_file.parent / "aligned.fasta"
    subprocess.run([
        "mafft", "--auto", str(combined_file)
    ], stdout=open(aligned_file, 'w'), stderr=subprocess.DEVNULL)

    # Calculate distance matrix
    alignment = AlignIO.read(aligned_file, "fasta")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Build tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Find closest reference
    results = {
        'closest_reference': None,
        'distance_to_reference': 999,
        'phylogenetic_clade': 'Unknown'
    }

    # Simple distance-based classification
    for i, seq_id in enumerate(distance_matrix.names):
        if 'consensus' in seq_id.lower():
            min_dist = 999
            closest_ref = None
            for j, ref_id in enumerate(distance_matrix.names):
                if 'reference' in ref_id.lower() and i != j:
                    dist = distance_matrix[i][j]
                    if dist < min_dist:
                        min_dist = dist
                        closest_ref = ref_id

            results['closest_reference'] = closest_ref
            results['distance_to_reference'] = round(min_dist, 4)

            # Determine clade based on distance
            if min_dist < 0.05:
                results['phylogenetic_clade'] = 'Highly similar to known strain'
            elif min_dist < 0.15:
                results['phylogenetic_clade'] = 'Moderately divergent'
            else:
                results['phylogenetic_clade'] = 'Highly divergent/Novel'

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--consensus', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--run-id', required=True)

    args = parser.parse_args()

    if args.consensus.exists():
        results = run_phylogenetic_analysis(args.consensus, args.reference)
    else:
        results = {
            'error': 'No consensus sequence available',
            'phylogenetic_clade': 'Not determined'
        }

    results['run_id'] = args.run_id

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == '__main__':
    main()