#!/usr/bin/env python3
"""
Duplicate Circular Genome References for Junction Read Mapping

Purpose: Create duplicated FASTA references for circular genomes (PCV2, PCV3, TTV)
         to properly handle sequencing reads spanning the circularization junction.

Version: 1.0 (2025-11-13)
Related: Protocol 12 v2.1, Circular Genome Handling Guide

Example:
    Original circular genome (1768 bp):
        Position 1 ────────────── Position 1768
             ↑                          ↓
             └──────────────────────────┘

    Duplicated reference (3536 bp):
        Position 1 ──→ Position 1768 ──→ Position 1769 ──→ Position 3536
        └─── First copy ───┘  └─── Second copy (duplicate) ───┘

Usage:
    python duplicate_circular_genomes.py \\
        --input PCV2.fasta PCV3.fasta TTV.fasta \\
        --output circular_genomes_duplicated.fasta \\
        --metadata circular_genomes_metadata.json
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Known circular genomes requiring duplication
CIRCULAR_GENOMES = {
    'PCV2': {
        'full_name': 'Porcine Circovirus 2',
        'expected_length': 1768,
        'genbank': 'AF027217.1',
        'topology': 'circular ssDNA'
    },
    'PCV3': {
        'full_name': 'Porcine Circovirus 3',
        'expected_length': 2000,
        'genbank': 'MF318988.1',
        'topology': 'circular ssDNA'
    },
    'TTV': {
        'full_name': 'Torque Teno Virus',
        'expected_length': 3800,
        'genbank': 'AB017610.1',
        'topology': 'circular ssDNA'
    }
}


def parse_fasta(fasta_file: Path) -> List[SeqRecord]:
    """Parse FASTA file and return list of SeqRecords."""
    records = []
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            records.append(record)
    return records


def identify_circular_genome(record: SeqRecord) -> Tuple[str, Dict]:
    """
    Identify if a sequence record is a known circular genome.

    Args:
        record: BioPython SeqRecord

    Returns:
        Tuple of (genome_code, metadata_dict) or (None, None) if not circular
    """
    record_id = record.id.upper()
    record_desc = record.description.upper()
    record_len = len(record.seq)

    for code, metadata in CIRCULAR_GENOMES.items():
        # Check by GenBank accession
        if metadata['genbank'] in record.id:
            return code, metadata

        # Check by name
        if code in record_id or code in record_desc:
            return code, metadata

        # Check by full name
        if metadata['full_name'].upper() in record_desc:
            return code, metadata

        # Check by length (±10% tolerance)
        expected_len = metadata['expected_length']
        if abs(record_len - expected_len) / expected_len < 0.10:
            # Length matches, likely the same genome
            print(f"  Identified {record.id} as {code} by length ({record_len} bp)")
            return code, metadata

    return None, None


def duplicate_circular_sequence(record: SeqRecord, genome_code: str) -> SeqRecord:
    """
    Duplicate a circular genome sequence.

    Args:
        record: Original SeqRecord
        genome_code: Genome code (e.g., 'PCV2')

    Returns:
        New SeqRecord with duplicated sequence
    """
    original_seq = str(record.seq)
    original_len = len(original_seq)
    duplicated_seq = original_seq + original_seq
    duplicated_len = len(duplicated_seq)

    # Create new SeqRecord
    new_record = SeqRecord(
        Seq(duplicated_seq),
        id=f"{genome_code}_circular_dup",
        name=f"{genome_code}_circular_dup",
        description=(
            f"{record.description} | Duplicated for junction read mapping | "
            f"Original length: {original_len} bp | Duplicated length: {duplicated_len} bp | "
            f"Topology: circular | Protocol 12 v2.1 compatible"
        )
    )

    return new_record


def main():
    parser = argparse.ArgumentParser(
        description='Duplicate circular genome references for junction read mapping',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Duplicate single genome
  python duplicate_circular_genomes.py --input PCV2.fasta --output PCV2_dup.fasta

  # Duplicate multiple genomes
  python duplicate_circular_genomes.py \\
      --input PCV2.fasta PCV3.fasta TTV.fasta \\
      --output circular_genomes_duplicated.fasta \\
      --metadata circular_metadata.json

Output:
  - Duplicated FASTA file ready for Minimap2 indexing
  - Metadata JSON with original and duplicated lengths
        """
    )
    parser.add_argument('--input', nargs='+', required=True, type=Path,
                       help='Input FASTA files (one or more)')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output FASTA file with duplicated sequences')
    parser.add_argument('--metadata', type=Path,
                       help='Output metadata JSON file (optional)')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite output files if they exist')

    args = parser.parse_args()

    # Check if output exists
    if args.output.exists() and not args.force:
        print(f"ERROR: Output file {args.output} already exists. Use --force to overwrite.")
        return 1

    # Check if BioPython is installed
    try:
        from Bio import SeqIO
    except ImportError:
        print("ERROR: BioPython not installed. Install with: pip install biopython")
        return 1

    print("=== Circular Genome Reference Duplication ===\n")
    print(f"Input files: {len(args.input)}")
    print(f"Output: {args.output}\n")

    duplicated_records = []
    metadata = {
        'version': '1.0',
        'date': '2025-11-13',
        'protocol': 'Protocol 12 v2.1',
        'genomes': {}
    }

    # Process each input file
    for input_file in args.input:
        if not input_file.exists():
            print(f"WARNING: Input file not found: {input_file}")
            continue

        print(f"Processing {input_file.name}...")
        records = parse_fasta(input_file)
        print(f"  Found {len(records)} sequence(s)")

        for record in records:
            # Identify if circular genome
            genome_code, genome_metadata = identify_circular_genome(record)

            if genome_code:
                print(f"  ✓ Identified: {genome_code} ({genome_metadata['full_name']})")
                print(f"    Original length: {len(record.seq)} bp")
                print(f"    Topology: {genome_metadata['topology']}")

                # Duplicate sequence
                duplicated_record = duplicate_circular_sequence(record, genome_code)
                duplicated_records.append(duplicated_record)

                print(f"    Duplicated length: {len(duplicated_record.seq)} bp")
                print(f"    New ID: {duplicated_record.id}\n")

                # Store metadata
                metadata['genomes'][genome_code] = {
                    'full_name': genome_metadata['full_name'],
                    'original_id': record.id,
                    'original_length': len(record.seq),
                    'duplicated_length': len(duplicated_record.seq),
                    'topology': genome_metadata['topology'],
                    'genbank': genome_metadata['genbank'],
                    'duplicated_id': duplicated_record.id
                }
            else:
                print(f"  ⚠ Warning: {record.id} not recognized as circular genome")
                print(f"    Length: {len(record.seq)} bp")
                print(f"    Skipping (not duplicated)\n")

    # Write duplicated sequences
    if duplicated_records:
        with open(args.output, 'w') as f:
            SeqIO.write(duplicated_records, f, 'fasta')

        print(f"\n=== Summary ===")
        print(f"Duplicated genomes: {len(duplicated_records)}")
        print(f"Output FASTA: {args.output}")
        print(f"Total sequences: {len(duplicated_records)}")

        # Write metadata if requested
        if args.metadata:
            with open(args.metadata, 'w') as f:
                json.dump(metadata, f, indent=2)
            print(f"Metadata JSON: {args.metadata}")

        print("\nNext steps:")
        print("1. Combine with other PMDA pathogen references:")
        print(f"   cat {args.output} other_pmda_references.fasta > pmda_all_91.fasta")
        print("2. Build Minimap2 index:")
        print("   minimap2 -d pmda_all_91.mmi pmda_all_91.fasta")
        print("3. Update pipeline configuration to use duplicated references")

        return 0
    else:
        print("\nERROR: No circular genomes found in input files")
        print("Expected genomes: PCV2, PCV3, TTV")
        print("Check input FASTA headers match GenBank accessions or contain genome names")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
