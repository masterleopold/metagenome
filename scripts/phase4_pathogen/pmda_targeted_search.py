#!/usr/bin/env python3
"""
PMDA Targeted Pathogen Search
Performs targeted alignment-based search for PMDA 91 pathogens
Uses Minimap2 for fast alignment against PMDA reference database
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from collections import defaultdict
import tempfile


def load_pmda_config(config_file: Path) -> dict:
    """Load PMDA pathogen configuration."""
    with open(config_file) as f:
        return json.load(f)


def run_minimap2_alignment(input_fastq: Path, database: Path,
                          threads: int, output_bam: Path) -> bool:
    """
    Run Minimap2 alignment against PMDA pathogen database.

    Returns True if successful.
    """
    print(f"Running Minimap2 alignment...")
    print(f"  Database: {database}")
    print(f"  Input: {input_fastq}")
    print(f"  Output: {output_bam}")

    try:
        # Run minimap2 with alignment parameters optimized for pathogen detection
        cmd = [
            'minimap2',
            '-ax', 'map-ont',  # Oxford Nanopore preset
            '-t', str(threads),
            '--secondary=no',  # No secondary alignments
            '-N', '10',  # Keep up to 10 alignments per read
            str(database),
            str(input_fastq)
        ]

        # Pipe to samtools for BAM conversion and sorting
        minimap_proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        samtools_cmd = [
            'samtools', 'sort',
            '-@', str(threads),
            '-o', str(output_bam),
            '-'
        ]

        samtools_proc = subprocess.Popen(
            samtools_cmd,
            stdin=minimap_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        minimap_proc.stdout.close()

        # Wait for completion
        samtools_stdout, samtools_stderr = samtools_proc.communicate()
        minimap_stderr = minimap_proc.stderr.read()

        if samtools_proc.returncode != 0:
            print(f"ERROR: Samtools failed: {samtools_stderr.decode()}", file=sys.stderr)
            return False

        if minimap_proc.returncode != 0:
            print(f"ERROR: Minimap2 failed: {minimap_stderr.decode()}", file=sys.stderr)
            return False

        # Index BAM file
        subprocess.run(['samtools', 'index', str(output_bam)], check=True)

        return True

    except Exception as e:
        print(f"ERROR: Alignment failed: {e}", file=sys.stderr)
        return False


def parse_bam_results(bam_file: Path, pmda_config: dict,
                     min_identity: float = 0.90,
                     min_length: int = 100) -> dict:
    """
    Parse BAM file and extract PMDA pathogen matches.
    """
    try:
        import pysam
    except ImportError:
        print("ERROR: pysam not installed. Install with: pip install pysam", file=sys.stderr)
        sys.exit(1)

    print("Parsing alignment results...")

    # Build reference name to pathogen mapping
    ref_to_pathogen = {}
    for category, data in pmda_config['categories'].items():
        for pathogen in data.get('pathogens', []):
            code = pathogen['code']
            name = pathogen['name']
            # Reference sequences might be named by code or name
            ref_to_pathogen[code] = pathogen
            ref_to_pathogen[name.lower()] = pathogen
            # Also try variations
            ref_to_pathogen[name.lower().replace(' ', '_')] = pathogen

    detections = defaultdict(lambda: {
        'reads': 0,
        'total_bases': 0,
        'mapped_positions': set(),
        'identity_scores': [],
        'read_names': []
    })

    bam = pysam.AlignmentFile(str(bam_file), "rb")

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Calculate alignment identity
        aligned_length = read.query_alignment_length
        if aligned_length < min_length:
            continue

        # Get identity from NM tag (edit distance)
        try:
            nm = read.get_tag('NM')
            identity = 1.0 - (nm / aligned_length)
        except KeyError:
            # Fallback: estimate from CIGAR
            matches = sum(count for op, count in read.cigartuples if op == 0)  # M operations
            identity = matches / aligned_length if aligned_length > 0 else 0

        if identity < min_identity:
            continue

        # Map reference to pathogen
        ref_name = read.reference_name
        pathogen = None

        # Try exact match first
        if ref_name in ref_to_pathogen:
            pathogen = ref_to_pathogen[ref_name]
        else:
            # Try partial matching
            ref_lower = ref_name.lower()
            for ref_key, path_data in ref_to_pathogen.items():
                if ref_key.lower() in ref_lower or ref_lower in ref_key.lower():
                    pathogen = path_data
                    break

        if not pathogen:
            continue

        code = pathogen['code']

        # Record detection
        detections[code]['reads'] += 1
        detections[code]['total_bases'] += aligned_length
        detections[code]['identity_scores'].append(identity)
        detections[code]['read_names'].append(read.query_name)

        # Track coverage positions
        for pos in range(read.reference_start, read.reference_end):
            detections[code]['mapped_positions'].add(pos)

        # Store pathogen metadata (first time only)
        if 'name' not in detections[code]:
            detections[code]['name'] = pathogen['name']
            detections[code]['category'] = None
            detections[code]['risk_level'] = pathogen['risk_level']

            # Find category
            for category, data in pmda_config['categories'].items():
                if pathogen in data.get('pathogens', []):
                    detections[code]['category'] = category
                    break

    bam.close()

    # Process detections
    results = {}
    for code, data in detections.items():
        if data['reads'] == 0:
            continue

        # Calculate statistics
        avg_identity = sum(data['identity_scores']) / len(data['identity_scores'])
        coverage_breadth = len(data['mapped_positions'])

        results[code] = {
            'code': code,
            'name': data['name'],
            'category': data['category'],
            'risk_level': data['risk_level'],
            'reads': data['reads'],
            'average_identity': round(avg_identity, 4),
            'total_bases_aligned': data['total_bases'],
            'coverage_breadth': coverage_breadth,
            'unique_reads': len(set(data['read_names'])),
            'detection_confidence': 'high' if data['reads'] >= 100 else 'medium' if data['reads'] >= 10 else 'low'
        }

    return results


def main():
    parser = argparse.ArgumentParser(
        description='PMDA targeted pathogen search using Minimap2'
    )
    parser.add_argument('--input', required=True, type=Path,
                       help='Input FASTQ file or directory')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output directory')
    parser.add_argument('--database', required=True, type=Path,
                       help='PMDA pathogen reference database (FASTA)')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('--config', type=Path,
                       default='/home/user/metagenome/templates/config/pmda_pathogens.json',
                       help='PMDA pathogen configuration file')
    parser.add_argument('--threads', type=int, default=16,
                       help='Number of threads (default: 16)')
    parser.add_argument('--min-identity', type=float, default=0.90,
                       help='Minimum alignment identity (default: 0.90)')
    parser.add_argument('--min-length', type=int, default=100,
                       help='Minimum alignment length (default: 100)')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("PMDA Targeted Pathogen Search")
    print("="*70)
    print(f"Run ID: {args.run_id}")
    print(f"Database: {args.database}")
    print(f"Min identity: {args.min_identity}")
    print(f"Min length: {args.min_length}")
    print(f"Threads: {args.threads}")
    print("="*70)

    # Validate inputs
    if not args.database.exists():
        print(f"ERROR: Database not found: {args.database}", file=sys.stderr)
        sys.exit(1)

    # Try multiple locations for config file
    config_file = args.config
    if not config_file.exists():
        script_dir = Path(__file__).parent.parent.parent
        config_file = script_dir / 'templates' / 'config' / 'pmda_pathogens.json'

    if not config_file.exists():
        print(f"ERROR: PMDA config not found: {config_file}", file=sys.stderr)
        sys.exit(1)

    # Load PMDA configuration
    pmda_config = load_pmda_config(config_file)
    print(f"Loaded PMDA pathogen database (version {pmda_config.get('version', 'unknown')})")

    # Find input FASTQ
    if args.input.is_dir():
        # Find first FASTQ file in directory
        fastq_files = list(args.input.glob('*.fastq')) + list(args.input.glob('*.fastq.gz'))
        if not fastq_files:
            print(f"ERROR: No FASTQ files found in {args.input}", file=sys.stderr)
            sys.exit(1)
        input_fastq = fastq_files[0]
    else:
        input_fastq = args.input

    if not input_fastq.exists():
        print(f"ERROR: Input file not found: {input_fastq}", file=sys.stderr)
        sys.exit(1)

    # Run alignment
    output_bam = args.output / 'pmda_alignments.bam'

    success = run_minimap2_alignment(
        input_fastq, args.database, args.threads, output_bam
    )

    if not success:
        print("ERROR: Alignment failed", file=sys.stderr)
        sys.exit(1)

    # Parse results
    detections = parse_bam_results(
        output_bam, pmda_config, args.min_identity, args.min_length
    )

    # Calculate summary statistics
    critical_detections = [d for d in detections.values() if d['risk_level'] == 'CRITICAL']
    high_risk_detections = [d for d in detections.values() if d['risk_level'] == 'HIGH']

    summary = {
        'run_id': args.run_id,
        'method': 'pmda_targeted_minimap2',
        'database': str(args.database),
        'total_pathogens_detected': len(detections),
        'critical_pathogens_detected': len(critical_detections),
        'high_risk_pathogens_detected': len(high_risk_detections),
        'pathogens': detections,
        'parameters': {
            'min_identity': args.min_identity,
            'min_length': args.min_length,
            'threads': args.threads
        }
    }

    # Write output
    output_file = args.output / 'pmda_targeted_results.json'
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print("\n" + "="*70)
    print("PMDA Targeted Search Summary")
    print("="*70)
    print(f"Total PMDA pathogens detected: {len(detections)}")
    print(f"Critical pathogens detected: {len(critical_detections)}")
    print(f"High-risk pathogens detected: {len(high_risk_detections)}")

    if critical_detections:
        print("\n[WARNING] CRITICAL PATHOGENS DETECTED:")
        for det in critical_detections:
            print(f"  - {det['code']}: {det['name']} ({det['reads']} reads, "
                  f"{det['average_identity']:.2%} identity)")

    print(f"\nResults written to: {output_file}")
    print("="*70)

    # Exit code
    if len(critical_detections) > 0:
        sys.exit(2)  # Critical pathogen detected
    elif len(detections) > 0:
        sys.exit(0)  # Normal detection
    else:
        sys.exit(0)  # No detections


if __name__ == '__main__':
    main()
