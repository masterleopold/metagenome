#!/usr/bin/env python3
"""
PMDA 4-Virus Specific Detection
Implements high-sensitivity detection for Polyomavirus, Hantavirus, EEEV, and Spumavirus

Handles virus-specific requirements:
- Polyomavirus: CpG-depleted DNA, coverage-based detection
- Hantavirus: 3-segment concordance (L AND M AND S)
- EEEV: Poly(A)+ RNA, alphavirus phylogeny
- Spumavirus: Nested PCR from PBMC DNA (NOT metagenomic)
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
import tempfile

# Virus-specific configuration
VIRUS_CONFIG = {
    "polyomavirus": {
        "database": "polyomavirus/minimap2/polyoma_all.mmi",
        "min_reads": 100,
        "min_coverage": 10,  # 10× mean depth
        "min_identity": 0.90,
        "genome_size": 5200,  # approximate
        "sample_type": "plasma_cfDNA",
        "detection_method": "minimap2_alignment"
    },
    "hantavirus": {
        "database": "hantavirus/minimap2/hantavirus_all.mmi",
        "segments": {
            "L": {"database": "hantavirus/minimap2/hantavirus_L.mmi", "min_reads": 50, "size": 6533},
            "M": {"database": "hantavirus/minimap2/hantavirus_M.mmi", "min_reads": 50, "size": 3651},
            "S": {"database": "hantavirus/minimap2/hantavirus_S.mmi", "min_reads": 50, "size": 1696}
        },
        "min_identity": 0.85,
        "sample_type": "plasma_cfRNA",
        "detection_method": "segment_concordance",
        "requires_all_segments": True
    },
    "eeev": {
        "database": "alphavirus/minimap2/alphavirus_all.mmi",
        "min_reads": 100,
        "min_coverage": 10,
        "min_identity": 0.90,
        "genome_size": 11841,
        "sample_type": "plasma_cfRNA",
        "detection_method": "minimap2_alignment",
        "phylogeny_required": True
    },
    "spumavirus": {
        "database": "spumavirus/blast/spumavirus_pol",
        "sample_type": "PBMC_genomic_DNA",
        "detection_method": "nested_pcr",
        "note": "Not detected by metagenomic sequencing - requires separate nested PCR workflow"
    }
}


def run_minimap2_alignment(fastq: Path, database: Path, threads: int, output_bam: Path) -> bool:
    """Run Minimap2 alignment against virus-specific database."""
    try:
        print(f"Running Minimap2 alignment...")
        print(f"  Database: {database}")
        print(f"  Input: {fastq}")

        # Minimap2 alignment
        minimap_cmd = [
            'minimap2',
            '-ax', 'map-ont',
            '-t', str(threads),
            '--secondary=no',
            '-N', '10',
            str(database),
            str(fastq)
        ]

        # Pipe to samtools
        samtools_cmd = [
            'samtools', 'sort',
            '-@', str(threads),
            '-o', str(output_bam),
            '-'
        ]

        minimap_proc = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        samtools_proc = subprocess.Popen(samtools_cmd, stdin=minimap_proc.stdout,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        minimap_proc.stdout.close()
        samtools_stdout, samtools_stderr = samtools_proc.communicate()

        if samtools_proc.returncode != 0:
            print(f"ERROR: Alignment failed: {samtools_stderr.decode()}", file=sys.stderr)
            return False

        # Index BAM
        subprocess.run(['samtools', 'index', str(output_bam)], check=True)

        return True

    except Exception as e:
        print(f"ERROR: Alignment failed: {e}", file=sys.stderr)
        return False


def detect_polyomavirus(fastq: Path, db_path: Path, threads: int, output_dir: Path) -> Dict:
    """
    Detect polyomavirus with coverage-based validation.

    Criteria:
    - ≥100 reads mapped
    - ≥10× mean coverage
    - ≥90% identity
    """
    print("\n" + "="*70)
    print("POLYOMAVIRUS DETECTION")
    print("="*70)

    config = VIRUS_CONFIG["polyomavirus"]
    database = db_path / config["database"]

    if not database.exists():
        return {"status": "ERROR", "message": f"Database not found: {database}"}

    # Run alignment
    bam_file = output_dir / "polyomavirus_alignments.bam"
    success = run_minimap2_alignment(fastq, database, threads, bam_file)

    if not success:
        return {"status": "ERROR", "message": "Alignment failed"}

    # Parse results
    try:
        import pysam

        bam = pysam.AlignmentFile(str(bam_file), "rb")

        stats = defaultdict(lambda: {
            'reads': 0,
            'bases_aligned': 0,
            'identity_scores': [],
            'coverage_positions': set()
        })

        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Calculate identity
            aligned_length = read.query_alignment_length
            if aligned_length < 100:
                continue

            try:
                nm = read.get_tag('NM')
                identity = 1.0 - (nm / aligned_length)
            except KeyError:
                matches = sum(count for op, count in read.cigartuples if op == 0)
                identity = matches / aligned_length if aligned_length > 0 else 0

            if identity < config["min_identity"]:
                continue

            ref_name = read.reference_name

            stats[ref_name]['reads'] += 1
            stats[ref_name]['bases_aligned'] += aligned_length
            stats[ref_name]['identity_scores'].append(identity)

            # Track coverage
            for pos in range(read.reference_start, read.reference_end):
                stats[ref_name]['coverage_positions'].add(pos)

        bam.close()

        # Evaluate results
        results = {}
        for ref_name, data in stats.items():
            if data['reads'] == 0:
                continue

            avg_identity = sum(data['identity_scores']) / len(data['identity_scores'])
            coverage_breadth = len(data['coverage_positions'])
            mean_depth = data['bases_aligned'] / config["genome_size"]

            # Detection criteria
            detected = (data['reads'] >= config["min_reads"] and
                       mean_depth >= config["min_coverage"])

            confidence = "HIGH" if (data['reads'] >= 200 and mean_depth >= 20) else \
                        "MEDIUM" if detected else "LOW"

            results[ref_name] = {
                'detected': detected,
                'reads': data['reads'],
                'average_identity': round(avg_identity, 4),
                'coverage_breadth': coverage_breadth,
                'mean_depth': round(mean_depth, 2),
                'confidence': confidence
            }

        # Determine final status
        any_detected = any(r['detected'] for r in results.values())

        return {
            'status': 'DETECTED' if any_detected else 'NOT_DETECTED',
            'virus': 'Polyomavirus',
            'references': results,
            'method': 'minimap2_alignment',
            'detection_criteria': f"≥{config['min_reads']} reads AND ≥{config['min_coverage']}× coverage"
        }

    except ImportError:
        return {"status": "ERROR", "message": "pysam not installed"}
    except Exception as e:
        return {"status": "ERROR", "message": str(e)}


def detect_hantavirus(fastq: Path, db_path: Path, threads: int, output_dir: Path) -> Dict:
    """
    Detect hantavirus with 3-segment concordance validation.

    Criteria:
    - L segment: ≥50 reads
    - M segment: ≥50 reads
    - S segment: ≥50 reads
    - All three segments MUST be detected simultaneously
    """
    print("\n" + "="*70)
    print("HANTAVIRUS DETECTION (3-Segment Concordance)")
    print("="*70)

    config = VIRUS_CONFIG["hantavirus"]
    segment_results = {}

    # Align to each segment database
    for segment, seg_config in config["segments"].items():
        print(f"\nAligning to {segment} segment...")

        database = db_path / seg_config["database"]
        if not database.exists():
            return {"status": "ERROR", "message": f"Database not found: {database}"}

        bam_file = output_dir / f"hantavirus_{segment}_alignments.bam"
        success = run_minimap2_alignment(fastq, database, threads, bam_file)

        if not success:
            return {"status": "ERROR", "message": f"Alignment failed for {segment} segment"}

        # Count reads for this segment
        try:
            import pysam

            bam = pysam.AlignmentFile(str(bam_file), "rb")

            segment_stats = {
                'reads': 0,
                'bases_aligned': 0,
                'species_reads': defaultdict(int),
                'identity_scores': []
            }

            for read in bam:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                aligned_length = read.query_alignment_length
                if aligned_length < 100:
                    continue

                try:
                    nm = read.get_tag('NM')
                    identity = 1.0 - (nm / aligned_length)
                except KeyError:
                    matches = sum(count for op, count in read.cigartuples if op == 0)
                    identity = matches / aligned_length if aligned_length > 0 else 0

                if identity < config["min_identity"]:
                    continue

                segment_stats['reads'] += 1
                segment_stats['bases_aligned'] += aligned_length
                segment_stats['identity_scores'].append(identity)

                # Track species (Hantaan, Seoul, Dobrava, Puumala)
                ref_name = read.reference_name.lower()
                if 'hantaan' in ref_name:
                    segment_stats['species_reads']['Hantaan'] += 1
                elif 'seoul' in ref_name:
                    segment_stats['species_reads']['Seoul'] += 1
                elif 'dobrava' in ref_name:
                    segment_stats['species_reads']['Dobrava'] += 1
                elif 'puumala' in ref_name:
                    segment_stats['species_reads']['Puumala'] += 1

            bam.close()

            # Calculate statistics
            avg_identity = (sum(segment_stats['identity_scores']) / len(segment_stats['identity_scores'])
                          if segment_stats['identity_scores'] else 0)
            mean_depth = segment_stats['bases_aligned'] / seg_config["size"]

            # Determine most likely species for this segment
            most_likely_species = max(segment_stats['species_reads'].items(),
                                     key=lambda x: x[1])[0] if segment_stats['species_reads'] else "Unknown"

            segment_results[segment] = {
                'reads': segment_stats['reads'],
                'detected': segment_stats['reads'] >= seg_config["min_reads"],
                'average_identity': round(avg_identity, 4),
                'mean_depth': round(mean_depth, 2),
                'most_likely_species': most_likely_species,
                'species_distribution': dict(segment_stats['species_reads'])
            }

        except ImportError:
            return {"status": "ERROR", "message": "pysam not installed"}
        except Exception as e:
            return {"status": "ERROR", "message": str(e)}

    # Check 3-segment concordance
    all_segments_detected = all(seg['detected'] for seg in segment_results.values())

    # Determine consensus species (most reads across all segments)
    species_total = defaultdict(int)
    for seg_data in segment_results.values():
        for species, count in seg_data['species_distribution'].items():
            species_total[species] += count

    consensus_species = max(species_total.items(), key=lambda x: x[1])[0] if species_total else "Unknown"

    # Final determination
    if all_segments_detected:
        status = "DETECTED"
        confidence = "HIGH"
        message = f"All 3 segments detected. Species: {consensus_species}"
    elif any(seg['detected'] for seg in segment_results.values()):
        status = "INCONCLUSIVE"
        confidence = "LOW"
        detected_segments = [seg for seg, data in segment_results.items() if data['detected']]
        message = f"Only {len(detected_segments)} segment(s) detected: {', '.join(detected_segments)}. Hantavirus requires ALL 3 segments."
    else:
        status = "NOT_DETECTED"
        confidence = "N/A"
        message = "No segments detected above threshold"

    return {
        'status': status,
        'virus': 'Hantavirus',
        'confidence': confidence,
        'segments': segment_results,
        'consensus_species': consensus_species,
        'species_distribution': dict(species_total),
        'method': '3-segment_concordance',
        'detection_criteria': 'L AND M AND S segments with ≥50 reads each',
        'message': message
    }


def detect_eeev(fastq: Path, db_path: Path, threads: int, output_dir: Path) -> Dict:
    """
    Detect EEEV (Eastern Equine Encephalitis Virus) from alphavirus database.

    Criteria:
    - ≥100 reads mapped to EEEV
    - ≥10× mean coverage
    - Phylogenetic assignment (EEEV vs other alphaviruses)
    """
    print("\n" + "="*70)
    print("EEEV (ALPHAVIRUS) DETECTION")
    print("="*70)

    config = VIRUS_CONFIG["eeev"]
    database = db_path / config["database"]

    if not database.exists():
        return {"status": "ERROR", "message": f"Database not found: {database}"}

    # Run alignment
    bam_file = output_dir / "alphavirus_alignments.bam"
    success = run_minimap2_alignment(fastq, database, threads, bam_file)

    if not success:
        return {"status": "ERROR", "message": "Alignment failed"}

    # Parse results
    try:
        import pysam

        bam = pysam.AlignmentFile(str(bam_file), "rb")

        alphavirus_stats = defaultdict(lambda: {
            'reads': 0,
            'bases_aligned': 0,
            'identity_scores': []
        })

        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            aligned_length = read.query_alignment_length
            if aligned_length < 100:
                continue

            try:
                nm = read.get_tag('NM')
                identity = 1.0 - (nm / aligned_length)
            except KeyError:
                matches = sum(count for op, count in read.cigartuples if op == 0)
                identity = matches / aligned_length if aligned_length > 0 else 0

            if identity < config["min_identity"]:
                continue

            ref_name = read.reference_name

            # Classify by virus
            virus_type = None
            if 'EEEV' in ref_name or 'Eastern' in ref_name:
                virus_type = 'EEEV'
            elif 'WEEV' in ref_name or 'Western' in ref_name:
                virus_type = 'WEEV'
            elif 'VEEV' in ref_name or 'Venezuelan' in ref_name:
                virus_type = 'VEEV'
            elif 'Getah' in ref_name:
                virus_type = 'Getah'
            else:
                virus_type = ref_name

            alphavirus_stats[virus_type]['reads'] += 1
            alphavirus_stats[virus_type]['bases_aligned'] += aligned_length
            alphavirus_stats[virus_type]['identity_scores'].append(identity)

        bam.close()

        # Evaluate results
        results = {}
        for virus, data in alphavirus_stats.items():
            if data['reads'] == 0:
                continue

            avg_identity = sum(data['identity_scores']) / len(data['identity_scores'])
            mean_depth = data['bases_aligned'] / config["genome_size"]

            detected = (data['reads'] >= config["min_reads"] and
                       mean_depth >= config["min_coverage"])

            results[virus] = {
                'detected': detected,
                'reads': data['reads'],
                'average_identity': round(avg_identity, 4),
                'mean_depth': round(mean_depth, 2)
            }

        # Determine which alphavirus is most likely
        if 'EEEV' in results and results['EEEV']['detected']:
            status = "DETECTED"
            primary_virus = "EEEV"
            confidence = "HIGH" if results['EEEV']['reads'] >= 200 else "MEDIUM"
        elif any(data['detected'] for data in results.values()):
            # Other alphavirus detected (not EEEV)
            primary_virus = max(results.items(), key=lambda x: x[1]['reads'])[0]
            status = "OTHER_ALPHAVIRUS"
            confidence = "MEDIUM"
        else:
            status = "NOT_DETECTED"
            primary_virus = None
            confidence = "N/A"

        return {
            'status': status,
            'virus': 'EEEV' if status == "DETECTED" else primary_virus,
            'confidence': confidence,
            'all_alphaviruses': results,
            'method': 'minimap2_alignment',
            'detection_criteria': f"≥{config['min_reads']} reads AND ≥{config['min_coverage']}× coverage",
            'note': 'Phylogenetic analysis recommended for lineage assignment (North American vs South American)'
        }

    except ImportError:
        return {"status": "ERROR", "message": "pysam not installed"}
    except Exception as e:
        return {"status": "ERROR", "message": str(e)}


def detect_spumavirus_metagenomic(fastq: Path, db_path: Path, threads: int, output_dir: Path) -> Dict:
    """
    Attempt spumavirus detection from metagenomic data.

    NOTE: This is expected to fail for most samples because:
    - No porcine spumavirus reference genome exists
    - Proviral DNA is ultra-rare (0.001% of genome)
    - Cross-genus identity too low (30-50%)

    Nested PCR approach is required for reliable detection.
    """
    print("\n" + "="*70)
    print("SPUMAVIRUS DETECTION (Metagenomic - Limited Sensitivity)")
    print("="*70)

    print("WARNING: Spumavirus detection from metagenomic sequencing has very low sensitivity.")
    print("         Nested PCR from PBMC genomic DNA is the recommended approach.")
    print("         See: scripts/phase4_pathogen/detect_spumavirus_nested_pcr.py")

    config = VIRUS_CONFIG["spumavirus"]

    # BLAST search against foamy virus pol genes
    database = db_path / config["database"]
    if not database.exists():
        return {
            "status": "NOT_APPLICABLE",
            "message": "Metagenomic detection not recommended for spumavirus",
            "recommendation": "Use nested PCR approach with PBMC genomic DNA"
        }

    # Convert FASTQ to FASTA for BLAST
    fasta_file = output_dir / "input.fasta"
    try:
        with open(fastq, 'r') as fq, open(fasta_file, 'w') as fa:
            line_num = 0
            for line in fq:
                if line_num % 4 == 0:  # Header line
                    fa.write('>' + line[1:])
                elif line_num % 4 == 1:  # Sequence line
                    fa.write(line)
                line_num += 1
    except Exception as e:
        return {"status": "ERROR", "message": f"FASTA conversion failed: {e}"}

    # Run BLAST
    blast_output = output_dir / "spumavirus_blast.txt"
    try:
        blast_cmd = [
            'blastn',
            '-query', str(fasta_file),
            '-db', str(database),
            '-out', str(blast_output),
            '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore',
            '-evalue', '1e-5',
            '-num_threads', str(threads),
            '-max_target_seqs', '10'
        ]

        subprocess.run(blast_cmd, check=True, capture_output=True)

    except subprocess.CalledProcessError as e:
        return {"status": "ERROR", "message": f"BLAST failed: {e.stderr.decode()}"}

    # Parse BLAST results
    hits = []
    try:
        with open(blast_output) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 10:
                    hits.append({
                        'query': fields[0],
                        'subject': fields[1],
                        'identity': float(fields[2]),
                        'length': int(fields[3]),
                        'evalue': float(fields[8]),
                        'bitscore': float(fields[9])
                    })
    except Exception as e:
        return {"status": "ERROR", "message": f"BLAST parsing failed: {e}"}

    # Filter hits
    significant_hits = [h for h in hits if h['identity'] >= 70 and h['length'] >= 100]

    # Classify hits
    sfv_hits = [h for h in significant_hits if 'SFV' in h['subject'] or 'Simian' in h['subject']]
    perv_hits = [h for h in significant_hits if 'PERV' in h['subject']]

    if sfv_hits:
        status = "POSSIBLE_DETECTION"
        confidence = "LOW"
        message = f"Found {len(sfv_hits)} hits to foamy virus pol gene (≥70% identity). Sanger confirmation required."
    elif significant_hits:
        status = "INCONCLUSIVE"
        confidence = "VERY_LOW"
        message = f"Found {len(significant_hits)} low-confidence hits. May be PERV or other retrovirus."
    else:
        status = "NOT_DETECTED"
        confidence = "N/A"
        message = "No significant BLAST hits to foamy virus references."

    return {
        'status': status,
        'virus': 'Porcine Spumavirus (metagenomic)',
        'confidence': confidence,
        'blast_hits': len(hits),
        'significant_hits': len(significant_hits),
        'sfv_hits': len(sfv_hits),
        'perv_hits': len(perv_hits),
        'method': 'blastn',
        'message': message,
        'recommendation': 'Use nested PCR from PBMC genomic DNA for definitive detection'
    }


def main():
    parser = argparse.ArgumentParser(
        description='PMDA 4-Virus Specific Detection with High Sensitivity',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Detect all 4 viruses
  python detect_pmda_4viruses.py --input reads.fastq --output results/ \\
      --db-path /mnt/efs/databases/pmda/2024.1/ --targets all

  # Detect specific viruses
  python detect_pmda_4viruses.py --input reads.fastq --output results/ \\
      --db-path /mnt/efs/databases/pmda/2024.1/ \\
      --targets polyomavirus hantavirus eeev

Note: Spumavirus requires separate nested PCR workflow (not metagenomic)
        """
    )

    parser.add_argument('--input', required=True, type=Path,
                       help='Input FASTQ file')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output directory')
    parser.add_argument('--db-path', required=True, type=Path,
                       help='PMDA database base path (e.g., /mnt/efs/databases/pmda/2024.1/)')
    parser.add_argument('--targets', nargs='+',
                       choices=['all', 'polyomavirus', 'hantavirus', 'eeev', 'spumavirus'],
                       default=['all'],
                       help='Target viruses to detect (default: all)')
    parser.add_argument('--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('--threads', type=int, default=16,
                       help='Number of threads (default: 16)')

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Determine targets
    if 'all' in args.targets:
        targets = ['polyomavirus', 'hantavirus', 'eeev', 'spumavirus']
    else:
        targets = args.targets

    print("="*70)
    print("PMDA 4-VIRUS SPECIFIC DETECTION")
    print("="*70)
    print(f"Run ID: {args.run_id}")
    print(f"Input: {args.input}")
    print(f"Database: {args.db_path}")
    print(f"Targets: {', '.join(targets)}")
    print(f"Threads: {args.threads}")
    print("="*70)

    # Validate input
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not args.db_path.exists():
        print(f"ERROR: Database path not found: {args.db_path}", file=sys.stderr)
        sys.exit(1)

    # Run detections
    results = {
        'run_id': args.run_id,
        'input_file': str(args.input),
        'database_path': str(args.db_path),
        'targets': targets,
        'detections': {}
    }

    if 'polyomavirus' in targets:
        results['detections']['polyomavirus'] = detect_polyomavirus(
            args.input, args.db_path, args.threads, args.output
        )

    if 'hantavirus' in targets:
        results['detections']['hantavirus'] = detect_hantavirus(
            args.input, args.db_path, args.threads, args.output
        )

    if 'eeev' in targets:
        results['detections']['eeev'] = detect_eeev(
            args.input, args.db_path, args.threads, args.output
        )

    if 'spumavirus' in targets:
        results['detections']['spumavirus'] = detect_spumavirus_metagenomic(
            args.input, args.db_path, args.threads, args.output
        )

    # Write results
    output_file = args.output / 'pmda_4virus_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print("\n" + "="*70)
    print("DETECTION SUMMARY")
    print("="*70)

    any_detected = False
    critical_detections = []

    for virus, result in results['detections'].items():
        status = result.get('status', 'UNKNOWN')
        print(f"\n{virus.upper()}: {status}")

        if status == "DETECTED":
            any_detected = True
            print(f"  Confidence: {result.get('confidence', 'N/A')}")
            if virus == 'hantavirus':
                print(f"  Species: {result.get('consensus_species', 'Unknown')}")
                for seg, data in result.get('segments', {}).items():
                    print(f"    {seg} segment: {data['reads']} reads")
            else:
                print(f"  Method: {result.get('method', 'N/A')}")
        elif status == "INCONCLUSIVE":
            print(f"  Message: {result.get('message', 'N/A')}")
        elif status == "ERROR":
            print(f"  ERROR: {result.get('message', 'Unknown error')}")

    print(f"\nResults written to: {output_file}")
    print("="*70)

    # Exit code
    if any_detected:
        sys.exit(0)  # Detection found
    else:
        sys.exit(0)  # No detection (not an error)


if __name__ == '__main__':
    main()
