#!/usr/bin/env python3
"""
MinION Metagenomics Pipeline - Phase 3: RNA Host Removal
Removes pig ribosomal RNA sequences to retain only viral/microbial RNA reads

This script complements the DNA host removal (remove_host.sh) by providing
RNA-specific filtering for:
1. Hantavirus (ssRNA-, no poly(A)) - requires rRNA depletion
2. EEEV (ssRNA+, with poly(A)) - poly(A) selected, but validate rRNA removal

PMDA 4-Virus Protocol Integration
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple
import tempfile
import gzip


class RNAHostRemover:
    """Removes host ribosomal RNA from metagenome samples."""

    def __init__(self, rrna_db: Path, threads: int = 8, min_identity: float = 85.0):
        """
        Initialize RNA host remover.

        Args:
            rrna_db: Path to pig rRNA reference database (minimap2 index)
            threads: Number of threads for alignment
            min_identity: Minimum identity threshold for rRNA match (%)
        """
        self.rrna_db = rrna_db
        self.threads = threads
        self.min_identity = min_identity

        if not self.rrna_db.exists():
            raise FileNotFoundError(f"rRNA database not found: {self.rrna_db}")

    def count_fastq_reads(self, fastq_path: Path) -> int:
        """Count reads in FASTQ file."""
        if str(fastq_path).endswith('.gz'):
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'

        with opener(fastq_path, mode) as f:
            return sum(1 for _ in f) // 4

    def align_to_rrna(self, input_fastq: Path, output_bam: Path) -> None:
        """
        Align reads to pig rRNA database using minimap2.

        Args:
            input_fastq: Input FASTQ file
            output_bam: Output BAM file path
        """
        print(f"Aligning reads to rRNA database: {self.rrna_db}")

        # Minimap2 for RNA alignment (splice-aware not needed for rRNA)
        cmd = [
            'minimap2',
            '-ax', 'map-ont',  # Oxford Nanopore preset
            '-t', str(self.threads),
            '--secondary=no',  # No secondary alignments
            str(self.rrna_db),
            str(input_fastq)
        ]

        # Pipe to samtools
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        cmd_samtools = ['samtools', 'view', '-bS', '-']
        with open(output_bam, 'wb') as f:
            p2 = subprocess.Popen(cmd_samtools, stdin=p1.stdout, stdout=f, stderr=subprocess.PIPE)
            p1.stdout.close()
            _, stderr2 = p2.communicate()

        _, stderr1 = p1.communicate()

        if p1.returncode != 0:
            raise RuntimeError(f"Minimap2 failed: {stderr1.decode()}")
        if p2.returncode != 0:
            raise RuntimeError(f"Samtools failed: {stderr2.decode()}")

    def extract_non_rrna_reads(self, aligned_bam: Path, output_fastq: Path) -> Tuple[int, int]:
        """
        Extract reads that did NOT align to rRNA (unmapped reads).

        Args:
            aligned_bam: BAM file with alignments to rRNA
            output_fastq: Output FASTQ file for non-rRNA reads

        Returns:
            Tuple of (total_reads, rrna_reads)
        """
        # Get statistics
        total_reads = int(subprocess.check_output(
            ['samtools', 'view', '-c', str(aligned_bam)]
        ).decode().strip())

        rrna_reads = int(subprocess.check_output(
            ['samtools', 'view', '-c', '-F', '4', str(aligned_bam)]
        ).decode().strip())

        print(f"Total reads: {total_reads:,}")
        print(f"rRNA reads (mapped): {rrna_reads:,} ({rrna_reads/total_reads*100:.2f}%)")
        print(f"Non-rRNA reads (unmapped): {total_reads - rrna_reads:,}")

        # Extract unmapped reads (non-rRNA) to FASTQ
        print(f"Extracting non-rRNA reads to {output_fastq}...")

        # samtools fastq with -f 4 flag (unmapped reads only)
        if str(output_fastq).endswith('.gz'):
            cmd = f"samtools fastq -f 4 {aligned_bam} | pigz -p {self.threads // 2} > {output_fastq}"
        else:
            cmd = f"samtools fastq -f 4 {aligned_bam} > {output_fastq}"

        subprocess.run(cmd, shell=True, check=True)

        return total_reads, rrna_reads

    def remove_rrna(self, input_fastq: Path, output_fastq: Path,
                    run_id: str, output_dir: Path) -> Dict:
        """
        Complete rRNA removal workflow.

        Args:
            input_fastq: Input FASTQ file
            output_fastq: Output FASTQ file (non-rRNA reads)
            run_id: Run identifier
            output_dir: Directory for intermediate files

        Returns:
            Dictionary with removal statistics
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Count input reads
        reads_before = self.count_fastq_reads(input_fastq)
        print(f"Input reads: {reads_before:,}")

        # Align to rRNA database
        aligned_bam = output_dir / "aligned_to_rrna.bam"
        self.align_to_rrna(input_fastq, aligned_bam)

        # Extract non-rRNA reads
        total_reads, rrna_reads = self.extract_non_rrna_reads(aligned_bam, output_fastq)
        non_rrna_reads = total_reads - rrna_reads

        # Calculate metrics
        rrna_depletion_rate = (rrna_reads / total_reads * 100) if total_reads > 0 else 0
        non_rrna_percentage = 100 - rrna_depletion_rate

        # PMDA 4-virus protocol requirement:
        # - For Hantavirus: rRNA should be <5% (95% depletion)
        # - For EEEV with poly(A) selection: rRNA should be <2% (98% depletion)
        metrics = {
            'run_id': run_id,
            'input_reads': reads_before,
            'total_reads_aligned': total_reads,
            'rrna_reads_removed': rrna_reads,
            'non_rrna_reads_retained': non_rrna_reads,
            'rrna_depletion_rate_percent': round(rrna_depletion_rate, 2),
            'non_rrna_percentage': round(non_rrna_percentage, 2),
            'hantavirus_protocol_compliant': non_rrna_percentage >= 95.0,
            'eeev_protocol_compliant': non_rrna_percentage >= 98.0,
            'min_identity_threshold': self.min_identity
        }

        return metrics

    def remove_host_genomic_rna(self, input_fastq: Path, output_fastq: Path,
                                host_genome: Path, output_dir: Path) -> Dict:
        """
        Remove pig genomic RNA sequences (transcriptome alignment).

        This is complementary to rRNA removal for comprehensive host depletion.

        Args:
            input_fastq: Input FASTQ (after rRNA removal)
            output_fastq: Output FASTQ (non-host RNA)
            host_genome: Path to pig genome reference
            output_dir: Directory for intermediate files

        Returns:
            Dictionary with host RNA removal statistics
        """
        print(f"\nRemoving host genomic RNA...")
        print(f"Host genome: {host_genome}")

        # Align to host genome
        aligned_bam = output_dir / "aligned_to_host_genome.bam"

        cmd = [
            'minimap2',
            '-ax', 'map-ont',
            '-t', str(self.threads),
            '--secondary=no',
            str(host_genome),
            str(input_fastq)
        ]

        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_samtools = ['samtools', 'view', '-bS', '-']
        with open(aligned_bam, 'wb') as f:
            p2 = subprocess.Popen(cmd_samtools, stdin=p1.stdout, stdout=f)
            p1.stdout.close()
            p2.communicate()
        p1.communicate()

        # Extract unmapped reads
        total_reads = int(subprocess.check_output(
            ['samtools', 'view', '-c', str(aligned_bam)]
        ).decode().strip())

        host_rna_reads = int(subprocess.check_output(
            ['samtools', 'view', '-c', '-F', '4', str(aligned_bam)]
        ).decode().strip())

        non_host_reads = total_reads - host_rna_reads

        print(f"Total reads: {total_reads:,}")
        print(f"Host genomic RNA: {host_rna_reads:,} ({host_rna_reads/total_reads*100:.2f}%)")
        print(f"Non-host RNA: {non_host_reads:,}")

        # Extract non-host reads
        if str(output_fastq).endswith('.gz'):
            cmd = f"samtools fastq -f 4 {aligned_bam} | pigz -p {self.threads // 2} > {output_fastq}"
        else:
            cmd = f"samtools fastq -f 4 {aligned_bam} > {output_fastq}"

        subprocess.run(cmd, shell=True, check=True)

        host_rna_depletion_rate = (host_rna_reads / total_reads * 100) if total_reads > 0 else 0

        return {
            'total_reads': total_reads,
            'host_rna_reads': host_rna_reads,
            'non_host_rna_reads': non_host_reads,
            'host_rna_depletion_rate_percent': round(host_rna_depletion_rate, 2),
            # PMDA requirement: >80% host RNA removal (RNA is less abundant than DNA)
            'pmda_compliant': host_rna_depletion_rate >= 80.0
        }


def main():
    parser = argparse.ArgumentParser(
        description='Remove pig host RNA (rRNA and genomic RNA) from metagenome samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # rRNA removal only (for Hantavirus protocol)
  %(prog)s -i input.fastq.gz -o output.fastq.gz -r RUN-001 --rrna-only

  # Complete host RNA removal (rRNA + genomic RNA)
  %(prog)s -i input.fastq.gz -o output.fastq.gz -r RUN-001 --complete

  # Custom rRNA database
  %(prog)s -i input.fastq.gz -o output.fastq.gz -r RUN-001 \\
      --rrna-db /custom/path/pig_rrna.mmi
        """
    )

    # Required arguments
    parser.add_argument('-i', '--input', required=True, type=Path,
                       help='Input FASTQ file (raw or post-QC)')
    parser.add_argument('-o', '--output', required=True, type=Path,
                       help='Output FASTQ file (non-host RNA reads)')
    parser.add_argument('-r', '--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('-d', '--output-dir', required=True, type=Path,
                       help='Output directory for intermediate files and metrics')

    # Database paths
    parser.add_argument('--rrna-db', type=Path,
                       default='/mnt/efs/databases/host/pig_rrna.mmi',
                       help='Pig rRNA database (minimap2 index)')
    parser.add_argument('--host-genome', type=Path,
                       default='/mnt/efs/host_genome/sus_scrofa_11.1.fa',
                       help='Pig genome reference for genomic RNA removal')

    # Workflow options
    parser.add_argument('--rrna-only', action='store_true',
                       help='Remove rRNA only (no genomic RNA removal)')
    parser.add_argument('--complete', action='store_true',
                       help='Complete host RNA removal (rRNA + genomic RNA)')

    # Parameters
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='Number of threads (default: 8)')
    parser.add_argument('--min-identity', type=float, default=85.0,
                       help='Minimum identity for rRNA match (default: 85)')

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not args.rrna_db.exists():
        print(f"ERROR: rRNA database not found: {args.rrna_db}", file=sys.stderr)
        print("Hint: Build with 'minimap2 -d pig_rrna.mmi pig_rrna.fasta'", file=sys.stderr)
        sys.exit(1)

    if args.complete and not args.host_genome.exists():
        print(f"ERROR: Host genome not found: {args.host_genome}", file=sys.stderr)
        sys.exit(1)

    # Default to complete removal if neither flag specified
    if not args.rrna_only and not args.complete:
        args.complete = True
        print("No workflow specified, defaulting to --complete")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("RNA Host Removal - PMDA 4-Virus Protocol")
    print("=" * 80)
    print(f"Run ID: {args.run_id}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Threads: {args.threads}")
    print(f"Workflow: {'rRNA only' if args.rrna_only else 'Complete (rRNA + genomic RNA)'}")
    print()

    # Initialize remover
    remover = RNAHostRemover(
        rrna_db=args.rrna_db,
        threads=args.threads,
        min_identity=args.min_identity
    )

    # Step 1: rRNA removal
    print("STEP 1: Ribosomal RNA Removal")
    print("-" * 80)

    if args.rrna_only:
        # Final output is after rRNA removal
        rrna_metrics = remover.remove_rrna(
            args.input, args.output, args.run_id, args.output_dir
        )
        final_metrics = rrna_metrics
    else:
        # Intermediate output after rRNA removal
        rrna_output = args.output_dir / "after_rrna_removal.fastq.gz"
        rrna_metrics = remover.remove_rrna(
            args.input, rrna_output, args.run_id, args.output_dir
        )

        # Step 2: Genomic host RNA removal
        print("\n" + "=" * 80)
        print("STEP 2: Genomic Host RNA Removal")
        print("-" * 80)

        genomic_metrics = remover.remove_host_genomic_rna(
            rrna_output, args.output, args.host_genome, args.output_dir
        )

        # Combine metrics
        final_metrics = {
            'run_id': args.run_id,
            'rrna_removal': rrna_metrics,
            'genomic_rna_removal': genomic_metrics,
            'workflow': 'complete'
        }

    # Save metrics
    metrics_file = args.output_dir / f"rna_host_removal_metrics.json"
    with open(metrics_file, 'w') as f:
        json.dump(final_metrics, f, indent=2)

    print("\n" + "=" * 80)
    print("RNA Host Removal Complete")
    print("=" * 80)

    # Print summary
    if args.rrna_only:
        print(f"rRNA depletion: {rrna_metrics['rrna_depletion_rate_percent']}%")
        print(f"Non-rRNA reads: {rrna_metrics['non_rrna_reads_retained']:,}")
        print(f"Hantavirus protocol compliant: {rrna_metrics['hantavirus_protocol_compliant']}")
        print(f"EEEV protocol compliant: {rrna_metrics['eeev_protocol_compliant']}")
    else:
        print(f"rRNA removal: {rrna_metrics['rrna_depletion_rate_percent']}%")
        print(f"Host RNA removal: {genomic_metrics['host_rna_depletion_rate_percent']}%")
        print(f"Final non-host reads: {genomic_metrics['non_host_rna_reads']:,}")

    print(f"\nMetrics saved to: {metrics_file}")
    print(f"Output saved to: {args.output}")

    # Upload to S3 if configured
    s3_bucket = subprocess.run(['printenv', 'S3_ANALYSIS_BUCKET'],
                              capture_output=True, text=True).stdout.strip()
    if s3_bucket:
        print(f"\nUploading results to S3...")
        subprocess.run([
            'aws', 's3', 'cp', str(args.output),
            f's3://{s3_bucket}/runs/{args.run_id}/host_removal_rna/'
        ], check=False)
        subprocess.run([
            'aws', 's3', 'cp', str(metrics_file),
            f's3://{s3_bucket}/runs/{args.run_id}/host_removal_rna/'
        ], check=False)
        print("Upload complete")


if __name__ == '__main__':
    main()
