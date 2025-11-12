#!/usr/bin/env python3
"""
MinION Metagenomics Pipeline - Phase 3: Host Removal Orchestrator
Routes samples to appropriate host removal workflow based on sample type

Supports:
- DNA samples: Minimap2 alignment to Sus scrofa genome (existing remove_host.sh)
- RNA samples: rRNA depletion + genomic RNA removal (new remove_host_rna.py)
- Dual samples: Both DNA and RNA protocols

PMDA 4-Virus Protocol Integration:
- Polyomavirus (dsDNA) → DNA host removal
- Hantavirus (ssRNA-) → RNA host removal (rRNA depletion)
- EEEV (ssRNA+) → RNA host removal (poly(A) selected, validate rRNA)
- Spumavirus (proviral DNA) → DNA host removal
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
from enum import Enum


class SampleType(Enum):
    """Sample type enumeration."""
    DNA = "dna"
    RNA = "rna"
    DUAL = "dual"  # Both DNA and RNA extracted


class HostRemovalOrchestrator:
    """Orchestrates host removal workflows for DNA/RNA samples."""

    def __init__(self, script_dir: Path, threads: int = 8):
        """
        Initialize orchestrator.

        Args:
            script_dir: Directory containing host removal scripts
            threads: Number of threads for parallel processing
        """
        self.script_dir = script_dir
        self.threads = threads

        # Script paths
        self.dna_removal_script = script_dir / "remove_host.sh"
        self.rna_removal_script = script_dir / "remove_host_rna.py"

        # Validate scripts exist
        if not self.dna_removal_script.exists():
            raise FileNotFoundError(f"DNA removal script not found: {self.dna_removal_script}")
        if not self.rna_removal_script.exists():
            raise FileNotFoundError(f"RNA removal script not found: {self.rna_removal_script}")

    def remove_dna_host(self, input_fastq: Path, output_dir: Path, run_id: str,
                       host_genome: Path) -> Dict:
        """
        Execute DNA host removal workflow.

        Args:
            input_fastq: Input FASTQ file
            output_dir: Output directory
            run_id: Run identifier
            host_genome: Path to Sus scrofa genome

        Returns:
            Dictionary with removal metrics
        """
        print("\n" + "=" * 80)
        print("DNA HOST REMOVAL WORKFLOW")
        print("=" * 80)

        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            'bash', str(self.dna_removal_script),
            '-i', str(input_fastq),
            '-o', str(output_dir),
            '-r', run_id,
            '-g', str(host_genome),
            '-t', str(self.threads)
        ]

        print(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"ERROR: DNA host removal failed", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(1)

        print(result.stdout)

        # Load metrics
        metrics_file = output_dir / "host_removal_stats.json"
        if metrics_file.exists():
            with open(metrics_file) as f:
                return json.load(f)
        else:
            return {'status': 'completed', 'metrics_file': 'not_found'}

    def remove_rna_host(self, input_fastq: Path, output_dir: Path, run_id: str,
                       rrna_db: Path, host_genome: Path, workflow: str = 'complete') -> Dict:
        """
        Execute RNA host removal workflow.

        Args:
            input_fastq: Input FASTQ file
            output_dir: Output directory
            run_id: Run identifier
            rrna_db: Path to pig rRNA database
            host_genome: Path to Sus scrofa genome
            workflow: 'rrna_only' or 'complete'

        Returns:
            Dictionary with removal metrics
        """
        print("\n" + "=" * 80)
        print("RNA HOST REMOVAL WORKFLOW")
        print("=" * 80)

        output_dir.mkdir(parents=True, exist_ok=True)
        output_fastq = output_dir / "non_host_rna_reads.fastq.gz"

        cmd = [
            'python3', str(self.rna_removal_script),
            '-i', str(input_fastq),
            '-o', str(output_fastq),
            '-r', run_id,
            '-d', str(output_dir),
            '--rrna-db', str(rrna_db),
            '--host-genome', str(host_genome),
            '-t', str(self.threads)
        ]

        if workflow == 'rrna_only':
            cmd.append('--rrna-only')
        else:
            cmd.append('--complete')

        print(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"ERROR: RNA host removal failed", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(1)

        print(result.stdout)

        # Load metrics
        metrics_file = output_dir / "rna_host_removal_metrics.json"
        if metrics_file.exists():
            with open(metrics_file) as f:
                return json.load(f)
        else:
            return {'status': 'completed', 'metrics_file': 'not_found'}

    def determine_workflow(self, target_viruses: Optional[List[str]] = None) -> SampleType:
        """
        Determine appropriate host removal workflow based on target viruses.

        Args:
            target_viruses: List of PMDA virus codes to detect

        Returns:
            SampleType (DNA, RNA, or DUAL)
        """
        if not target_viruses:
            # Default: assume DNA workflow for general pathogen screening
            return SampleType.DNA

        # PMDA 4-virus classification
        dna_viruses = {'POLYOMA', 'SPUMV'}  # Polyomavirus (dsDNA), Spumavirus (proviral DNA)
        rna_viruses = {'HANTV', 'EEEV'}  # Hantavirus (ssRNA-), EEEV (ssRNA+)

        has_dna = any(v in dna_viruses for v in target_viruses)
        has_rna = any(v in rna_viruses for v in target_viruses)

        if has_dna and has_rna:
            return SampleType.DUAL
        elif has_rna:
            return SampleType.RNA
        else:
            return SampleType.DNA

    def execute_workflow(self, input_fastq: Path, output_base_dir: Path, run_id: str,
                        sample_type: SampleType, host_genome: Path,
                        rrna_db: Optional[Path] = None) -> Dict:
        """
        Execute complete host removal workflow.

        Args:
            input_fastq: Input FASTQ file
            output_base_dir: Base output directory
            run_id: Run identifier
            sample_type: DNA, RNA, or DUAL
            host_genome: Path to Sus scrofa genome
            rrna_db: Path to pig rRNA database (required for RNA/DUAL)

        Returns:
            Dictionary with all metrics
        """
        results = {
            'run_id': run_id,
            'sample_type': sample_type.value,
            'workflows_executed': []
        }

        if sample_type == SampleType.DNA:
            # DNA-only workflow
            dna_output = output_base_dir / "dna"
            dna_metrics = self.remove_dna_host(input_fastq, dna_output, run_id, host_genome)
            results['dna_host_removal'] = dna_metrics
            results['workflows_executed'].append('dna')

        elif sample_type == SampleType.RNA:
            # RNA-only workflow
            if not rrna_db:
                raise ValueError("RNA workflow requires --rrna-db parameter")

            rna_output = output_base_dir / "rna"
            rna_metrics = self.remove_rna_host(
                input_fastq, rna_output, run_id, rrna_db, host_genome, workflow='complete'
            )
            results['rna_host_removal'] = rna_metrics
            results['workflows_executed'].append('rna')

        elif sample_type == SampleType.DUAL:
            # Both DNA and RNA workflows
            if not rrna_db:
                raise ValueError("Dual workflow requires --rrna-db parameter")

            # DNA workflow
            dna_output = output_base_dir / "dna"
            dna_metrics = self.remove_dna_host(input_fastq, dna_output, run_id, host_genome)
            results['dna_host_removal'] = dna_metrics

            # RNA workflow
            rna_output = output_base_dir / "rna"
            rna_metrics = self.remove_rna_host(
                input_fastq, rna_output, run_id, rrna_db, host_genome, workflow='complete'
            )
            results['rna_host_removal'] = rna_metrics

            results['workflows_executed'] = ['dna', 'rna']

        return results


def main():
    parser = argparse.ArgumentParser(
        description='Orchestrate host removal workflows for DNA/RNA samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect workflow based on target viruses
  %(prog)s -i input.fastq.gz -o output/ -r RUN-001 \\
      --target-viruses POLYOMA HANTV

  # Force DNA workflow
  %(prog)s -i input.fastq.gz -o output/ -r RUN-001 --sample-type dna

  # Force RNA workflow
  %(prog)s -i input.fastq.gz -o output/ -r RUN-001 --sample-type rna \\
      --rrna-db /mnt/efs/databases/host/pig_rrna.mmi

  # Dual workflow (both DNA and RNA)
  %(prog)s -i input.fastq.gz -o output/ -r RUN-001 --sample-type dual \\
      --rrna-db /mnt/efs/databases/host/pig_rrna.mmi
        """
    )

    # Required arguments
    parser.add_argument('-i', '--input', required=True, type=Path,
                       help='Input FASTQ file')
    parser.add_argument('-o', '--output', required=True, type=Path,
                       help='Output base directory')
    parser.add_argument('-r', '--run-id', required=True,
                       help='Run identifier')

    # Workflow selection
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--sample-type', type=str, choices=['dna', 'rna', 'dual'],
                      help='Force specific workflow type')
    group.add_argument('--target-viruses', nargs='+',
                      help='PMDA virus codes to detect (auto-determines workflow)')

    # Database paths
    parser.add_argument('--host-genome', type=Path,
                       default='/mnt/efs/host_genome/sus_scrofa_11.1.fa',
                       help='Sus scrofa genome reference')
    parser.add_argument('--rrna-db', type=Path,
                       default='/mnt/efs/databases/host/pig_rrna.mmi',
                       help='Pig rRNA database (required for RNA/DUAL workflows)')

    # Parameters
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='Number of threads')

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not args.host_genome.exists():
        print(f"ERROR: Host genome not found: {args.host_genome}", file=sys.stderr)
        sys.exit(1)

    # Determine workflow
    script_dir = Path(__file__).parent
    orchestrator = HostRemovalOrchestrator(script_dir, threads=args.threads)

    if args.sample_type:
        sample_type = SampleType(args.sample_type)
    elif args.target_viruses:
        sample_type = orchestrator.determine_workflow(args.target_viruses)
        print(f"Auto-detected workflow: {sample_type.value.upper()}")
    else:
        # Default to DNA
        sample_type = SampleType.DNA
        print("No workflow specified, defaulting to DNA")

    # Validate RNA/DUAL workflow requirements
    if sample_type in [SampleType.RNA, SampleType.DUAL]:
        if not args.rrna_db.exists():
            print(f"ERROR: rRNA database required for {sample_type.value} workflow", file=sys.stderr)
            print(f"Specified: {args.rrna_db}", file=sys.stderr)
            sys.exit(1)

    # Execute workflow
    print("=" * 80)
    print("HOST REMOVAL ORCHESTRATOR - PMDA 4-Virus Protocol")
    print("=" * 80)
    print(f"Run ID: {args.run_id}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Workflow: {sample_type.value.upper()}")
    print(f"Threads: {args.threads}")

    if args.target_viruses:
        print(f"Target viruses: {', '.join(args.target_viruses)}")

    try:
        results = orchestrator.execute_workflow(
            input_fastq=args.input,
            output_base_dir=args.output,
            run_id=args.run_id,
            sample_type=sample_type,
            host_genome=args.host_genome,
            rrna_db=args.rrna_db if sample_type != SampleType.DNA else None
        )

        # Save combined results
        results_file = args.output / "host_removal_orchestrator_results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        print("\n" + "=" * 80)
        print("HOST REMOVAL COMPLETE")
        print("=" * 80)
        print(f"Workflows executed: {', '.join(results['workflows_executed'])}")
        print(f"Results saved to: {results_file}")

        # Print summary statistics
        if 'dna_host_removal' in results:
            dna = results['dna_host_removal']
            if 'host_removal_rate' in dna:
                print(f"\nDNA Host Removal: {dna['host_removal_rate']}%")
                print(f"Non-host DNA reads: {dna.get('non_host_reads', 'N/A'):,}")

        if 'rna_host_removal' in results:
            rna = results['rna_host_removal']
            if 'rrna_removal' in rna:
                rrna = rna['rrna_removal']
                print(f"\nrRNA Removal: {rrna['rrna_depletion_rate_percent']}%")
            if 'genomic_rna_removal' in rna:
                grna = rna['genomic_rna_removal']
                print(f"Genomic RNA Removal: {grna['host_rna_depletion_rate_percent']}%")

        # Upload to S3 if configured
        s3_bucket = subprocess.run(['printenv', 'S3_ANALYSIS_BUCKET'],
                                  capture_output=True, text=True).stdout.strip()
        if s3_bucket:
            print(f"\nUploading results to S3...")
            subprocess.run([
                'aws', 's3', 'sync', str(args.output),
                f's3://{s3_bucket}/runs/{args.run_id}/host_removal/'
            ], check=False)

    except Exception as e:
        print(f"\nERROR: Workflow failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
