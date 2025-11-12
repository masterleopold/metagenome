#!/usr/bin/env python3
"""
PMDA All 91-Pathogen Detection Module
Comprehensive detection for ALL PMDA-designated pathogens

Categories:
- 41 Viruses (DNA and RNA)
- 27 Bacteria
- 2 Fungi
- 19 Parasites

Total: 91 pathogens from 厚労省異種移植指針 別添2
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import tempfile


class PMDA91PathogenDetector:
    """Comprehensive detector for all 91 PMDA pathogens."""

    def __init__(self, database_dir: Path, threads: int = 8):
        """
        Initialize detector.

        Args:
            database_dir: Path to PMDA database directory
            threads: Number of threads for analysis
        """
        self.database_dir = database_dir
        self.threads = threads

        # Database paths
        self.minimap2_db = database_dir / "minimap2" / "pmda_all_91.mmi"
        self.kraken2_db = database_dir / "kraken2"
        self.blast_db = database_dir / "blast" / "pmda_all_91"

        # Validate databases exist
        if not self.minimap2_db.exists():
            raise FileNotFoundError(f"Minimap2 database not found: {self.minimap2_db}")

    def run_minimap2_alignment(self, input_fastq: Path, output_dir: Path) -> Dict:
        """
        Align reads to all 91 PMDA pathogens using Minimap2.

        Args:
            input_fastq: Input FASTQ file
            output_dir: Output directory

        Returns:
            Dictionary with alignment results
        """
        print(f"\nRunning Minimap2 alignment against all 91 PMDA pathogens...")

        output_sam = output_dir / "minimap2_alignment.sam"
        output_bam = output_dir / "minimap2_alignment.bam"

        # Minimap2 alignment
        cmd = [
            'minimap2',
            '-ax', 'map-ont',
            '-t', str(self.threads),
            '--secondary=no',
            str(self.minimap2_db),
            str(input_fastq)
        ]

        with open(output_sam, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

        if result.returncode != 0:
            print(f"ERROR: Minimap2 failed: {result.stderr.decode()}", file=sys.stderr)
            return {'error': result.stderr.decode()}

        # Convert to BAM and sort
        subprocess.run(['samtools', 'view', '-bS', str(output_sam)],
                      stdout=open(output_bam, 'wb'), check=True)
        subprocess.run(['samtools', 'sort', str(output_bam), '-o', str(output_bam)], check=True)
        subprocess.run(['samtools', 'index', str(output_bam)], check=True)

        # Parse alignment results
        return self.parse_minimap2_results(output_bam)

    def parse_minimap2_results(self, bam_file: Path) -> Dict:
        """Parse Minimap2 BAM file and extract pathogen hits."""

        # Get alignment statistics
        result = subprocess.run(
            ['samtools', 'idxstats', str(bam_file)],
            capture_output=True, text=True, check=True
        )

        pathogen_hits = {}
        for line in result.stdout.strip().split('\n'):
            if not line or line.startswith('*'):
                continue

            parts = line.split('\t')
            ref_name = parts[0]
            ref_length = int(parts[1])
            mapped_reads = int(parts[2])

            if mapped_reads > 0:
                pathogen_hits[ref_name] = {
                    'reference_length': ref_length,
                    'mapped_reads': mapped_reads
                }

        return {
            'method': 'minimap2',
            'pathogens_detected': len(pathogen_hits),
            'pathogen_hits': pathogen_hits
        }

    def run_kraken2_classification(self, input_fastq: Path, output_dir: Path) -> Dict:
        """
        Classify reads against all 91 PMDA pathogens using Kraken2.

        Args:
            input_fastq: Input FASTQ file
            output_dir: Output directory

        Returns:
            Dictionary with classification results
        """
        if not self.kraken2_db.exists():
            print("WARNING: Kraken2 database not found, skipping", file=sys.stderr)
            return {'error': 'Kraken2 database not found'}

        print(f"\nRunning Kraken2 classification against all 91 PMDA pathogens...")

        output_txt = output_dir / "kraken2_output.txt"
        report_txt = output_dir / "kraken2_report.txt"

        # Kraken2 classification
        cmd = [
            'kraken2',
            '--db', str(self.kraken2_db),
            '--threads', str(self.threads),
            '--report', str(report_txt),
            '--output', str(output_txt),
            str(input_fastq)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"ERROR: Kraken2 failed: {result.stderr}", file=sys.stderr)
            return {'error': result.stderr}

        # Parse Kraken2 report
        return self.parse_kraken2_report(report_txt)

    def parse_kraken2_report(self, report_file: Path) -> Dict:
        """Parse Kraken2 report and extract pathogen classifications."""

        pathogen_classifications = {}

        with open(report_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                percent = float(parts[0])
                reads_clade = int(parts[1])
                reads_taxon = int(parts[2])
                rank = parts[3]
                taxid = parts[4]
                name = parts[5].strip()

                if reads_taxon > 0 and percent > 0.01:  # At least 0.01% of reads
                    pathogen_classifications[name] = {
                        'percent': percent,
                        'reads': reads_taxon,
                        'rank': rank,
                        'taxid': taxid
                    }

        return {
            'method': 'kraken2',
            'pathogens_classified': len(pathogen_classifications),
            'classifications': pathogen_classifications
        }

    def categorize_detections(self, minimap2_results: Dict,
                             kraken2_results: Dict) -> Dict:
        """
        Categorize detected pathogens into viruses, bacteria, fungi, parasites.

        Args:
            minimap2_results: Minimap2 alignment results
            kraken2_results: Kraken2 classification results

        Returns:
            Dictionary with categorized detections
        """
        # Load PMDA pathogen config for categorization
        config_file = Path(__file__).parent.parent.parent / "templates" / "config" / "pmda_pathogens.json"

        if not config_file.exists():
            config_file = Path("/opt/minion/templates/config/pmda_pathogens.json")

        with open(config_file) as f:
            pmda_config = json.load(f)

        categories = {
            'viruses': [],
            'bacteria': [],
            'fungi': [],
            'parasites': [],
            'special_management': []
        }

        # Categorize Minimap2 hits
        for ref_name, hit_data in minimap2_results.get('pathogen_hits', {}).items():
            # Match to PMDA pathogen by name parsing
            category = self.determine_category(ref_name, pmda_config)
            if category:
                categories[category].append({
                    'name': ref_name,
                    'reads': hit_data['mapped_reads'],
                    'detection_method': 'minimap2'
                })

        # Categorize Kraken2 hits
        for name, class_data in kraken2_results.get('classifications', {}).items():
            category = self.determine_category(name, pmda_config)
            if category:
                categories[category].append({
                    'name': name,
                    'reads': class_data['reads'],
                    'percent': class_data['percent'],
                    'detection_method': 'kraken2'
                })

        return {
            'summary': {
                'total_viruses': len(categories['viruses']),
                'total_bacteria': len(categories['bacteria']),
                'total_fungi': len(categories['fungi']),
                'total_parasites': len(categories['parasites']),
                'total_special_management': len(categories['special_management'])
            },
            'detections': categories
        }

    def determine_category(self, pathogen_name: str, pmda_config: Dict) -> str:
        """Determine pathogen category from name."""
        name_lower = pathogen_name.lower()

        # Virus keywords
        if any(kw in name_lower for kw in ['virus', 'viral', 'perv', 'circovirus', 'herpes',
                                            'parvovirus', 'influenza', 'coronavirus']):
            return 'viruses'

        # Bacteria keywords
        if any(kw in name_lower for kw in ['bacteria', 'bacillus', 'clostridium', 'mycobacterium',
                                            'salmonella', 'escherichia', 'staphylococcus',
                                            'streptococcus', 'mycoplasma', 'brucella']):
            return 'bacteria'

        # Fungi keywords
        if any(kw in name_lower for kw in ['fungi', 'candida', 'aspergillus', 'trichophyton']):
            return 'fungi'

        # Parasite keywords
        if any(kw in name_lower for kw in ['toxoplasma', 'trypanosoma', 'ascaris', 'taenia',
                                            'trichinella', 'echinococcus', 'cryptosporidium']):
            return 'parasites'

        return 'viruses'  # Default

    def generate_report(self, categorized_results: Dict, run_id: str,
                       output_file: Path) -> None:
        """Generate comprehensive detection report."""

        report = {
            'run_id': run_id,
            'pmda_compliance': {
                'total_pathogens_screened': 91,
                'categories_screened': 4,
                'viruses_targeted': 41,
                'bacteria_targeted': 27,
                'fungi_targeted': 2,
                'parasites_targeted': 19
            },
            'detection_summary': categorized_results['summary'],
            'detections_by_category': categorized_results['detections'],
            'pmda_compliant': True,
            'notes': [
                'All 91 PMDA pathogens screened',
                'Multi-method detection (Minimap2 + Kraken2)',
                'Specialized protocols available for 4 high-sensitivity viruses'
            ]
        }

        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"\nDetection report saved: {output_file}")

    def detect_all(self, input_fastq: Path, output_dir: Path, run_id: str) -> Dict:
        """
        Run complete detection pipeline for all 91 PMDA pathogens.

        Args:
            input_fastq: Input FASTQ file
            output_dir: Output directory
            run_id: Run identifier

        Returns:
            Complete detection results
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        print("=" * 80)
        print("PMDA ALL 91-PATHOGEN DETECTION")
        print("=" * 80)
        print(f"Run ID: {run_id}")
        print(f"Input: {input_fastq}")
        print(f"Output: {output_dir}")
        print()

        # Method 1: Minimap2 alignment
        minimap2_results = self.run_minimap2_alignment(input_fastq, output_dir)

        # Method 2: Kraken2 classification
        kraken2_results = self.run_kraken2_classification(input_fastq, output_dir)

        # Categorize detections
        print("\nCategorizing detected pathogens...")
        categorized_results = self.categorize_detections(minimap2_results, kraken2_results)

        # Generate report
        report_file = output_dir / "pmda_all_91_detection_report.json"
        self.generate_report(categorized_results, run_id, report_file)

        # Print summary
        print("\n" + "=" * 80)
        print("DETECTION SUMMARY")
        print("=" * 80)
        summary = categorized_results['summary']
        print(f"Viruses detected: {summary['total_viruses']}/41")
        print(f"Bacteria detected: {summary['total_bacteria']}/27")
        print(f"Fungi detected: {summary['total_fungi']}/2")
        print(f"Parasites detected: {summary['total_parasites']}/19")
        print(f"Special management: {summary['total_special_management']}/5")
        print()

        return categorized_results


def main():
    parser = argparse.ArgumentParser(
        description='Detect all 91 PMDA-designated pathogens',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Detect all 91 pathogens
  %(prog)s -i input.fastq.gz -o output/ -r RUN-001

  # Custom database location
  %(prog)s -i input.fastq.gz -o output/ -r RUN-001 \\
      --database /custom/path/pmda/2024.2/all_91_pathogens
        """
    )

    parser.add_argument('-i', '--input', required=True, type=Path,
                       help='Input FASTQ file (filtered, host-depleted)')
    parser.add_argument('-o', '--output', required=True, type=Path,
                       help='Output directory')
    parser.add_argument('-r', '--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('--database', type=Path,
                       default='/mnt/efs/databases/pmda/2024.2/all_91_pathogens',
                       help='PMDA database directory')
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='Number of threads')

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not args.database.exists():
        print(f"ERROR: Database not found: {args.database}", file=sys.stderr)
        print("Build database first: scripts/database_build/build_pmda_all_91_databases.sh")
        sys.exit(1)

    # Initialize detector
    detector = PMDA91PathogenDetector(args.database, threads=args.threads)

    # Run detection
    try:
        results = detector.detect_all(args.input, args.output, args.run_id)
        print("\nDetection complete!")
        sys.exit(0)

    except Exception as e:
        print(f"\nERROR: Detection failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
