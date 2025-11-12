#!/usr/bin/env python3
"""
Spumavirus Nested PCR Detection and Confirmation Workflow
PMDA 4-Virus Protocol - Spumavirus Detection Module

CRITICAL: No porcine spumavirus reference genome exists in NCBI.
Detection requires cross-genus amplification with degenerate primers.

Workflow:
1. Nested PCR from PBMC genomic DNA (NOT plasma)
2. Agarose gel electrophoresis (expected: ~400 bp)
3. Sanger sequencing confirmation (both strands)
4. BLAST analysis against NCBI nt database
5. Phylogenetic analysis (Spumaretrovirinae vs Gammaretrovirus/PERV)
6. Final determination: Spumavirus vs PERV discrimination

This script handles the bioinformatics components (steps 4-6).
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import tempfile
from Bio import SeqIO, Phylo, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


# Degenerate primer sequences for nested PCR
NESTED_PCR_PRIMERS = {
    'outer': {
        'forward': 'GGNCARATHGGNATGTTYGG',  # FV-pol-F1 (96-fold degeneracy)
        'reverse': 'CCRTCNCCRAANCCRTC',     # FV-pol-R1 (64-fold degeneracy)
        'expected_product': 800  # bp
    },
    'inner': {
        'forward': 'ATHGGNCARGGNTTYACNAC',  # FV-pol-F2 (96-fold degeneracy)
        'reverse': 'GTRTCNGTYTTRTCNCC',     # FV-pol-R2 (64-fold degeneracy)
        'expected_product': 400  # bp
    }
}

# Reference sequences for phylogenetic analysis
REFERENCE_SEQUENCES = {
    'SFV': {
        'name': 'Simian foamy virus pol gene',
        'accession': 'NC_001364',
        'region': '2000-5000',
        'genus': 'Spumavirus'
    },
    'FFV': {
        'name': 'Feline foamy virus pol gene',
        'accession': 'NC_001871',
        'region': '2100-5100',
        'genus': 'Spumavirus'
    },
    'BFV': {
        'name': 'Bovine foamy virus pol gene',
        'accession': 'NC_001831',
        'region': '2050-5050',
        'genus': 'Spumavirus'
    },
    'PERV': {
        'name': 'Porcine endogenous retrovirus pol gene',
        'accession': 'AF038600',
        'region': '4500-6500',
        'genus': 'Gammaretrovirus'
    }
}


class SpumavirusDetector:
    """Spumavirus detection and discrimination from PERV."""

    def __init__(self, work_dir: Path, blast_db: Optional[Path] = None):
        """
        Initialize detector.

        Args:
            work_dir: Working directory for intermediate files
            blast_db: Local BLAST database (optional, uses NCBI if None)
        """
        self.work_dir = work_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.blast_db = blast_db

    def validate_sequence_length(self, sequence_file: Path, min_length: int = 300) -> Tuple[bool, int]:
        """
        Validate Sanger sequence length.

        Args:
            sequence_file: FASTA file with Sanger sequence
            min_length: Minimum expected length (default: 300 bp)

        Returns:
            Tuple of (is_valid, sequence_length)
        """
        for record in SeqIO.parse(sequence_file, "fasta"):
            seq_len = len(record.seq)
            is_valid = seq_len >= min_length
            return is_valid, seq_len

        return False, 0

    def blast_sequence(self, sequence_file: Path, output_xml: Path) -> Dict:
        """
        BLAST sequence against NCBI nt database.

        Args:
            sequence_file: FASTA file with query sequence
            output_xml: Output XML file path

        Returns:
            Dictionary with top BLAST hits
        """
        print("Running BLAST search against NCBI nt database...")
        print("(This may take 5-10 minutes)")

        # Read query sequence
        query_record = SeqIO.read(sequence_file, "fasta")
        query_seq = str(query_record.seq)

        if self.blast_db:
            # Local BLAST
            cmd = [
                'blastn',
                '-query', str(sequence_file),
                '-db', str(self.blast_db),
                '-out', str(output_xml),
                '-outfmt', '5',  # XML format
                '-max_target_seqs', '20',
                '-evalue', '1e-5'
            ]
            subprocess.run(cmd, check=True)
        else:
            # NCBI BLAST (remote)
            result_handle = NCBIWWW.qblast("blastn", "nt", query_seq, hitlist_size=20)
            with open(output_xml, 'w') as f:
                f.write(result_handle.read())
            result_handle.close()

        print(f"BLAST results saved to: {output_xml}")

        # Parse BLAST XML
        return self.parse_blast_results(output_xml)

    def parse_blast_results(self, blast_xml: Path) -> Dict:
        """
        Parse BLAST XML results.

        Args:
            blast_xml: BLAST XML output file

        Returns:
            Dictionary with top hits
        """
        with open(blast_xml) as f:
            blast_records = NCBIXML.parse(f)
            blast_record = next(blast_records)

        hits = []
        for alignment in blast_record.alignments[:10]:  # Top 10 hits
            for hsp in alignment.hsps[:1]:  # Best HSP per alignment
                hit = {
                    'title': alignment.title,
                    'accession': alignment.accession,
                    'length': alignment.length,
                    'e_value': hsp.expect,
                    'identity': hsp.identities / hsp.align_length * 100,
                    'align_length': hsp.align_length,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end
                }
                hits.append(hit)

        return {
            'query_length': blast_record.query_length,
            'hits': hits
        }

    def classify_virus(self, blast_results: Dict) -> Dict:
        """
        Classify as Spumavirus vs PERV based on BLAST results.

        Classification logic:
        - Spumavirus: Top hit to Spumaretrovirinae with ≥70% identity
        - PERV: Top hit to Gammaretrovirus/PERV with ≥80% identity
        - Ambiguous: Neither criteria met

        Args:
            blast_results: Parsed BLAST results

        Returns:
            Classification results
        """
        if not blast_results['hits']:
            return {
                'classification': 'unknown',
                'confidence': 'none',
                'reason': 'No BLAST hits found'
            }

        top_hit = blast_results['hits'][0]
        title_lower = top_hit['title'].lower()
        identity = top_hit['identity']

        # Check for Spumavirus keywords
        spumavirus_keywords = ['spumavirus', 'foamy virus', 'spumaretrovir']
        is_spumavirus = any(kw in title_lower for kw in spumavirus_keywords)

        # Check for PERV keywords
        perv_keywords = ['perv', 'porcine endogenous retrovirus', 'gammaretrovirus']
        is_perv = any(kw in title_lower for kw in perv_keywords)

        if is_spumavirus and identity >= 70:
            return {
                'classification': 'spumavirus',
                'confidence': 'high' if identity >= 80 else 'medium',
                'top_hit': top_hit['title'],
                'identity': round(identity, 2),
                'reason': f'Top BLAST hit to Spumavirus with {identity:.1f}% identity'
            }
        elif is_perv and identity >= 80:
            return {
                'classification': 'PERV',
                'confidence': 'high',
                'top_hit': top_hit['title'],
                'identity': round(identity, 2),
                'reason': f'Top BLAST hit to PERV with {identity:.1f}% identity',
                'warning': 'PERV is endogenous to all pig genomes - NOT a contamination risk'
            }
        else:
            return {
                'classification': 'ambiguous',
                'confidence': 'low',
                'top_hit': top_hit['title'],
                'identity': round(identity, 2),
                'reason': f'Top hit identity ({identity:.1f}%) below thresholds'
            }

    def download_reference_sequences(self, output_dir: Path) -> Path:
        """
        Download reference sequences for phylogenetic analysis.

        Args:
            output_dir: Directory to save reference sequences

        Returns:
            Path to combined reference FASTA file
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        ref_sequences = []

        for ref_code, ref_info in REFERENCE_SEQUENCES.items():
            accession = ref_info['accession']
            region = ref_info['region']

            # Download from NCBI
            cmd = f"esearch -db nucleotide -query '{accession}' | efetch -format fasta"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode == 0 and result.stdout:
                ref_file = output_dir / f"{ref_code}_pol.fasta"
                with open(ref_file, 'w') as f:
                    f.write(result.stdout)
                ref_sequences.append(ref_file)
                print(f"Downloaded {ref_code}: {ref_info['name']}")
            else:
                print(f"WARNING: Failed to download {ref_code}", file=sys.stderr)

        # Combine references
        combined = output_dir / "references_combined.fasta"
        with open(combined, 'w') as outf:
            for ref_file in ref_sequences:
                with open(ref_file) as inf:
                    outf.write(inf.read())

        return combined

    def build_phylogenetic_tree(self, query_fasta: Path, reference_fasta: Path,
                               output_tree: Path) -> Dict:
        """
        Build phylogenetic tree with query and references.

        Args:
            query_fasta: Query sequence FASTA
            reference_fasta: Reference sequences FASTA
            output_tree: Output Newick tree file

        Returns:
            Tree analysis results
        """
        print("\nBuilding phylogenetic tree...")

        # Combine query and references
        combined_fasta = self.work_dir / "combined_for_alignment.fasta"
        with open(combined_fasta, 'w') as outf:
            with open(query_fasta) as f:
                outf.write(f.read())
            with open(reference_fasta) as f:
                outf.write(f.read())

        # Multiple sequence alignment with MUSCLE
        aligned_fasta = self.work_dir / "aligned.fasta"
        muscle_cline = MuscleCommandline(input=str(combined_fasta), out=str(aligned_fasta))

        try:
            stdout, stderr = muscle_cline()
            print("Multiple sequence alignment completed")
        except Exception as e:
            print(f"ERROR: MUSCLE alignment failed: {e}", file=sys.stderr)
            return {'error': str(e)}

        # Build distance tree
        try:
            alignment = AlignIO.read(aligned_fasta, "fasta")
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(alignment)

            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)  # Neighbor-joining tree

            # Save tree
            Phylo.write(tree, output_tree, "newick")
            print(f"Phylogenetic tree saved to: {output_tree}")

            # Analyze tree topology
            query_clade = [c for c in tree.get_terminals() if 'query' in c.name.lower()][0]
            nearest_neighbor = self.find_nearest_neighbor(tree, query_clade)

            return {
                'tree_file': str(output_tree),
                'nearest_neighbor': nearest_neighbor.name if nearest_neighbor else 'unknown',
                'clustering': 'Spumavirus' if 'SFV' in nearest_neighbor.name or 'FFV' in nearest_neighbor.name
                                            or 'BFV' in nearest_neighbor.name else 'PERV'
            }

        except Exception as e:
            print(f"ERROR: Tree construction failed: {e}", file=sys.stderr)
            return {'error': str(e)}

    def find_nearest_neighbor(self, tree, query_clade):
        """Find nearest neighbor to query in phylogenetic tree."""
        min_distance = float('inf')
        nearest = None

        for terminal in tree.get_terminals():
            if terminal != query_clade:
                distance = tree.distance(query_clade, terminal)
                if distance < min_distance:
                    min_distance = distance
                    nearest = terminal

        return nearest

    def run_complete_workflow(self, sanger_fasta: Path, run_id: str,
                            sample_id: str) -> Dict:
        """
        Run complete Spumavirus detection workflow.

        Args:
            sanger_fasta: Sanger sequencing FASTA file (~400 bp)
            run_id: Run identifier
            sample_id: Sample identifier (PBMC source)

        Returns:
            Complete analysis results
        """
        print("=" * 80)
        print("SPUMAVIRUS NESTED PCR DETECTION - Complete Workflow")
        print("=" * 80)
        print(f"Run ID: {run_id}")
        print(f"Sample ID: {sample_id}")
        print(f"Input sequence: {sanger_fasta}")
        print()

        results = {
            'run_id': run_id,
            'sample_id': sample_id,
            'input_file': str(sanger_fasta),
            'workflow_steps': {}
        }

        # Step 1: Validate sequence
        print("STEP 1: Sequence Validation")
        print("-" * 80)
        is_valid, seq_len = self.validate_sequence_length(sanger_fasta, min_length=300)

        results['workflow_steps']['validation'] = {
            'sequence_length': seq_len,
            'valid': is_valid,
            'expected_length': '350-450 bp'
        }

        if not is_valid:
            print(f"ERROR: Sequence too short ({seq_len} bp). Expected 300-450 bp.", file=sys.stderr)
            results['final_determination'] = 'invalid_sequence'
            return results

        print(f"Sequence length: {seq_len} bp - VALID")

        # Step 2: BLAST analysis
        print("\nSTEP 2: BLAST Analysis")
        print("-" * 80)
        blast_xml = self.work_dir / "blast_results.xml"
        blast_results = self.blast_sequence(sanger_fasta, blast_xml)

        results['workflow_steps']['blast'] = blast_results

        if blast_results['hits']:
            print(f"\nTop 3 BLAST hits:")
            for i, hit in enumerate(blast_results['hits'][:3], 1):
                print(f"  {i}. {hit['title'][:80]}")
                print(f"     Identity: {hit['identity']:.1f}%, E-value: {hit['e_value']:.2e}")

        # Step 3: Classification
        print("\nSTEP 3: Virus Classification")
        print("-" * 80)
        classification = self.classify_virus(blast_results)

        results['workflow_steps']['classification'] = classification

        print(f"Classification: {classification['classification'].upper()}")
        print(f"Confidence: {classification['confidence']}")
        print(f"Reason: {classification['reason']}")

        if 'warning' in classification:
            print(f"\nWARNING: {classification['warning']}")

        # Step 4: Phylogenetic analysis (if ambiguous)
        if classification['confidence'] == 'low' or classification['classification'] == 'ambiguous':
            print("\nSTEP 4: Phylogenetic Analysis (disambiguation)")
            print("-" * 80)

            ref_dir = self.work_dir / "references"
            ref_fasta = self.download_reference_sequences(ref_dir)

            tree_file = self.work_dir / "phylogenetic_tree.nwk"
            tree_results = self.build_phylogenetic_tree(sanger_fasta, ref_fasta, tree_file)

            results['workflow_steps']['phylogenetic_analysis'] = tree_results

            if 'error' not in tree_results:
                print(f"Nearest neighbor: {tree_results['nearest_neighbor']}")
                print(f"Clustering: {tree_results['clustering']}")

                # Update classification based on phylogeny
                if tree_results['clustering'] == 'Spumavirus':
                    results['final_determination'] = 'spumavirus'
                    results['confidence'] = 'medium'
                else:
                    results['final_determination'] = 'PERV'
                    results['confidence'] = 'high'
            else:
                results['final_determination'] = 'indeterminate'
                results['confidence'] = 'none'
        else:
            # High confidence classification
            results['final_determination'] = classification['classification']
            results['confidence'] = classification['confidence']

        # Final report
        print("\n" + "=" * 80)
        print("FINAL DETERMINATION")
        print("=" * 80)
        print(f"Virus: {results['final_determination'].upper()}")
        print(f"Confidence: {results['confidence']}")

        if results['final_determination'] == 'spumavirus':
            print("\n⚠️  POSITIVE DETECTION: Porcine Spumavirus")
            print("Action required: Report to PMDA, reject donor")
        elif results['final_determination'] == 'PERV':
            print("\n✓ PERV detected (endogenous, expected)")
            print("Action: No action required - PERV is normal in all pigs")
        else:
            print("\n? Indeterminate result")
            print("Action: Repeat nested PCR, sequence additional clones")

        return results


def main():
    parser = argparse.ArgumentParser(
        description='Spumavirus nested PCR detection and confirmation workflow',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
1. Nested PCR from PBMC genomic DNA
2. Agarose gel electrophoresis (~400 bp band)
3. Sanger sequencing (both strands, consensus)
4. BLAST analysis (this script)
5. Phylogenetic analysis (this script)
6. Spumavirus vs PERV discrimination (this script)

Example:
  %(prog)s -i sanger_sequence.fasta -r RUN-001 -s PIG-2024-042-PBMC
        """
    )

    parser.add_argument('-i', '--input', required=True, type=Path,
                       help='Sanger sequence FASTA file (~400 bp)')
    parser.add_argument('-r', '--run-id', required=True,
                       help='Run identifier')
    parser.add_argument('-s', '--sample-id', required=True,
                       help='Sample identifier (PBMC source)')
    parser.add_argument('-o', '--output', required=True, type=Path,
                       help='Output JSON file')
    parser.add_argument('-w', '--work-dir', type=Path,
                       default=Path('./spumavirus_analysis'),
                       help='Working directory for intermediate files')
    parser.add_argument('--blast-db', type=Path,
                       help='Local BLAST database (optional, uses NCBI if not specified)')

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Initialize detector
    detector = SpumavirusDetector(work_dir=args.work_dir, blast_db=args.blast_db)

    # Run workflow
    try:
        results = detector.run_complete_workflow(args.input, args.run_id, args.sample_id)

        # Save results
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"\nResults saved to: {args.output}")

        # Exit code based on determination
        if results['final_determination'] == 'spumavirus':
            sys.exit(1)  # Positive detection - alert exit code
        else:
            sys.exit(0)  # Negative or PERV - normal exit

    except Exception as e:
        print(f"\nERROR: Workflow failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(2)


if __name__ == '__main__':
    main()
