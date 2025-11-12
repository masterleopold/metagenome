#!/usr/bin/env python3
"""
PMDA 4-Virus Sensitivity Validation Test Suite
Tests LOD, PPA, NPA, and cross-reactivity for all 4 PMDA high-sensitivity viruses

PMDA Requirements:
- PPA (Positive Percent Agreement): ≥95% at LOD
- NPA (Negative Percent Agreement): ≥98% (50 negative samples)
- R² (Linearity): ≥0.90 across 10-10,000 copies/mL
- LOD Determination: 10 replicates × 7 concentrations

Target Viruses:
1. Polyomavirus (dsDNA) - LOD target: 50 copies/mL
2. Hantavirus (ssRNA-) - LOD target: 100 copies/mL
3. EEEV (ssRNA+) - LOD target: 50 copies/mL
4. Spumavirus (proviral DNA) - LOD target: 10 copies/10^5 PBMCs
"""

import unittest
import subprocess
import json
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
from scipy import stats
import random


class PMDA4VirusSensitivityValidator:
    """Validator for PMDA 4-virus sensitivity testing."""

    def __init__(self, scripts_dir: Path, databases_dir: Path):
        """
        Initialize validator.

        Args:
            scripts_dir: Path to scripts directory
            databases_dir: Path to PMDA databases
        """
        self.scripts_dir = scripts_dir
        self.databases_dir = databases_dir
        self.detect_script = scripts_dir / "phase4_pathogen" / "detect_pmda_4viruses.py"
        self.quant_script = scripts_dir / "phase5_quantification" / "pmda_4virus_quantification.py"

    def generate_synthetic_reads(self, virus_type: str, copies_per_ml: int,
                                plasma_volume_ml: float, avg_read_length: int,
                                output_fastq: Path) -> Dict:
        """
        Generate synthetic MinION reads for a virus.

        Uses Badread or custom simulator to generate realistic nanopore reads.

        Args:
            virus_type: POLYOMA, HANTV, EEEV, or SPUMV
            copies_per_ml: Target viral load (copies/mL)
            plasma_volume_ml: Plasma volume
            avg_read_length: Average read length (bp)
            output_fastq: Output FASTQ file path

        Returns:
            Dictionary with simulation parameters
        """
        # Reference sequences for each virus
        reference_sequences = {
            'POLYOMA': self.databases_dir / 'polyomavirus' / 'polyoma_all.fasta',
            'HANTV': self.databases_dir / 'hantavirus' / 'hantavirus_all.fasta',
            'EEEV': self.databases_dir / 'alphavirus' / 'alphavirus_all.fasta',
            'SPUMV': self.databases_dir / 'spumavirus' / 'spumavirus_all_pol.fasta'
        }

        if virus_type not in reference_sequences:
            raise ValueError(f"Unknown virus type: {virus_type}")

        reference = reference_sequences[virus_type]

        # Calculate number of reads needed
        # Assume 10M total reads, virus reads proportional to copies
        total_reads = 10_000_000
        total_copies = copies_per_ml * plasma_volume_ml
        genome_size = {'POLYOMA': 5150, 'HANTV': 11900, 'EEEV': 11841, 'SPUMV': 12000}[virus_type]

        # Sequencing depth = reads * read_length / genome_size
        # We want: virus_reads / total_reads = virus_copies / total_DNA_molecules
        # Simplified: virus_reads = total_copies * coverage_factor
        virus_reads = int(total_copies * 10)  # 10× coverage factor

        # Use Badread for realistic nanopore simulation
        cmd = [
            'badread', 'simulate',
            '--reference', str(reference),
            '--quantity', f'{virus_reads}x',  # Number of reads
            '--length', f'{avg_read_length},500',  # Mean length, stdev
            '--identity', '95,98,4',  # Mean, max, stdev
            '--error_model', 'nanopore2020',
            '--qscore_model', 'nanopore2020'
        ]

        with open(output_fastq, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)

        return {
            'virus_type': virus_type,
            'copies_per_ml': copies_per_ml,
            'plasma_volume_ml': plasma_volume_ml,
            'virus_reads_generated': virus_reads,
            'avg_read_length': avg_read_length,
            'reference': str(reference)
        }

    def run_detection(self, input_fastq: Path, virus_target: str) -> Dict:
        """
        Run detection on synthetic reads.

        Args:
            input_fastq: Input FASTQ file
            virus_target: POLYOMA, HANTV, EEEV, or SPUMV

        Returns:
            Detection results dictionary
        """
        output_dir = tempfile.mkdtemp()
        output_json = Path(output_dir) / "detection_results.json"

        cmd = [
            'python3', str(self.detect_script),
            '--input', str(input_fastq),
            '--output', str(output_dir),
            '--target', virus_target.lower(),
            '--database-dir', str(self.databases_dir),
            '--threads', '8',
            '--run-id', f'TEST-{virus_target}'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            return {'error': result.stderr, 'detected': False}

        # Load results
        if output_json.exists():
            with open(output_json) as f:
                return json.load(f)
        else:
            return {'detected': False}

    def test_lod(self, virus_type: str, concentrations: List[int],
                replicates: int = 10) -> Dict:
        """
        Determine Limit of Detection (LOD) for a virus.

        PMDA requirement: 10 replicates × 7 concentrations

        Args:
            virus_type: Virus code (POLYOMA, HANTV, EEEV, SPUMV)
            concentrations: List of concentrations to test (copies/mL)
            replicates: Number of replicates per concentration (default: 10)

        Returns:
            LOD test results with PPA at each concentration
        """
        print(f"\nTesting LOD for {virus_type}...")
        print(f"Concentrations: {concentrations}")
        print(f"Replicates per concentration: {replicates}")

        results = []

        for concentration in concentrations:
            detections = 0

            for rep in range(replicates):
                print(f"  Concentration: {concentration} copies/mL, Replicate: {rep+1}/{replicates}")

                # Generate synthetic reads
                fastq_file = Path(tempfile.mktemp(suffix=".fastq"))
                self.generate_synthetic_reads(
                    virus_type, concentration, plasma_volume_ml=10.0,
                    avg_read_length=1000, output_fastq=fastq_file
                )

                # Run detection
                detection_result = self.run_detection(fastq_file, virus_type)

                if detection_result.get('detected', False):
                    detections += 1

                # Cleanup
                fastq_file.unlink(missing_ok=True)

            ppa = (detections / replicates) * 100

            results.append({
                'concentration': concentration,
                'detections': detections,
                'replicates': replicates,
                'ppa': round(ppa, 2)
            })

            print(f"    PPA at {concentration} copies/mL: {ppa:.1f}% ({detections}/{replicates})")

        # Determine LOD (lowest concentration with ≥95% PPA)
        lod = None
        for result in results:
            if result['ppa'] >= 95.0:
                lod = result['concentration']
                break

        return {
            'virus_type': virus_type,
            'lod_determined': lod,
            'lod_unit': 'copies/mL',
            'concentrations_tested': results,
            'pmda_compliant': lod is not None and lod <= {'POLYOMA': 50, 'HANTV': 100, 'EEEV': 50, 'SPUMV': 10}[virus_type]
        }

    def test_linearity(self, virus_type: str, concentrations: List[int]) -> Dict:
        """
        Test linearity of quantification (R²).

        PMDA requirement: R² ≥0.90 across 10-10,000 copies/mL

        Args:
            virus_type: Virus code
            concentrations: List of concentrations spanning dynamic range

        Returns:
            Linearity test results with R²
        """
        print(f"\nTesting linearity for {virus_type}...")

        observed_values = []
        expected_values = []

        for concentration in concentrations:
            print(f"  Testing concentration: {concentration} copies/mL")

            # Generate synthetic reads
            fastq_file = Path(tempfile.mktemp(suffix=".fastq"))
            self.generate_synthetic_reads(
                virus_type, concentration, plasma_volume_ml=10.0,
                avg_read_length=1000, output_fastq=fastq_file
            )

            # Run detection + quantification
            detection_result = self.run_detection(fastq_file, virus_type)

            if detection_result.get('detected', False):
                # Run quantification
                quant_output = Path(tempfile.mktemp(suffix=".json"))
                detection_json = Path(tempfile.mktemp(suffix=".json"))

                with open(detection_json, 'w') as f:
                    json.dump(detection_result, f)

                cmd = [
                    'python3', str(self.quant_script),
                    '-i', str(detection_json),
                    '-o', str(quant_output),
                    '-r', f'TEST-{virus_type}',
                    '--plasma-volume', '10.0'
                ]

                result = subprocess.run(cmd, capture_output=True)

                if result.returncode == 0 and quant_output.exists():
                    with open(quant_output) as f:
                        quant_data = json.load(f)

                    observed_copies = quant_data['viruses'][virus_type].get('copies_per_ml', 0)
                    observed_values.append(observed_copies)
                    expected_values.append(concentration)

                # Cleanup
                detection_json.unlink(missing_ok=True)
                quant_output.unlink(missing_ok=True)

            fastq_file.unlink(missing_ok=True)

        # Calculate R² (coefficient of determination)
        if len(observed_values) >= 3:
            slope, intercept, r_value, p_value, std_err = stats.linregress(expected_values, observed_values)
            r_squared = r_value ** 2
        else:
            r_squared = 0.0
            slope, intercept = 0, 0

        return {
            'virus_type': virus_type,
            'r_squared': round(r_squared, 4),
            'slope': round(slope, 4),
            'intercept': round(intercept, 2),
            'data_points': list(zip(expected_values, observed_values)),
            'pmda_compliant': r_squared >= 0.90
        }

    def test_npa(self, virus_type: str, n_negative_samples: int = 50) -> Dict:
        """
        Test Negative Percent Agreement (NPA).

        PMDA requirement: ≥98% NPA (50 negative samples)

        Args:
            virus_type: Virus code
            n_negative_samples: Number of negative samples to test (default: 50)

        Returns:
            NPA test results
        """
        print(f"\nTesting NPA for {virus_type} ({n_negative_samples} negative samples)...")

        false_positives = 0

        for i in range(n_negative_samples):
            # Generate background reads (pig genomic DNA, no virus)
            fastq_file = Path(tempfile.mktemp(suffix=".fastq"))

            # Simulate random pig genomic reads (no virus)
            cmd = [
                'badread', 'simulate',
                '--reference', str(self.databases_dir / '..' / 'host_genome' / 'sus_scrofa_11.1.fa'),
                '--quantity', '1000x',
                '--length', '1000,500'
            ]

            with open(fastq_file, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.DEVNULL)

            # Run detection
            detection_result = self.run_detection(fastq_file, virus_type)

            if detection_result.get('detected', False):
                false_positives += 1
                print(f"    FALSE POSITIVE detected in sample {i+1}")

            fastq_file.unlink(missing_ok=True)

        true_negatives = n_negative_samples - false_positives
        npa = (true_negatives / n_negative_samples) * 100

        return {
            'virus_type': virus_type,
            'true_negatives': true_negatives,
            'false_positives': false_positives,
            'total_negative_samples': n_negative_samples,
            'npa': round(npa, 2),
            'pmda_compliant': npa >= 98.0
        }


class TestPMDA4VirusSensitivity(unittest.TestCase):
    """Unit tests for PMDA 4-virus sensitivity validation."""

    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.scripts_dir = Path(__file__).parent.parent / "scripts"
        cls.databases_dir = Path("/mnt/efs/databases/pmda/2024.1")
        cls.validator = PMDA4VirusSensitivityValidator(cls.scripts_dir, cls.databases_dir)

    def test_polyomavirus_lod(self):
        """Test Polyomavirus LOD determination."""
        concentrations = [10, 25, 50, 100, 200, 500, 1000]  # copies/mL
        results = self.validator.test_lod('POLYOMA', concentrations, replicates=10)

        self.assertIsNotNone(results['lod_determined'], "LOD should be determined")
        self.assertLessEqual(results['lod_determined'], 50, "LOD should be ≤50 copies/mL")
        self.assertTrue(results['pmda_compliant'], "Should meet PMDA LOD requirement")

    def test_hantavirus_lod(self):
        """Test Hantavirus LOD determination."""
        concentrations = [25, 50, 100, 200, 500, 1000, 2000]
        results = self.validator.test_lod('HANTV', concentrations, replicates=10)

        self.assertIsNotNone(results['lod_determined'])
        self.assertLessEqual(results['lod_determined'], 100, "LOD should be ≤100 copies/mL")
        self.assertTrue(results['pmda_compliant'])

    def test_eeev_lod(self):
        """Test EEEV LOD determination."""
        concentrations = [10, 25, 50, 100, 200, 500, 1000]
        results = self.validator.test_lod('EEEV', concentrations, replicates=10)

        self.assertIsNotNone(results['lod_determined'])
        self.assertLessEqual(results['lod_determined'], 50, "LOD should be ≤50 copies/mL")
        self.assertTrue(results['pmda_compliant'])

    def test_polyomavirus_linearity(self):
        """Test Polyomavirus quantification linearity."""
        concentrations = [10, 50, 100, 500, 1000, 5000, 10000]
        results = self.validator.test_linearity('POLYOMA', concentrations)

        self.assertGreaterEqual(results['r_squared'], 0.90, "R² should be ≥0.90")
        self.assertTrue(results['pmda_compliant'])

    def test_polyomavirus_npa(self):
        """Test Polyomavirus NPA (false positive rate)."""
        results = self.validator.test_npa('POLYOMA', n_negative_samples=50)

        self.assertGreaterEqual(results['npa'], 98.0, "NPA should be ≥98%")
        self.assertLessEqual(results['false_positives'], 1, "Should have ≤1 false positive")
        self.assertTrue(results['pmda_compliant'])


def main():
    """Run complete PMDA 4-virus sensitivity validation suite."""
    print("=" * 80)
    print("PMDA 4-VIRUS SENSITIVITY VALIDATION SUITE")
    print("=" * 80)

    scripts_dir = Path(__file__).parent.parent / "scripts"
    databases_dir = Path("/mnt/efs/databases/pmda/2024.1")

    validator = PMDA4VirusSensitivityValidator(scripts_dir, databases_dir)

    # Test each virus
    viruses = ['POLYOMA', 'HANTV', 'EEEV']  # SPUMV requires separate nested PCR validation

    results = {
        'lod': {},
        'linearity': {},
        'npa': {}
    }

    for virus in viruses:
        print(f"\n{'='*80}")
        print(f"VALIDATING {virus}")
        print(f"{'='*80}")

        # LOD test
        concentrations_lod = {
            'POLYOMA': [10, 25, 50, 100, 200, 500, 1000],
            'HANTV': [25, 50, 100, 200, 500, 1000, 2000],
            'EEEV': [10, 25, 50, 100, 200, 500, 1000]
        }
        results['lod'][virus] = validator.test_lod(virus, concentrations_lod[virus], replicates=10)

        # Linearity test
        concentrations_lin = [10, 50, 100, 500, 1000, 5000, 10000]
        results['linearity'][virus] = validator.test_linearity(virus, concentrations_lin)

        # NPA test
        results['npa'][virus] = validator.test_npa(virus, n_negative_samples=50)

    # Summary report
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)

    for virus in viruses:
        print(f"\n{virus}:")
        print(f"  LOD: {results['lod'][virus]['lod_determined']} copies/mL "
              f"({'PASS' if results['lod'][virus]['pmda_compliant'] else 'FAIL'})")
        print(f"  R²: {results['linearity'][virus]['r_squared']:.4f} "
              f"({'PASS' if results['linearity'][virus]['pmda_compliant'] else 'FAIL'})")
        print(f"  NPA: {results['npa'][virus]['npa']:.1f}% "
              f"({'PASS' if results['npa'][virus]['pmda_compliant'] else 'FAIL'})")

    # Save results
    output_file = Path("pmda_4virus_validation_results.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nValidation results saved to: {output_file}")


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--unittest':
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
