#!/usr/bin/env python3
"""
Unit tests for pathogen detection module
"""

import unittest
import tempfile
import json
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd

class TestPathogenDetection(unittest.TestCase):
    """Tests for pathogen detection phase operations."""

    def setUp(self):
        """Set up test environment."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.sample_fastq = self.test_dir / 'sample.fastq'

        # Create sample FASTQ data
        with open(self.sample_fastq, 'w') as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCGATCG\n")
            f.write("+\n")
            f.write("IIIIIIIIIIIIIIII\n")

        # PMDA pathogens for testing
        self.pmda_pathogens = {
            'PERV-A': {'name': 'Porcine Endogenous Retrovirus A', 'risk': 'CRITICAL'},
            'PCV2': {'name': 'Porcine circovirus 2', 'risk': 'MEDIUM'},
            'HEV': {'name': 'Hepatitis E virus', 'risk': 'HIGH'}
        }

    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_kraken2_classification(self):
        """Test Kraken2 taxonomic classification."""
        from scripts.phase4_pathogen import kraken2_classifier

        # Mock Kraken2 output
        mock_output = """
        U\tread1\t0\t100\t0:66
        C\tread2\t10239\t150\t10239:116
        C\tread3\t12345\t200\t12345:166
        """

        with patch('scripts.phase4_pathogen.kraken2_classifier.run_kraken2') as mock_kraken:
            mock_kraken.return_value = mock_output

            results = kraken2_classifier.classify(
                input_fastq=str(self.sample_fastq),
                database='/data/kraken2_db'
            )

            self.assertIn('classified_reads', results)
            self.assertIn('unclassified_reads', results)
            self.assertEqual(results['classified_reads'], 2)
            self.assertEqual(results['unclassified_reads'], 1)

    def test_perv_detection(self):
        """Test PERV-specific detection."""
        from scripts.phase4_pathogen import perv_detector

        # Mock BLAST results with PERV hit
        mock_blast_results = [
            {
                'query': 'read1',
                'subject': 'PERV-A_ref',
                'identity': 98.5,
                'length': 1500,
                'evalue': 1e-100
            }
        ]

        with patch('scripts.phase4_pathogen.perv_detector.run_blast') as mock_blast:
            mock_blast.return_value = mock_blast_results

            perv_results = perv_detector.detect_perv(
                input_fastq=str(self.sample_fastq),
                perv_database='/data/perv_db'
            )

            self.assertTrue(perv_results['perv_detected'])
            self.assertIn('PERV-A', perv_results['variants_found'])
            self.assertEqual(perv_results['risk_level'], 'CRITICAL')

    def test_rvdb_virus_search(self):
        """Test RVDB viral database search."""
        from scripts.phase4_pathogen import rvdb_searcher

        # Mock Diamond output
        mock_diamond_output = pd.DataFrame({
            'qseqid': ['read1', 'read2'],
            'sseqid': ['NC_001422.1|Coliphage_phiX174', 'NC_002549.1|Zaire_ebolavirus'],
            'pident': [99.0, 85.0],
            'length': [100, 150],
            'evalue': [1e-50, 1e-30]
        })

        with patch('scripts.phase4_pathogen.rvdb_searcher.run_diamond') as mock_diamond:
            mock_diamond.return_value = mock_diamond_output

            virus_results = rvdb_searcher.search_viruses(
                input_fastq=str(self.sample_fastq),
                rvdb_path='/data/rvdb.faa'
            )

            self.assertEqual(len(virus_results['viruses_found']), 2)
            self.assertIn('phiX174', str(virus_results['viruses_found']))

    def test_pmda_pathogen_check(self):
        """Test PMDA 91 pathogen checklist validation."""
        from scripts.phase4_pathogen import pmda_validator

        # Mock detection results
        detected_pathogens = {
            'PERV-A': {'reads': 10, 'confidence': 0.99},
            'PCV2': {'reads': 50, 'confidence': 0.95},
            'HEV': {'reads': 5, 'confidence': 0.85}
        }

        # Run validation
        pmda_report = pmda_validator.generate_pmda_report(
            detected_pathogens=detected_pathogens,
            pmda_list=self.pmda_pathogens
        )

        self.assertEqual(pmda_report['total_pmda_pathogens'], 91)
        self.assertEqual(pmda_report['detected_count'], 3)
        self.assertIn('PERV-A', pmda_report['critical_pathogens_detected'])
        self.assertTrue(pmda_report['immediate_alert_required'])

    def test_quantification_with_spike_in(self):
        """Test absolute quantification using spike-in control."""
        from scripts.phase5_quantification import spike_in_quantifier

        # Mock data
        spike_in_reads = 1000
        spike_in_copies = 1e6  # Known concentration
        pathogen_reads = 50

        # Calculate copies/mL
        pathogen_copies = spike_in_quantifier.calculate_absolute_copies(
            pathogen_reads=pathogen_reads,
            spike_in_reads=spike_in_reads,
            spike_in_copies_per_ml=spike_in_copies
        )

        expected = (pathogen_reads / spike_in_reads) * spike_in_copies
        self.assertEqual(pathogen_copies, expected)
        self.assertEqual(pathogen_copies, 5e4)

    def test_confidence_score_calculation(self):
        """Test pathogen detection confidence scoring."""
        from scripts.phase4_pathogen import confidence_scorer

        # Test high confidence (multiple methods agree)
        detections = {
            'kraken2': {'PCV2': 100},
            'blast': {'PCV2': 95},
            'diamond': {'PCV2': 90}
        }

        confidence = confidence_scorer.calculate_confidence(
            pathogen='PCV2',
            detections=detections
        )

        self.assertGreater(confidence, 0.9)

        # Test low confidence (single method)
        detections_low = {
            'kraken2': {'PCV2': 10}
        }

        confidence_low = confidence_scorer.calculate_confidence(
            pathogen='PCV2',
            detections=detections_low
        )

        self.assertLess(confidence_low, 0.5)

    def test_host_contamination_detection(self):
        """Test detection of host genome contamination."""
        from scripts.phase3_host_removal import contamination_checker

        # Mock alignment stats
        alignment_stats = {
            'total_reads': 100000,
            'aligned_reads': 95000,
            'unaligned_reads': 5000
        }

        contamination_level = contamination_checker.calculate_contamination(
            alignment_stats
        )

        self.assertEqual(contamination_level, 0.95)

        # Check if depletion is insufficient
        is_sufficient = contamination_checker.check_depletion_efficiency(
            contamination_level,
            threshold=0.90
        )

        self.assertFalse(is_sufficient)  # 95% contamination means only 5% depletion

    def test_multi_database_integration(self):
        """Test integration of results from multiple databases."""
        from scripts.phase4_pathogen import result_integrator

        # Mock results from different databases
        kraken_results = {
            'PCV2': {'reads': 100, 'confidence': 0.95},
            'PRRSV': {'reads': 50, 'confidence': 0.80}
        }

        blast_results = {
            'PCV2': {'reads': 95, 'confidence': 0.98},
            'HEV': {'reads': 20, 'confidence': 0.75}
        }

        rvdb_results = {
            'PCV2': {'reads': 90, 'confidence': 0.92},
            'PERV-A': {'reads': 5, 'confidence': 0.99}
        }

        # Integrate results
        integrated = result_integrator.integrate_all_results(
            kraken=kraken_results,
            blast=blast_results,
            rvdb=rvdb_results
        )

        # PCV2 should have highest confidence (detected by all 3)
        self.assertIn('PCV2', integrated)
        self.assertGreater(integrated['PCV2']['confidence'], 0.95)

        # PERV-A should trigger critical alert despite low reads
        self.assertIn('PERV-A', integrated)
        self.assertEqual(integrated['PERV-A']['risk_level'], 'CRITICAL')

    def test_unknown_pathogen_detection(self):
        """Test detection of unknown/novel pathogens."""
        from scripts.phase4_pathogen import novel_pathogen_detector

        # Mock unclassified reads with high similarity cluster
        unclassified_reads = [
            {'id': 'read1', 'seq': 'ATCGATCGATCG' * 10},
            {'id': 'read2', 'seq': 'ATCGATCGATCG' * 10},  # Similar
            {'id': 'read3', 'seq': 'GCTAGCTAGCTA' * 10}   # Different
        ]

        # Detect potential novel pathogen
        novel_clusters = novel_pathogen_detector.find_novel_clusters(
            unclassified_reads,
            min_cluster_size=2,
            similarity_threshold=0.90
        )

        self.assertGreaterEqual(len(novel_clusters), 1)
        self.assertIn('potential_novel_pathogen', novel_clusters[0])

    def test_detection_limit_validation(self):
        """Test limit of detection (LOD) for pathogen detection."""
        from scripts.phase4_pathogen import lod_validator

        # Test detection at different dilutions
        dilutions = [1e6, 1e5, 1e4, 1e3, 1e2, 1e1]
        detected = [True, True, True, True, False, False]

        lod = lod_validator.calculate_lod(
            dilutions=dilutions,
            detected=detected,
            confidence_level=0.95
        )

        self.assertEqual(lod, 1e3)  # LOD is lowest detected concentration

    def test_critical_pathogen_alert(self):
        """Test immediate alert for critical pathogen detection."""
        from scripts.phase4_pathogen import alert_manager

        with patch('scripts.phase4_pathogen.alert_manager.send_sns_alert') as mock_sns:
            # Detect PERV
            alert_manager.check_and_alert(
                pathogen='PERV-A',
                risk_level='CRITICAL',
                run_id='TEST-001'
            )

            # Verify immediate SNS alert was sent
            mock_sns.assert_called_once()
            call_args = mock_sns.call_args[0]
            self.assertIn('CRITICAL', call_args[0])  # Subject
            self.assertIn('PERV-A', call_args[1])    # Message

    def test_pmda_compliance_validation(self):
        """Test full PMDA 91 pathogen compliance check."""
        from scripts.phase4_pathogen import pmda_compliance

        # Mock detection covering all categories
        detection_results = {
            'viruses': ['PERV-A', 'HEV', 'JEV'],  # 3/50
            'bacteria': ['MT', 'BA'],              # 2/35
            'parasites': ['TG'],                   # 1/5
            'fungi': [],                           # 0/5
            'prions': []                           # 0/1
        }

        compliance_report = pmda_compliance.validate_coverage(
            detection_results,
            pmda_pathogen_list=self.pmda_pathogens
        )

        self.assertEqual(compliance_report['total_pmda'], 91)
        self.assertEqual(compliance_report['total_detected'], 6)
        self.assertTrue(compliance_report['perv_analysis_included'])
        self.assertIn('coverage_by_category', compliance_report)

if __name__ == '__main__':
    unittest.main()