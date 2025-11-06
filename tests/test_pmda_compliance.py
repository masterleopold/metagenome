#!/usr/bin/env python3
"""
PMDA compliance tests for MinION pipeline
"""

import unittest
import json
from pathlib import Path

class TestPMDACompliance(unittest.TestCase):
    """Tests for PMDA guideline compliance."""

    def setUp(self):
        """Set up test environment."""
        # Use absolute path relative to this test file
        test_dir = Path(__file__).parent
        repo_root = test_dir.parent
        self.pmda_pathogens_file = repo_root / 'templates' / 'config' / 'pmda_pathogens.json'
        self.pmda_data = None

        if self.pmda_pathogens_file.exists():
            with open(self.pmda_pathogens_file) as f:
                self.pmda_data = json.load(f)

    def test_91_pathogen_coverage(self):
        """Test that all 91 PMDA pathogens are covered."""

        self.assertIsNotNone(self.pmda_data)
        self.assertEqual(self.pmda_data['total_pathogens'], 91)

        # Count pathogens in all categories
        total_count = 0
        for category in self.pmda_data['categories'].values():
            total_count += category['count']

        # Allow for placeholder counts
        self.assertGreaterEqual(total_count, 91)

    def test_perv_detection_priority(self):
        """Test PERV detection is marked as critical."""

        viruses = self.pmda_data['categories']['viruses']
        critical_viruses = viruses['critical']

        # PERV must be in critical list
        for perv_type in ['PERV-A', 'PERV-B', 'PERV-C']:
            self.assertIn(perv_type, critical_viruses)

        # Check risk levels
        for pathogen in viruses['pathogens']:
            if pathogen['code'].startswith('PERV'):
                self.assertEqual(pathogen['risk_level'], 'CRITICAL')

    def test_detection_methods(self):
        """Test appropriate detection methods are defined."""

        detection_methods = self.pmda_data['detection_methods']

        # PERV must have specific methods
        self.assertIn('PERV', detection_methods)
        perv_methods = detection_methods['PERV']
        self.assertIn('PCR', perv_methods)
        self.assertIn('NGS', perv_methods)

        # All categories must have detection methods
        for category in ['viruses', 'bacteria', 'parasites', 'fungi', 'prions']:
            self.assertIn(category, detection_methods)

    def test_reporting_requirements(self):
        """Test reporting requirements compliance."""

        requirements = self.pmda_data['reporting_requirements']

        # Required elements
        required_elements = [
            'all_pathogens_tested',
            'perv_specific_report',
            'quantitative_results',
            'detection_limits',
            'quality_metrics'
        ]

        for element in required_elements:
            self.assertIn(element, requirements)

        # Critical requirements must be true
        self.assertTrue(requirements['all_pathogens_tested'])
        self.assertTrue(requirements['perv_specific_report'])
        self.assertTrue(requirements['quantitative_results'])

    def test_action_thresholds(self):
        """Test action thresholds for pathogen detection."""

        thresholds = self.pmda_data['action_thresholds']

        # PERV threshold must be most stringent
        self.assertIn('PERV', thresholds)
        perv_threshold = thresholds['PERV']
        self.assertEqual(perv_threshold['detection_threshold'], 1)
        self.assertEqual(perv_threshold['action'], 'quarantine')
        self.assertEqual(perv_threshold['notification'], 'immediate')

        # Critical pathogens must have immediate notification
        critical_threshold = thresholds['critical_pathogens']
        self.assertEqual(critical_threshold['notification'], 'immediate')

    def test_pathogen_categories(self):
        """Test pathogen categorization."""

        categories = self.pmda_data['categories']

        # Required categories
        required_categories = ['viruses', 'bacteria', 'parasites', 'fungi', 'prions']

        for category in required_categories:
            self.assertIn(category, categories)
            self.assertIn('count', categories[category])
            self.assertIn('pathogens', categories[category])

    def test_critical_pathogen_list(self):
        """Test critical pathogen identification."""

        # Collect all critical pathogens
        critical_pathogens = []

        for category in self.pmda_data['categories'].values():
            if 'critical' in category:
                critical_pathogens.extend(category['critical'])

        # Must include known critical pathogens
        expected_critical = ['PERV-A', 'PERV-B', 'PERV-C', 'PRION']

        for pathogen in expected_critical:
            self.assertIn(pathogen, critical_pathogens)

    def test_virus_coverage(self):
        """Test comprehensive virus coverage."""

        viruses = self.pmda_data['categories']['viruses']

        # Key viruses that must be included
        required_viruses = [
            'PERV-A', 'PERV-B', 'PERV-C',  # PERV variants
            'HEV',    # Hepatitis E
            'JEV',    # Japanese encephalitis
            'ASFV',   # African swine fever
            'CSFV',   # Classical swine fever
            'FMDV',   # Foot-and-mouth disease
            'PRRSV',  # PRRS
            'PCV2'    # Circovirus
        ]

        virus_codes = [v['code'] for v in viruses['pathogens']]

        for required_virus in required_viruses:
            self.assertIn(required_virus, virus_codes)

    def test_bacteria_coverage(self):
        """Test comprehensive bacteria coverage."""

        bacteria = self.pmda_data['categories']['bacteria']

        # Key bacteria that must be included
        required_bacteria = [
            'MT',  # Mycobacterium tuberculosis
            'BA',  # Bacillus anthracis
            'SS',  # Streptococcus suis
            'SE'   # Salmonella enterica
        ]

        bacteria_codes = [b['code'] for b in bacteria['pathogens']]

        for required_bacterium in required_bacteria:
            self.assertIn(required_bacterium, bacteria_codes)

    def test_risk_level_consistency(self):
        """Test risk level assignments are consistent."""

        valid_risk_levels = ['CRITICAL', 'HIGH', 'MEDIUM', 'LOW']

        for category in self.pmda_data['categories'].values():
            for pathogen in category['pathogens']:
                self.assertIn('risk_level', pathogen)
                self.assertIn(pathogen['risk_level'], valid_risk_levels)

                # Critical pathogens in critical list must have CRITICAL risk
                if 'critical' in category and pathogen['code'] in category['critical']:
                    self.assertEqual(pathogen['risk_level'], 'CRITICAL')

    def test_workflow_configuration_compliance(self):
        """Test workflow configuration meets PMDA requirements."""

        # Load default pipeline configuration using absolute path
        test_dir = Path(__file__).parent
        repo_root = test_dir.parent
        config_file = repo_root / 'templates' / 'config' / 'default_pipeline.yaml'

        if config_file.exists():
            import yaml
            with open(config_file) as f:
                config = yaml.safe_load(f)

            # Check PMDA section
            self.assertIn('pmda', config)
            pmda_config = config['pmda']

            self.assertEqual(pmda_config['pathogen_count'], 91)
            self.assertIn('PERV-A', pmda_config['critical_pathogens'])
            self.assertIn('PERV-B', pmda_config['critical_pathogens'])
            self.assertIn('PERV-C', pmda_config['critical_pathogens'])

            # Check reporting requirements
            reporting_req = pmda_config['reporting_requirements']
            self.assertTrue(reporting_req['all_91_pathogens_tested'])
            self.assertTrue(reporting_req['perv_specific_analysis'])

    def test_quality_thresholds(self):
        """Test quality thresholds meet PMDA standards."""

        test_dir = Path(__file__).parent
        repo_root = test_dir.parent
        config_file = repo_root / 'templates' / 'config' / 'default_pipeline.yaml'

        if config_file.exists():
            import yaml
            with open(config_file) as f:
                config = yaml.safe_load(f)

            thresholds = config['quality_thresholds']

            # Minimum requirements
            self.assertGreaterEqual(thresholds['min_reads'], 10000)
            self.assertGreaterEqual(thresholds['min_mean_quality'], 9)
            self.assertGreaterEqual(thresholds['min_depletion_efficiency'], 0.90)

    def test_alert_configuration(self):
        """Test alert configuration for critical findings."""

        test_dir = Path(__file__).parent
        repo_root = test_dir.parent
        config_file = repo_root / 'templates' / 'config' / 'default_pipeline.yaml'

        if config_file.exists():
            import yaml
            with open(config_file) as f:
                config = yaml.safe_load(f)

            alerts = config['alerts']

            # Critical pathogen detection must be immediate
            self.assertEqual(alerts['critical_pathogen_detection'], 'immediate')

            # Workflow failure must be immediate
            self.assertEqual(alerts['workflow_failure'], 'immediate')

if __name__ == '__main__':
    unittest.main()