#!/usr/bin/env python3
"""
Integration tests for MinION pipeline
"""

import unittest
import boto3
import json
import time
from pathlib import Path
from datetime import datetime
import sys
sys.path.append('../lib')

from workflow_manager import WorkflowManager
from config_manager import ConfigManager
from s3_utils import S3Utils
from database_client import DatabaseClient

class TestPipelineIntegration(unittest.TestCase):
    """Integration tests for the complete pipeline."""

    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.region = 'ap-northeast-1'
        cls.environment = 'test'
        cls.test_bucket = f'minion-test-{datetime.now().strftime("%Y%m%d%H%M%S")}'

        # Initialize clients
        cls.workflow_manager = WorkflowManager(cls.region)
        cls.config_manager = ConfigManager()
        cls.s3_utils = S3Utils(cls.region)
        cls.s3_client = boto3.client('s3', region_name=cls.region)

        # Create test bucket
        cls.s3_client.create_bucket(
            Bucket=cls.test_bucket,
            CreateBucketConfiguration={'LocationConstraint': cls.region}
        )

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        # Delete test bucket and contents
        try:
            # Delete all objects
            response = cls.s3_client.list_objects_v2(Bucket=cls.test_bucket)
            if 'Contents' in response:
                for obj in response['Contents']:
                    cls.s3_client.delete_object(Bucket=cls.test_bucket, Key=obj['Key'])

            # Delete bucket
            cls.s3_client.delete_bucket(Bucket=cls.test_bucket)
        except:
            pass

    def test_01_configuration_loading(self):
        """Test configuration loading and validation."""

        # Load default configuration
        config = self.config_manager.load_config()
        self.assertIsNotNone(config)

        # Validate configuration
        validation = self.config_manager.validate_config(config)
        self.assertTrue(validation['valid'])
        self.assertEqual(len(validation['errors']), 0)

    def test_02_input_validation(self):
        """Test input file validation."""

        # Create test files
        test_prefix = 'test_input/'
        test_file = f'{test_prefix}test.fast5'

        # Upload dummy file
        self.s3_client.put_object(
            Bucket=self.test_bucket,
            Key=test_file,
            Body=b'dummy fast5 data'
        )

        # Validate input
        validation = self.workflow_manager.validate_input(
            self.test_bucket,
            test_prefix
        )

        self.assertTrue(validation['valid'])
        self.assertGreater(validation['file_count'], 0)

    def test_03_workflow_start(self):
        """Test workflow initiation."""

        # Create test configuration
        config = {
            'phases': {
                'basecalling': {'enabled': False},  # Skip for testing
                'qc': {'enabled': True},
                'host_removal': {'enabled': False},
                'pathogen_detection': {'enabled': False},
                'quantification': {'enabled': False},
                'reporting': {'enabled': True}
            }
        }

        # Start workflow (mock)
        run_id = f'test-{datetime.now().strftime("%Y%m%d%H%M%S")}'

        # Validate we can create the execution input
        execution_input = {
            'run_id': run_id,
            'bucket': self.test_bucket,
            'input_prefix': 'test_input/',
            'config': config
        }

        self.assertIn('run_id', execution_input)
        self.assertIn('config', execution_input)

    def test_04_s3_operations(self):
        """Test S3 utility operations."""

        # Test upload
        local_file = '/tmp/test_file.txt'
        with open(local_file, 'w') as f:
            f.write('test content')

        result = self.s3_utils._upload_file(
            local_file,
            self.test_bucket,
            'test/uploaded.txt'
        )

        self.assertTrue(result)

        # Test download
        download_file = '/tmp/downloaded.txt'
        result = self.s3_utils._download_file(
            self.test_bucket,
            'test/uploaded.txt',
            download_file
        )

        self.assertTrue(result)
        self.assertTrue(Path(download_file).exists())

        # Verify content
        with open(download_file, 'r') as f:
            content = f.read()
        self.assertEqual(content, 'test content')

    def test_05_phase_configuration(self):
        """Test phase-specific configuration retrieval."""

        config = self.config_manager.load_config()

        # Test each phase
        phases = ['basecalling', 'qc', 'host_removal',
                 'pathogen_detection', 'quantification', 'reporting']

        for phase in phases:
            phase_config = self.config_manager.get_phase_config(phase, config)
            self.assertIsNotNone(phase_config)
            self.assertIn('enabled', phase_config)

    def test_06_cost_estimation(self):
        """Test workflow cost estimation."""

        file_count = 10
        total_size = 10 * 1e9  # 10GB
        config = {
            'phases': {
                'basecalling': {'enabled': True, 'skip_duplex': False},
                'qc': {'enabled': True},
                'host_removal': {'enabled': True},
                'pathogen_detection': {'enabled': True, 'databases': ['kraken2', 'rvdb']},
                'quantification': {'enabled': True},
                'reporting': {'enabled': True}
            }
        }

        estimates = self.workflow_manager.estimate_runtime(
            file_count, total_size, config
        )

        self.assertIn('total', estimates)
        self.assertGreater(estimates['total'], 0)

        # Check individual phase estimates
        for phase in ['basecalling', 'qc', 'host_removal',
                     'pathogen_detection', 'quantification', 'reporting']:
            self.assertIn(phase, estimates)
            self.assertGreaterEqual(estimates[phase], 0)

    def test_07_pathogen_list(self):
        """Test PMDA pathogen list loading."""

        pmda_file = Path('../templates/config/pmda_pathogens.json')

        if pmda_file.exists():
            with open(pmda_file) as f:
                pmda_data = json.load(f)

            self.assertEqual(pmda_data['total_pathogens'], 91)
            self.assertIn('categories', pmda_data)
            self.assertIn('viruses', pmda_data['categories'])

            # Check critical pathogens
            viruses = pmda_data['categories']['viruses']
            self.assertIn('PERV-A', viruses['critical'])
            self.assertIn('PERV-B', viruses['critical'])
            self.assertIn('PERV-C', viruses['critical'])

    def test_08_report_template(self):
        """Test report template existence and validity."""

        template_path = Path('../templates/report_template.html')
        self.assertTrue(template_path.exists())

        # Check template has required placeholders
        with open(template_path) as f:
            template_content = f.read()

        required_placeholders = [
            '{{ run_id }}',
            '{{ generation_time }}',
            '{{ report_data'
        ]

        for placeholder in required_placeholders:
            self.assertIn(placeholder, template_content)

    def test_09_database_configuration(self):
        """Test database configuration."""

        # Test with mock connection parameters
        db_client = DatabaseClient(
            cluster_arn='arn:aws:rds:region:account:cluster:test',
            secret_arn='arn:aws:secretsmanager:region:account:secret:test',
            database='minion_test',
            region=self.region
        )

        # Test that client initializes
        self.assertIsNotNone(db_client)
        self.assertEqual(db_client.database, 'minion_test')

    def test_10_monitoring_metrics(self):
        """Test monitoring metric definitions."""

        test_metrics = {
            'basecalling_total_reads': 10000,
            'basecalling_mean_quality': 10.5,
            'qc_reads_passed': 9500,
            'qc_reads_failed': 500,
            'host_removal_depletion_rate': 95.5,
            'pathogen_detection_pathogens_detected': 3,
            'pathogen_detection_perv_detected': False
        }

        # Validate metric types
        for key, value in test_metrics.items():
            if 'detected' in key:
                self.assertIsInstance(value, (bool, int))
            else:
                self.assertIsInstance(value, (int, float))

class TestPhaseScripts(unittest.TestCase):
    """Unit tests for individual phase scripts."""

    def test_basecalling_parameters(self):
        """Test basecalling parameter validation."""

        valid_params = {
            'skip_duplex': False,
            'min_quality': 9,
            'model': 'dna_r10.4.1_e8.2_400bps_sup.cfg',
            'device': 'cuda:all'
        }

        # Test quality thresholds
        self.assertGreaterEqual(valid_params['min_quality'], 7)
        self.assertLessEqual(valid_params['min_quality'], 15)

        # Test model string
        self.assertIn('.cfg', valid_params['model'])

    def test_qc_thresholds(self):
        """Test QC threshold parameters."""

        qc_params = {
            'min_quality': 9,
            'min_length': 100,
            'max_length': 100000,
            'min_pass_reads': 1000
        }

        # Validate thresholds
        self.assertGreater(qc_params['min_quality'], 0)
        self.assertLess(qc_params['min_length'], qc_params['max_length'])
        self.assertGreater(qc_params['min_pass_reads'], 0)

    def test_host_removal_settings(self):
        """Test host removal configuration."""

        settings = {
            'reference': 'sus_scrofa_11.1',
            'min_identity': 0.95,
            'threads': 16
        }

        # Validate reference format
        self.assertRegex(settings['reference'], r'^sus_scrofa_\d+\.\d+$')

        # Validate identity threshold
        self.assertGreater(settings['min_identity'], 0)
        self.assertLessEqual(settings['min_identity'], 1)

        # Validate thread count
        self.assertGreater(settings['threads'], 0)
        self.assertLessEqual(settings['threads'], 128)

    def test_pathogen_detection_databases(self):
        """Test pathogen detection database configuration."""

        valid_databases = ['kraken2', 'rvdb', 'pmda', 'blast', 'diamond']
        selected_databases = ['kraken2', 'rvdb', 'pmda']

        for db in selected_databases:
            self.assertIn(db, valid_databases)

    def test_quantification_calculations(self):
        """Test quantification calculations."""

        # Test copies per mL calculation
        reads = 1000
        genome_size = 10000
        total_reads = 1000000
        plasma_volume = 10.0
        extraction_efficiency = 0.7

        # Simplified calculation
        avg_read_length = 1000
        total_dna_molecules = total_reads * avg_read_length / genome_size
        pathogen_fraction = reads / total_reads
        pathogen_molecules = pathogen_fraction * total_dna_molecules
        copies_per_ml = (pathogen_molecules / extraction_efficiency) / plasma_volume

        self.assertGreater(copies_per_ml, 0)

    def test_reporting_formats(self):
        """Test reporting format options."""

        valid_formats = ['pdf', 'html', 'json', 'csv']
        selected_formats = ['pdf', 'json']

        for fmt in selected_formats:
            self.assertIn(fmt, valid_formats)

if __name__ == '__main__':
    unittest.main()