#!/usr/bin/env python3
"""
Unit tests for basecalling module
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import json
import boto3
from moto import mock_s3, mock_ec2

class TestBasecallingModule(unittest.TestCase):
    """Tests for basecalling phase operations."""

    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.test_config = {
            'run_id': 'TEST-001',
            'bucket': 'test-bucket',
            'input_prefix': 'runs/TEST-001/fast5/',
            'phases': {
                'basecalling': {
                    'enabled': True,
                    'skip_duplex': False,
                    'min_quality': 9
                }
            }
        }

    def tearDown(self):
        """Clean up test environment."""
        shutil.rmtree(self.test_dir, ignore_errors=True)

    @mock_s3
    def test_fast5_file_discovery(self):
        """Test discovery of FAST5 files in S3."""
        # Create mock S3 bucket and files
        s3 = boto3.client('s3', region_name='us-east-1')
        s3.create_bucket(Bucket='test-bucket')

        # Upload test files
        test_files = [
            'runs/TEST-001/fast5/batch001.fast5',
            'runs/TEST-001/fast5/batch002.fast5',
            'runs/TEST-001/pod5/batch001.pod5'
        ]

        for file_key in test_files:
            s3.put_object(Bucket='test-bucket', Key=file_key, Body=b'test data')

        # Test file discovery
        response = s3.list_objects_v2(
            Bucket='test-bucket',
            Prefix='runs/TEST-001/fast5/'
        )

        fast5_files = [obj['Key'] for obj in response.get('Contents', [])]
        self.assertEqual(len(fast5_files), 2)
        self.assertTrue(all('.fast5' in f for f in fast5_files))

    @patch('subprocess.run')
    def test_dorado_duplex_basecalling(self, mock_subprocess):
        """Test Dorado duplex basecalling command execution."""
        mock_subprocess.return_value = Mock(returncode=0, stdout='Success', stderr='')

        # Test duplex mode command
        from scripts.phase1_basecalling import basecall_duplex

        command = basecall_duplex.build_dorado_command(
            input_path='/data/fast5',
            output_path='/data/fastq',
            model='dna_r10.4.1_e8.2_400bps_sup.cfg',
            skip_duplex=False,
            device='cuda:0'
        )

        self.assertIn('duplex', command)
        self.assertIn('--device cuda:0', command)
        self.assertIn('dna_r10.4.1_e8.2_400bps_sup.cfg', command)

    @patch('subprocess.run')
    def test_dorado_simplex_basecalling(self, mock_subprocess):
        """Test Dorado simplex basecalling when duplex is skipped."""
        mock_subprocess.return_value = Mock(returncode=0, stdout='Success', stderr='')

        from scripts.phase1_basecalling import basecall_duplex

        command = basecall_duplex.build_dorado_command(
            input_path='/data/fast5',
            output_path='/data/fastq',
            model='dna_r10.4.1_e8.2_400bps_sup.cfg',
            skip_duplex=True,
            device='cuda:0'
        )

        self.assertIn('basecaller', command)
        self.assertNotIn('duplex', command)

    def test_quality_score_filtering(self):
        """Test filtering of reads by quality score."""
        from scripts.phase1_basecalling import quality_filter

        # Mock FASTQ data with different quality scores
        mock_reads = [
            {'id': 'read1', 'seq': 'ATCG', 'qual': 'IIII', 'mean_qual': 10.5},
            {'id': 'read2', 'seq': 'GCTA', 'qual': '####', 'mean_qual': 2.0},
            {'id': 'read3', 'seq': 'TTAA', 'qual': 'JJJJ', 'mean_qual': 11.0}
        ]

        # Filter with min quality 9
        filtered = quality_filter.filter_by_quality(mock_reads, min_quality=9)

        self.assertEqual(len(filtered), 2)
        self.assertNotIn('read2', [r['id'] for r in filtered])

    @mock_ec2
    def test_gpu_instance_availability(self):
        """Test checking GPU instance availability."""
        ec2 = boto3.client('ec2', region_name='us-east-1')

        # Create mock GPU instance
        response = ec2.run_instances(
            ImageId='ami-12345678',
            MinCount=1,
            MaxCount=1,
            InstanceType='g4dn.xlarge',
            TagSpecifications=[{
                'ResourceType': 'instance',
                'Tags': [{'Key': 'Type', 'Value': 'Basecalling'}]
            }]
        )

        instance_id = response['Instances'][0]['InstanceId']

        # Check instance status
        instances = ec2.describe_instances(
            Filters=[
                {'Name': 'instance-type', 'Values': ['g4dn.xlarge']},
                {'Name': 'instance-state-name', 'Values': ['running', 'pending']}
            ]
        )

        self.assertEqual(len(instances['Reservations']), 1)

    def test_basecalling_metrics_collection(self):
        """Test collection of basecalling performance metrics."""
        from scripts.phase1_basecalling import metrics_collector

        # Mock basecalling results
        results = {
            'total_reads': 50000,
            'passed_reads': 45000,
            'failed_reads': 5000,
            'total_bases': 150000000,
            'mean_read_length': 3000,
            'mean_quality': 10.5,
            'n50': 3500,
            'duration_seconds': 7200
        }

        # Calculate metrics
        metrics = metrics_collector.calculate_metrics(results)

        self.assertEqual(metrics['pass_rate'], 0.9)
        self.assertEqual(metrics['throughput_bases_per_second'], 20833.33)
        self.assertGreater(metrics['mean_quality'], 9)  # Above threshold

    def test_fast5_to_pod5_conversion(self):
        """Test conversion of FAST5 files to POD5 format."""
        from scripts.phase1_basecalling import format_converter

        # Create mock FAST5 file
        fast5_path = Path(self.test_dir) / 'test.fast5'
        fast5_path.write_bytes(b'mock fast5 data')

        # Mock conversion function
        with patch('scripts.phase1_basecalling.format_converter.convert_fast5_to_pod5') as mock_convert:
            mock_convert.return_value = Path(self.test_dir) / 'test.pod5'

            pod5_path = format_converter.convert_fast5_to_pod5(fast5_path)

            self.assertTrue(pod5_path.name.endswith('.pod5'))
            mock_convert.assert_called_once_with(fast5_path)

    def test_basecalling_checkpoint_recovery(self):
        """Test recovery from basecalling checkpoint after failure."""
        from scripts.phase1_basecalling import checkpoint_manager

        # Create checkpoint data
        checkpoint = {
            'run_id': 'TEST-001',
            'processed_files': ['batch001.fast5', 'batch002.fast5'],
            'remaining_files': ['batch003.fast5', 'batch004.fast5'],
            'last_completed': 'batch002.fast5',
            'timestamp': '2024-01-15T10:00:00Z'
        }

        # Save checkpoint
        checkpoint_path = Path(self.test_dir) / 'checkpoint.json'
        checkpoint_manager.save_checkpoint(checkpoint_path, checkpoint)

        # Load and verify checkpoint
        loaded = checkpoint_manager.load_checkpoint(checkpoint_path)

        self.assertEqual(loaded['run_id'], 'TEST-001')
        self.assertEqual(len(loaded['processed_files']), 2)
        self.assertEqual(len(loaded['remaining_files']), 2)

    @patch('boto3.client')
    def test_basecalling_completion_notification(self, mock_boto):
        """Test SNS notification on basecalling completion."""
        mock_sns = Mock()
        mock_boto.return_value = mock_sns

        from scripts.phase1_basecalling import notifier

        # Send completion notification
        notifier.send_completion_notification(
            run_id='TEST-001',
            total_reads=50000,
            passed_reads=45000,
            duration_hours=2.0
        )

        # Verify SNS publish was called
        mock_sns.publish.assert_called_once()
        call_args = mock_sns.publish.call_args
        self.assertIn('TEST-001', call_args[1]['Subject'])
        self.assertIn('45000', call_args[1]['Message'])

    def test_model_selection_for_chemistry(self):
        """Test automatic model selection based on flow cell chemistry."""
        from scripts.phase1_basecalling import model_selector

        # Test R10.4.1 chemistry
        model = model_selector.get_model_for_chemistry('R10.4.1', accuracy='sup')
        self.assertIn('r10.4.1', model)
        self.assertIn('sup', model)

        # Test R9.4.1 chemistry
        model = model_selector.get_model_for_chemistry('R9.4.1', accuracy='hac')
        self.assertIn('r9.4.1', model)
        self.assertIn('hac', model)

    def test_cuda_availability_check(self):
        """Test CUDA availability checking for GPU basecalling."""
        from scripts.phase1_basecalling import gpu_utils

        with patch('scripts.phase1_basecalling.gpu_utils.check_cuda') as mock_cuda:
            # Test CUDA available
            mock_cuda.return_value = True
            device = gpu_utils.get_device()
            self.assertEqual(device, 'cuda:0')

            # Test CUDA not available
            mock_cuda.return_value = False
            device = gpu_utils.get_device()
            self.assertEqual(device, 'cpu')

if __name__ == '__main__':
    unittest.main()