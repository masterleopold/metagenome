# Code Patterns & Conventions

Essential code patterns and conventions for the MinION Pathogen Screening Pipeline.

## File Handling Patterns

### File Validation
Always validate files before processing:

```python
from pathlib import Path

def process_file(file_path: str):
    file = Path(file_path)

    # Check existence
    if not file.exists():
        raise FileNotFoundError(f"File not found: {file}")

    # Check readability
    if not file.is_file():
        raise ValueError(f"Not a file: {file}")

    # Check size
    if file.stat().st_size == 0:
        raise ValueError(f"Empty file: {file}")
```

### BAM File Handling
```python
import pysam

def read_bam(bam_file: Path) -> pysam.AlignmentFile:
    """Safely open and validate BAM file."""
    if not bam_file.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    try:
        bam = pysam.AlignmentFile(str(bam_file), "rb")
        # Validate index exists
        if not Path(f"{bam_file}.bai").exists():
            pysam.index(str(bam_file))
        return bam
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file {bam_file}: {e}")
```

### FASTQ Processing
```python
from Bio import SeqIO

def process_fastq(fastq_path: Path):
    """Process FASTQ with error handling."""
    try:
        with open(fastq_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # Process record
                yield record
    except Exception as e:
        raise RuntimeError(f"Error processing FASTQ {fastq_path}: {e}")
```

## AWS Patterns

### S3 Operations
```python
import boto3
from botocore.exceptions import ClientError

s3 = boto3.client('s3', region_name='ap-northeast-1')

def safe_s3_download(bucket: str, key: str, local_path: str):
    """Download from S3 with retry logic."""
    max_retries = 3
    for attempt in range(max_retries):
        try:
            s3.download_file(bucket, key, local_path)
            return
        except ClientError as e:
            if attempt == max_retries - 1:
                raise
            time.sleep(2 ** attempt)  # Exponential backoff
```

### Lambda Function Pattern
```python
import json
import logging
from typing import Dict, Any

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def lambda_handler(event: Dict[str, Any], context: Any) -> Dict[str, Any]:
    """Standard Lambda handler pattern."""
    try:
        # Input validation
        if 'Records' not in event:
            raise ValueError("No Records in event")

        # Process event
        result = process_event(event)

        return {
            'statusCode': 200,
            'body': json.dumps(result)
        }
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        return {
            'statusCode': 500,
            'body': json.dumps({'error': str(e)})
        }
```

### EC2 Script Pattern
```python
#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path

# Add shared libraries
sys.path.append('/opt/minion/lib')
from workflow_manager import WorkflowManager

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, type=str)
    parser.add_argument('--output', required=True, type=str)
    args = parser.parse_args()

    # Validate inputs
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input not found: {input_path}")

    # Process
    workflow = WorkflowManager()
    workflow.process(input_path, args.output)

if __name__ == "__main__":
    main()
```

## Testing Patterns

### Test File Organization
```python
from pathlib import Path
import pytest

# Standard test setup
@pytest.fixture
def test_data_dir():
    """Get test data directory."""
    test_dir = Path(__file__).parent
    return test_dir / 'data'

@pytest.fixture
def config_file():
    """Get config file path."""
    test_dir = Path(__file__).parent
    repo_root = test_dir.parent
    return repo_root / 'templates' / 'config' / 'pmda_pathogens.json'
```

### Mocking AWS Services
```python
import boto3
from moto import mock_s3, mock_ec2
import pytest

@mock_s3
def test_s3_operations():
    """Test S3 operations with moto."""
    # Create mock S3
    s3 = boto3.client('s3', region_name='ap-northeast-1')
    s3.create_bucket(
        Bucket='test-bucket',
        CreateBucketConfiguration={'LocationConstraint': 'ap-northeast-1'}
    )

    # Test operations
    s3.put_object(Bucket='test-bucket', Key='test.txt', Body=b'test')
    response = s3.get_object(Bucket='test-bucket', Key='test.txt')
    assert response['Body'].read() == b'test'
```

### Integration Test Pattern
```python
@pytest.mark.integration
def test_full_pipeline():
    """Integration test - skipped in unit tests."""
    # Full pipeline test
    pass

# Run with: pytest tests/ -k "not integration"  # Skip integration tests
# Or: pytest tests/ -m integration  # Only integration tests
```

## PERV Detection Patterns

### PERV Markers
```python
# Standard PERV marker definition
PERV_MARKERS = {
    'PERV-A': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC']
    },
    'PERV-B': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAACCACCGTAGC', 'TGGAAACCTGGAAAACC']
    },
    'PERV-C': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCATAGG', 'TGGAGACCTGGAAGAAC']
    }
}
```

### Alert Pattern
```python
import boto3

sns = boto3.client('sns', region_name='ap-northeast-1')

def send_perv_alert(perv_type: str, abundance: float, run_id: str):
    """Send immediate alert for PERV detection."""
    message = {
        'alert_type': 'CRITICAL_PERV_DETECTION',
        'perv_type': perv_type,
        'abundance_copies_ml': abundance,
        'run_id': run_id,
        'timestamp': datetime.utcnow().isoformat()
    }

    sns.publish(
        TopicArn=os.environ['SNS_ALERT_TOPIC'],
        Subject=f'⚠️ CRITICAL: {perv_type} Detected',
        Message=json.dumps(message)
    )
```

## Configuration Patterns

### Loading Config
```python
import json
from pathlib import Path

def load_config(config_name: str = 'pmda_pathogens.json'):
    """Load configuration file."""
    config_path = Path(__file__).parent.parent / 'templates' / 'config' / config_name

    with open(config_path, 'r') as f:
        config = json.load(f)

    # Validate required fields
    required = ['version', 'total_pathogens', 'pathogens']
    for field in required:
        if field not in config:
            raise ValueError(f"Missing required field: {field}")

    return config
```

### Environment Variables
```python
import os

# Required environment variables for Lambda
REQUIRED_ENV = [
    'STATE_MACHINE_ARN',
    'CLUSTER_ARN',
    'SECRET_ARN',
    'DATABASE',
    'SNS_TOPIC_ARN'
]

def validate_environment():
    """Validate all required environment variables are set."""
    missing = []
    for var in REQUIRED_ENV:
        if var not in os.environ:
            missing.append(var)

    if missing:
        raise EnvironmentError(f"Missing environment variables: {missing}")
```

## Logging Patterns

### Standard Logging Setup
```python
import logging
import sys

def setup_logging(level=logging.INFO):
    """Setup standard logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('/var/log/minion-pipeline.log')
        ]
    )

    # Suppress noisy libraries
    logging.getLogger('boto3').setLevel(logging.WARNING)
    logging.getLogger('botocore').setLevel(logging.WARNING)

    return logging.getLogger(__name__)
```

## Error Handling Patterns

### Retry Decorator
```python
from functools import wraps
import time

def retry(max_attempts=3, delay=1, backoff=2):
    """Retry decorator with exponential backoff."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            attempt = 0
            while attempt < max_attempts:
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    attempt += 1
                    if attempt == max_attempts:
                        raise
                    time.sleep(delay * (backoff ** (attempt - 1)))
            return None
        return wrapper
    return decorator

# Usage
@retry(max_attempts=3, delay=2)
def unstable_operation():
    # Operation that might fail
    pass
```

## Database Patterns

### RDS Aurora Connection
```python
import boto3
import json

rds = boto3.client('rds-data', region_name='ap-northeast-1')

def execute_query(sql: str, parameters: list = None):
    """Execute query on Aurora Serverless."""
    response = rds.execute_statement(
        resourceArn=os.environ['CLUSTER_ARN'],
        secretArn=os.environ['SECRET_ARN'],
        database=os.environ['DATABASE'],
        sql=sql,
        parameters=parameters or []
    )
    return response['records']
```