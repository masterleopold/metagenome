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

---

## v2.0 Patterns (Type Safety & Repository)

> **New in v2.0**: See [NEW_PATTERNS_GUIDE.md](NEW_PATTERNS_GUIDE.md) for complete guide

### Type-Safe Models (Pydantic)

#### PERV Detection with Validation
```python
from lib.models.pathogen import (
    PERVTypingOutput,
    PERVDetectionResult,
    PERVSubtype,
    PathogenConfidence
)
from pathlib import Path

# ❌ OLD WAY (v1.0) - error-prone dicts
def identify_perv_old(bam_file):
    return {
        'PERV-A': {'reads': 42, 'coverage': 0.85},
        'PERV-B': {'reads': 0, 'coverage': 0.0}
    }

# ✅ NEW WAY (v2.0) - type-safe with auto-validation
def identify_perv_v2(bam_file: Path) -> PERVTypingOutput:
    """Type-safe PERV detection with automatic confidence calculation."""
    return PERVTypingOutput(
        run_id=os.environ['RUN_ID'],
        bam_file=bam_file,
        detections={
            PERVSubtype.PERV_A: PERVDetectionResult(
                subtype=PERVSubtype.PERV_A,
                reads_aligned=42,
                coverage=0.85,  # ✅ Auto-validated: 0.0 ≤ coverage ≤ 1.0
                specific_motifs_found=['ATGGCAGCCACCACAGC'],
                mean_identity=96.5,  # ✅ Auto-validated: 0.0 ≤ identity ≤ 100.0
                confidence=PathogenConfidence.HIGH  # ✅ Auto-calculated from metrics
            )
        }
    )

# Type-safe property access
result = identify_perv_v2(Path("/data/aligned.bam"))
if result.requires_sns_alert:  # ✅ Type-safe boolean property
    send_alert(result.to_audit_log())  # ✅ Type-safe method
```

#### PMDA 91 Pathogen Screening
```python
from lib.models.pathogen import PMDA91PathogenResult, PathogenDetectionResult

def screen_91_pathogens() -> PMDA91PathogenResult:
    """Type-safe 91 pathogen screening with validation."""
    detections = [
        PathogenDetectionResult(
            pathogen_name="Hantavirus",
            method=PathogenDetectionMethod.KRAKEN2,
            reads_classified=150,
            confidence=PathogenConfidence.HIGH,
            database_version="PMDA_2024_v3"
        )
    ]

    result = PMDA91PathogenResult(
        run_id="RUN-2025-001",
        total_pathogens_screened=91,  # ✅ Auto-validated: must equal 91
        detections=detections
    )

    # ✅ Auto-validated counts
    assert result.positive_count == 1
    assert result.negative_count == 90

    return result
```

#### Workflow Execution with QC
```python
from lib.models.workflow import WorkflowExecution, WorkflowStatus, QCMetrics

def create_workflow() -> WorkflowExecution:
    """Type-safe workflow with automatic QC validation."""
    qc = QCMetrics(
        total_reads=5_000_000,
        total_bases_gb=15.5,
        mean_qscore=12.5,  # ✅ Auto-validated: 0.0 ≤ qscore ≤ 50.0
        median_qscore=13.0,
        n50=8500,
        median_read_length=7200,
        active_channels=450  # ✅ Auto-validated: 0 ≤ channels ≤ 512
    )

    # ✅ Automatic pass/fail logic
    if not qc.passes_qc:
        raise ValueError(f"QC failed: {qc}")

    return WorkflowExecution(
        run_id="RUN-2025-001",
        status=WorkflowStatus.BASECALLING,
        qc_metrics=qc
    )
```

### Repository Pattern

#### Production: RDS Repository
```python
from lib.repositories.rds_repository import RDSWorkflowRepository, RDSPathogenDetectionRepository
from lib.models.workflow import WorkflowExecution, WorkflowStatus
import os

# Initialize repository
workflow_repo = RDSWorkflowRepository(
    cluster_arn=os.environ['RDS_CLUSTER_ARN'],
    secret_arn=os.environ['RDS_SECRET_ARN'],
    database='minion_metadata',
    region='ap-northeast-1'
)

# ❌ OLD WAY (v1.0) - raw SQL everywhere
def save_workflow_old(run_id, status):
    rds = boto3.client('rds-data')
    rds.execute_statement(
        resourceArn=os.environ['CLUSTER_ARN'],
        secretArn=os.environ['SECRET_ARN'],
        sql=f"INSERT INTO workflows (run_id, status) VALUES ('{run_id}', '{status}')"
    )

# ✅ NEW WAY (v2.0) - type-safe repository
def save_workflow_v2(workflow: WorkflowExecution):
    workflow_repo.create(workflow)  # ✅ Type-safe, SQL injection safe

# Query with type safety
workflow = workflow_repo.get("RUN-2025-001")
if workflow and workflow.status == WorkflowStatus.COMPLETED:
    print(f"Workflow completed: {workflow.run_id}")

# Update status
workflow_repo.update_status("RUN-2025-001", WorkflowStatus.PATHOGEN_DETECTION)
```

#### Testing: SQLite Repository (No AWS!)
```python
from lib.repositories.sqlite_repository import SQLiteWorkflowRepository
from lib.models.workflow import WorkflowExecution, WorkflowStatus

# ❌ OLD WAY (v1.0) - complex boto3 mocking
@mock_rds_data
def test_workflow_old():
    # Setup complex mocks...
    pass

# ✅ NEW WAY (v2.0) - in-memory SQLite, 10x faster
def test_workflow_v2():
    """Test with in-memory SQLite - no AWS credentials needed!"""
    # Same interface as RDS repository!
    test_repo = SQLiteWorkflowRepository(db_path=":memory:")

    workflow = WorkflowExecution(
        run_id="TEST-001",
        status=WorkflowStatus.BASECALLING
    )

    # Same methods as production
    test_repo.create(workflow)
    retrieved = test_repo.get("TEST-001")

    assert retrieved.run_id == "TEST-001"
    assert retrieved.status == WorkflowStatus.BASECALLING
```

#### Pathogen Detection Repository
```python
from lib.repositories.rds_repository import RDSPathogenDetectionRepository
from lib.models.pathogen import PERVTypingOutput

pathogen_repo = RDSPathogenDetectionRepository(
    cluster_arn=os.environ['RDS_CLUSTER_ARN'],
    secret_arn=os.environ['RDS_SECRET_ARN']
)

# Save type-safe PERV detection
perv_result: PERVTypingOutput = identify_perv_v2(bam_file)
pathogen_repo.create_perv_detection(perv_result)

# Query PERV detections
recent_perv = pathogen_repo.get_recent_perv_detections(days=30)
for detection in recent_perv:
    print(f"{detection.run_id}: {detection.positive_subtypes}")
```

### Unified Logging (AWS Lambda Powertools)

#### Lambda Handler with Full Observability
```python
from aws_lambda_powertools import Logger, Tracer, Metrics
from aws_lambda_powertools.metrics import MetricUnit
from aws_lambda_powertools.utilities.typing import LambdaContext
from lib.logging.logger import get_logger, AuditLogger

# Initialize observability tools
logger = get_logger("pathogen-detection")
tracer = Tracer(service="pathogen-detection")
metrics = Metrics(namespace="MinION/Pipeline", service="pathogen-detection")
audit = AuditLogger(service="pathogen-detection")

# ❌ OLD WAY (v1.0) - basic print statements
def lambda_handler_old(event, context):
    print(f"Starting run {event['run_id']}")
    # No correlation IDs, no structured logs, no metrics

# ✅ NEW WAY (v2.0) - full observability
@logger.inject_lambda_context(log_event=True)
@tracer.capture_lambda_handler
@metrics.log_metrics(capture_cold_start_metric=True)
def lambda_handler(event: dict, context: LambdaContext) -> dict:
    """Lambda handler with structured logging and metrics."""
    run_id = event['run_id']

    # ✅ Structured JSON logging with correlation ID
    logger.info(
        "Phase 4 pathogen detection started",
        extra={
            "run_id": run_id,
            "workflow_id": event.get('workflow_id'),
            "instance_type": "c5.4xlarge",
            "databases": ["kraken2", "blast"]
        }
    )

    try:
        # Process detection
        perv_result = detect_perv(run_id)

        # ✅ PMDA audit logging (CRITICAL)
        if perv_result.requires_sns_alert:
            audit.log_perv_detection(
                run_id=run_id,
                subtypes_detected=[s.value for s in perv_result.positive_subtypes],
                confidence_levels={s.value: "HIGH" for s in perv_result.positive_subtypes},
                operator_email=event.get('operator_email', 'unknown')
            )

        # ✅ CloudWatch Metrics
        metrics.add_metric(name="PERVDetected", unit=MetricUnit.Count, value=1)
        metrics.add_metric(name="ProcessingTime", unit=MetricUnit.Seconds, value=120)

        # ✅ X-Ray tracing
        tracer.put_annotation(key="RunID", value=run_id)
        tracer.put_metadata(key="result", value=perv_result.model_dump())

        logger.info("Detection completed successfully", extra={"run_id": run_id})

        return {
            'statusCode': 200,
            'body': perv_result.model_dump_json()
        }

    except Exception as e:
        logger.error(
            "Detection failed",
            extra={"run_id": run_id, "error": str(e)},
            exc_info=True
        )
        metrics.add_metric(name="DetectionErrors", unit=MetricUnit.Count, value=1)
        raise
```

#### PMDA Audit Logging
```python
from lib.logging.logger import AuditLogger

audit = AuditLogger(service="pathogen-detection")

# Log PERV detection (CRITICAL for PMDA)
audit.log_perv_detection(
    run_id="RUN-2025-001",
    subtypes_detected=["PERV-A"],
    confidence_levels={"PERV-A": "HIGH"},
    operator_email="operator@example.com",
    additional_context={
        "reads_aligned": 42,
        "coverage": 0.95,
        "mean_identity": 98.5
    }
)

# Log pathogen screening
audit.log_pathogen_screening(
    run_id="RUN-2025-001",
    total_screened=91,
    positive_count=2,
    pathogens_detected=["Hantavirus", "Polyomavirus"],
    operator_email="operator@example.com"
)

# Log workflow completion
audit.log_workflow_completion(
    run_id="RUN-2025-001",
    final_status="COMPLETED",
    duration_seconds=7200,
    qc_passed=True,
    operator_email="operator@example.com"
)
```

### CloudWatch Audit Queries

#### Execute Pre-Built Queries
```python
from lib.audit.cloudwatch_queries import execute_query, get_query_by_name
import time

# Get all PERV detections in last 30 days
end_time = int(time.time())
start_time = end_time - (30 * 24 * 60 * 60)

results = execute_query(
    log_group_name="/aws/lambda/pathogen-detection",
    query_name="perv_detections",
    start_time=start_time,
    end_time=end_time
)

# Process results
for result in results['results']:
    fields = {r['field']: r['value'] for r in result}
    print(f"PERV detection: {fields['run_id']} - {fields['subtypes_detected']}")
```

#### Available Audit Queries
```python
# 1. PERV detections (CRITICAL for PMDA)
execute_query(log_group_name, "perv_detections", start_time, end_time)

# 2. Failed workflows by phase
execute_query(log_group_name, "failed_workflows", start_time, end_time)

# 3. Workflow duration metrics
execute_query(log_group_name, "workflow_durations", start_time, end_time)

# 4. 91 pathogen screening results
execute_query(log_group_name, "pathogen_screening", start_time, end_time)

# 5. Complete audit trail for specific run
execute_query(
    log_group_name,
    "run_audit_trail",
    start_time,
    end_time,
    run_id="RUN-2025-001"
)

# 6. QC failures
execute_query(log_group_name, "qc_failures", start_time, end_time)

# 7. Operator activity tracking
execute_query(log_group_name, "operator_activity", start_time, end_time)

# 8. High-priority pathogen alerts
execute_query(log_group_name, "high_priority_pathogens", start_time, end_time)
```

#### Query PMDA Compliance Report
```python
from lib.audit.cloudwatch_queries import execute_query
import time

def generate_pmda_compliance_report(days: int = 30):
    """Generate PMDA compliance report for last N days."""
    end_time = int(time.time())
    start_time = end_time - (days * 24 * 60 * 60)
    log_group = "/aws/lambda/pathogen-detection"

    # Get all workflow completions
    workflows = execute_query(
        log_group_name=log_group,
        query_name="recent_workflows",
        start_time=start_time,
        end_time=end_time
    )

    # Get PERV detections
    perv_alerts = execute_query(
        log_group_name=log_group,
        query_name="perv_detections",
        start_time=start_time,
        end_time=end_time
    )

    # Get QC failures
    qc_failures = execute_query(
        log_group_name=log_group,
        query_name="qc_failures",
        start_time=start_time,
        end_time=end_time
    )

    return {
        "total_workflows": len(workflows['results']),
        "perv_detections": len(perv_alerts['results']),
        "qc_failures": len(qc_failures['results']),
        "compliance_period_days": days
    }
```

### Migration Pattern (v1.0 → v2.0)

#### Gradual Migration Strategy
```python
# STEP 1: Update function signature to return Pydantic model
# Keep old dict logic for now
def identify_perv_type(bam_file: Path) -> PERVTypingOutput:
    """Migrated to v2.0 - returns type-safe model."""
    # Old logic (keep temporarily)
    old_result = _identify_perv_type_legacy(bam_file)

    # Convert to new model
    detections = {}
    for subtype_str, data in old_result.items():
        subtype = PERVSubtype(subtype_str)
        detections[subtype] = PERVDetectionResult(
            subtype=subtype,
            reads_aligned=data['reads'],
            coverage=data['coverage'],
            specific_motifs_found=data.get('motifs', []),
            mean_identity=data.get('identity', 0.0),
            confidence=PathogenConfidence.HIGH  # Will auto-calculate
        )

    return PERVTypingOutput(
        run_id=os.environ['RUN_ID'],
        bam_file=bam_file,
        detections=detections
    )

# STEP 2: Update callers to use new model
result = identify_perv_type(bam_file)
if result.requires_sns_alert:  # ✅ Use type-safe property
    send_alert(result.to_audit_log())

# STEP 3: Remove legacy code after validation
```

---

**See also:**
- [NEW_PATTERNS_GUIDE.md](NEW_PATTERNS_GUIDE.md) - Complete v2.0 guide
- [API_REFERENCE_V2.md](API_REFERENCE_V2.md) - Full API reference
- [REFACTORING_SUMMARY.md](REFACTORING_SUMMARY.md) - Implementation summary