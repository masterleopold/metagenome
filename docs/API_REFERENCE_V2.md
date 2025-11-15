# API Reference: Type-Safe Models & Repositories (v2.0)

**Version**: 2.0
**Date**: 2025-01-15
**Status**: Production Ready

---

## Table of Contents

1. [Pydantic Models](#pydantic-models)
2. [Repository Interfaces](#repository-interfaces)
3. [Logging API](#logging-api)
4. [CloudWatch Audit Queries](#cloudwatch-audit-queries)

---

## Pydantic Models

### lib.models.pathogen

#### PERVTypingOutput

Complete PERV typing results with automatic alert logic.

**Fields**:
- `run_id: str` - Unique run identifier (pattern: `^[A-Z0-9][A-Z0-9_-]{0,63}$`)
- `bam_file: Path` - Path to aligned BAM file
- `detections: dict[PERVSubtype, PERVDetectionResult]` - Detection results by subtype
- `analysis_timestamp: datetime` - Analysis completion time (default: `datetime.utcnow()`)
- `pipeline_version: str` - Pipeline version (default: "2.1.0")

**Properties**:
- `requires_sns_alert: bool` - Returns `True` if any PERV detection requires alert
- `positive_subtypes: list[PERVSubtype]` - List of detected PERV subtypes

**Methods**:
- `to_audit_log() -> dict` - Convert to PMDA compliance audit log format

**Example**:
```python
from lib.models.pathogen import PERVTypingOutput, PERVSubtype

result = PERVTypingOutput(
    run_id="RUN-001",
    bam_file=Path("test.bam"),
    detections={
        PERVSubtype.PERV_A: PERVDetectionResult(...)
    }
)

if result.requires_sns_alert:
    sns.publish(Message=json.dumps(result.to_audit_log()))
```

---

#### PERVDetectionResult

Individual PERV subtype detection result.

**Fields**:
- `subtype: PERVSubtype` - PERV subtype (PERV-A, PERV-B, or PERV-C)
- `reads_aligned: int` - Number of aligned reads (ge=0)
- `coverage: float` - Genome coverage fraction (0.0-1.0)
- `specific_motifs_found: list[str]` - PERV-specific motifs detected
- `confidence: PathogenConfidence` - Detection confidence (auto-calculated if not provided)
- `mean_identity: float` - Mean alignment identity percentage (0.0-100.0)

**Properties**:
- `requires_alert: bool` - Returns `True` if confidence is MEDIUM or HIGH

**Auto-Validation**:
- Confidence auto-calculated based on reads + coverage + identity
- `HIGH`: ≥10 reads AND ≥80% coverage AND ≥95% identity
- `MEDIUM`: ≥3 reads OR ≥50% coverage
- `LOW`: >0 reads
- `NEGATIVE`: 0 reads

---

#### PMDA91PathogenResult

Complete 91 pathogen screening results with validation.

**Fields**:
- `run_id: str` - Unique run identifier
- `all_91_tested: bool` - Confirmation all 91 pathogens tested
- `positive_detections: list[PathogenDetectionResult]` - Positive detections
- `negative_count: int` - Number of negative results (0-91)
- `kraken2_db_version: str` - Kraken2 database version
- `blast_db_version: str` - BLAST database version
- `analysis_timestamp: datetime` - Analysis completion time

**Validation**:
- Automatically validates: `len(positive_detections) + negative_count == 91`
- Raises `ValueError` if count mismatch

**Properties**:
- `requires_alert: bool` - Returns `True` if high-priority pathogen detected

**Example**:
```python
from lib.models.pathogen import PMDA91PathogenResult

result = PMDA91PathogenResult(
    run_id="RUN-001",
    all_91_tested=True,
    positive_detections=[detection1, detection2],  # 2 pathogens
    negative_count=89,  # 91 - 2 = 89 ✅
    kraken2_db_version="pmda_2024.2",
    blast_db_version="rvdb_2024.1"
)
# Validation passes: 2 + 89 = 91
```

---

### lib.models.workflow

#### WorkflowExecution

Complete workflow execution tracking.

**Fields** (Core):
- `run_id: str` - Unique run identifier (pattern validated)
- `sample_id: str` - Sample identifier
- `donor_id: Optional[str]` - Anonymized donor ID
- `status: WorkflowStatus` - Current workflow status
- `started_at: datetime` - Workflow start time
- `completed_at: Optional[datetime]` - Workflow completion time
- `operator_name: str` - Operator who initiated run
- `operator_email: str` - Operator email

**Fields** (Metrics):
- `qc_metrics: Optional[QCMetrics]` - Quality control metrics
- `host_removal_metrics: Optional[HostRemovalMetrics]` - Host removal stats

**Properties**:
- `is_completed: bool` - Returns `True` if status is COMPLETED
- `is_failed: bool` - Returns `True` if status is FAILED
- `duration_hours: Optional[float]` - Workflow duration in hours

**Methods**:
- `to_audit_log() -> dict` - Convert to PMDA audit log format

---

#### QCMetrics

Quality control metrics with pass/fail logic.

**Fields**:
- `total_reads: int` - Total reads after basecalling (ge=0)
- `total_bases_gb: float` - Total bases in gigabases (ge=0.0)
- `mean_qscore: float` - Mean quality score (0.0-50.0)
- `median_qscore: float` - Median quality score (0.0-50.0)
- `n50: int` - N50 read length in bp (ge=0)
- `median_read_length: int` - Median read length in bp (ge=0)
- `active_channels: int` - Active sequencing channels (0-512)

**Properties**:
- `passes_qc: bool` - Returns `True` if meets thresholds:
  - Mean Q-score ≥ 9.0
  - Total bases ≥ 0.1 GB
  - Active channels ≥ 100

**Example**:
```python
from lib.models.workflow import QCMetrics

qc = QCMetrics(
    total_reads=1_000_000,
    total_bases_gb=1.5,
    mean_qscore=12.5,  # ✅ Passes (≥9.0)
    median_qscore=13.0,
    n50=5000,
    median_read_length=3000,
    active_channels=450  # ✅ Passes (≥100)
)

assert qc.passes_qc == True
```

---

## Repository Interfaces

### lib.repositories.interfaces

#### WorkflowRepository (Protocol)

Interface for workflow execution storage.

**Methods**:

##### `create(workflow: WorkflowExecution) -> str`

Create new workflow record.

**Parameters**:
- `workflow: WorkflowExecution` - Workflow model instance

**Returns**:
- `str` - Created workflow run_id

**Raises**:
- `ValueError` - If run_id already exists
- `RuntimeError` - On database error

**Example**:
```python
workflow = WorkflowExecution(run_id="RUN-001", ...)
run_id = repo.create(workflow)
```

---

##### `get(run_id: str) -> Optional[WorkflowExecution]`

Retrieve workflow by run_id.

**Parameters**:
- `run_id: str` - Unique run identifier

**Returns**:
- `WorkflowExecution` - Workflow instance if found
- `None` - If not found

---

##### `update_status(run_id: str, status: WorkflowStatus, error_message: Optional[str] = None) -> None`

Update workflow status.

**Parameters**:
- `run_id: str` - Workflow identifier
- `status: WorkflowStatus` - New status
- `error_message: Optional[str]` - Error message if status is FAILED

**Raises**:
- `ValueError` - If run_id not found

---

##### `list_by_status(status: WorkflowStatus, limit: int = 100, offset: int = 0) -> list[WorkflowExecution]`

List workflows by status.

**Parameters**:
- `status: WorkflowStatus` - Filter by status
- `limit: int` - Maximum results (default: 100)
- `offset: int` - Pagination offset (default: 0)

**Returns**:
- `list[WorkflowExecution]` - List of matching workflows

---

##### `list_recent(limit: int = 50, days: int = 7) -> list[WorkflowExecution]`

List recent workflows.

**Parameters**:
- `limit: int` - Maximum results (default: 50)
- `days: int` - Look back days (default: 7)

**Returns**:
- `list[WorkflowExecution]` - Recent workflows, newest first

---

### Implementations

#### RDSWorkflowRepository

Production implementation using RDS Data API.

**Constructor**:
```python
from lib.repositories.rds_repository import RDSWorkflowRepository

repo = RDSWorkflowRepository(
    cluster_arn="arn:aws:rds:ap-northeast-1:...",
    secret_arn="arn:aws:secretsmanager:ap-northeast-1:...",
    database="minion_metadata",  # Optional, default shown
    region="ap-northeast-1"  # Optional, default shown
)
```

---

#### SQLiteWorkflowRepository

Testing implementation using local SQLite.

**Constructor**:
```python
from lib.repositories.sqlite_repository import SQLiteWorkflowRepository

# In-memory (for tests)
repo = SQLiteWorkflowRepository(db_path=":memory:")

# Persistent file
repo = SQLiteWorkflowRepository(db_path="/tmp/test.db")
```

**Benefits**:
- No AWS credentials needed
- 10x faster than RDS mocking
- Same interface as production repository

---

## Logging API

### lib.logging.logger

#### get_logger(service: str, level: str = "INFO") -> Logger

Get or create AWS Lambda Powertools logger.

**Parameters**:
- `service: str` - Service name (e.g., "pathogen-detection")
- `level: str` - Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

**Returns**:
- `Logger` - Configured logger instance

**Example**:
```python
from lib.logging.logger import get_logger

logger = get_logger("pathogen-detection", level="INFO")

logger.info("Processing started", extra={"run_id": "RUN-001"})
logger.warning("High memory usage", extra={"memory_gb": 120})
logger.error("Database connection failed", extra={"error": str(e)})
```

---

#### AuditLogger

Specialized logger for PMDA audit trail.

**Constructor**:
```python
from lib.logging.logger import AuditLogger

audit = AuditLogger(service="pathogen-detection")
```

**Methods**:

##### `log_perv_detection(run_id, subtypes_detected, confidence_levels, operator_email, **extra)`

Log PERV detection (CRITICAL for PMDA).

**Parameters**:
- `run_id: str` - Workflow run ID
- `subtypes_detected: list[str]` - PERV subtypes detected
- `confidence_levels: dict[str, str]` - Confidence per subtype
- `operator_email: str` - Operator email
- `**extra` - Additional metadata

**Example**:
```python
audit.log_perv_detection(
    run_id="RUN-001",
    subtypes_detected=["PERV-A"],
    confidence_levels={"PERV-A": "HIGH"},
    operator_email="tanaka@example.com"
)
```

---

##### `log_pathogen_screening(run_id, all_91_tested, positive_count, positive_pathogens, operator_email, **extra)`

Log PMDA 91 pathogen screening results.

---

##### `log_qc_metrics(run_id, qc_passed, metrics, **extra)`

Log QC analysis results.

---

##### `log_error(run_id, phase, error_message, error_type, operator_email, **extra)`

Log pipeline error for audit trail.

---

## CloudWatch Audit Queries

### lib.audit.cloudwatch_queries

#### execute_query(log_group_name, query_name, start_time, end_time, **params) -> dict

Execute pre-built CloudWatch Logs Insights query.

**Parameters**:
- `log_group_name: str` - CloudWatch log group
- `query_name: str` - Query name (see list below)
- `start_time: int` - Start time (Unix timestamp)
- `end_time: int` - End time (Unix timestamp)
- `**params` - Query parameters for substitution

**Returns**:
- `dict` - Query results from CloudWatch

**Available Queries**:
1. `perv_detections` - All PERV alerts
2. `failed_workflows` - Error analysis
3. `workflow_durations` - Performance metrics
4. `pathogen_screening` - 91 pathogen results
5. `run_audit_trail` - Complete run history (requires `run_id` parameter)
6. `qc_failures` - QC threshold failures
7. `operator_activity` - User actions (requires `operator_email` parameter)
8. `high_priority_pathogens` - PERV/HantV/EEEV
9. `database_operations` - DB success rates
10. `phase_performance` - Duration by phase
11. `recent_workflows` - Last 7 days
12. `error_rate_by_phase` - SLA monitoring

**Example**:
```python
from lib.audit.cloudwatch_queries import execute_query
import time

# Get PERV detections in last 30 days
end_time = int(time.time())
start_time = end_time - (30 * 24 * 60 * 60)

results = execute_query(
    log_group_name="/aws/lambda/pathogen-detection",
    query_name="perv_detections",
    start_time=start_time,
    end_time=end_time
)

for result in results['results']:
    print(f"Run: {result['run_id']}, Subtypes: {result['subtypes_detected']}")
```

---

#### get_query_by_name(query_name: str, **params) -> str

Get query string with parameter substitution.

**Parameters**:
- `query_name: str` - Query name
- `**params` - Parameters for substitution

**Returns**:
- `str` - Query string with substituted parameters

**Example**:
```python
from lib.audit.cloudwatch_queries import get_query_by_name

query = get_query_by_name("run_audit_trail", run_id="RUN-2025-001")
print(query)  # Query with {run_id} replaced
```

---

## Enums

### PathogenConfidence

```python
class PathogenConfidence(str, Enum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"
    NEGATIVE = "NEGATIVE"
```

### WorkflowStatus

```python
class WorkflowStatus(str, Enum):
    INITIATED = "initiated"
    UPLOADING = "uploading"
    BASECALLING = "basecalling"
    QC_ANALYSIS = "qc_analysis"
    HOST_REMOVAL = "host_removal"
    PATHOGEN_DETECTION = "pathogen_detection"
    QUANTIFICATION = "quantification"
    REPORT_GENERATION = "report_generation"
    COMPLETED = "completed"
    FAILED = "failed"
```

### PERVSubtype

```python
class PERVSubtype(str, Enum):
    PERV_A = "PERV-A"
    PERV_B = "PERV-B"
    PERV_C = "PERV-C"
```

---

## Complete Usage Example

```python
#!/usr/bin/env python3
"""
Complete example: PERV detection with type safety, logging, and audit.
"""

from pathlib import Path
from lib.models.pathogen import (
    PERVTypingOutput,
    PERVDetectionResult,
    PERVSubtype,
    PathogenConfidence
)
from lib.repositories.rds_repository import RDSWorkflowRepository
from lib.models.workflow import WorkflowExecution, WorkflowStatus
from lib.logging.logger import get_logger, AuditLogger

# Initialize logging
logger = get_logger("perv-detection")
audit = AuditLogger("perv-detection")

# Initialize repository
repo = RDSWorkflowRepository(
    cluster_arn=os.environ['RDS_CLUSTER_ARN'],
    secret_arn=os.environ['RDS_SECRET_ARN']
)

def process_perv_detection(run_id: str, bam_file: Path, operator_email: str):
    """Process PERV detection with full type safety and audit logging."""

    logger.info("Starting PERV detection", extra={"run_id": run_id})

    # Update workflow status
    repo.update_status(run_id=run_id, status=WorkflowStatus.PATHOGEN_DETECTION)

    # Perform PERV typing (type-safe)
    result = PERVTypingOutput(
        run_id=run_id,
        bam_file=bam_file,
        detections={
            PERVSubtype.PERV_A: PERVDetectionResult(
                subtype=PERVSubtype.PERV_A,
                reads_aligned=42,
                coverage=0.92,
                specific_motifs_found=['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC'],
                mean_identity=97.5,
                confidence=PathogenConfidence.HIGH  # Auto-validated
            )
        }
    )

    # Audit logging (PMDA compliance)
    if result.requires_sns_alert:
        audit.log_perv_detection(
            run_id=run_id,
            subtypes_detected=[s.value for s in result.positive_subtypes],
            confidence_levels={
                subtype.value: det.confidence.value
                for subtype, det in result.detections.items()
            },
            operator_email=operator_email
        )

        logger.critical(
            "PERV DETECTION ALERT",
            extra=result.to_audit_log()
        )

    # Complete workflow
    repo.update_status(run_id=run_id, status=WorkflowStatus.COMPLETED)

    logger.info("PERV detection complete", extra={"run_id": run_id})

    return result


# Usage
if __name__ == "__main__":
    result = process_perv_detection(
        run_id="RUN-2025-001",
        bam_file=Path("/data/aligned.bam"),
        operator_email="tanaka@example.com"
    )

    print(f"Alert required: {result.requires_sns_alert}")
    print(f"Positive subtypes: {result.positive_subtypes}")
```

---

**Last Updated**: 2025-01-15
**Version**: 2.0.0
