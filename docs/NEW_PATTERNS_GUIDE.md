# MinION Pipeline Refactoring: New Patterns Guide

**Version**: 2.0
**Date**: 2025-01-15
**Status**: ✅ Implemented

---

## Table of Contents

1. [Overview](#overview)
2. [Type Safety with Pydantic](#type-safety-with-pydantic)
3. [Repository Pattern](#repository-pattern)
4. [Unified Logging](#unified-logging)
5. [CloudWatch Audit Queries](#cloudwatch-audit-queries)
6. [Migration Guide](#migration-guide)
7. [Testing](#testing)

---

## Overview

This guide documents the refactored codebase patterns implemented based on expert recommendations. These improvements provide:

- ✅ **Type Safety**: Pydantic models with automatic validation
- ✅ **Database Abstraction**: Repository pattern for testability
- ✅ **PMDA Compliance**: Structured audit logging
- ✅ **Maintainability**: Clear interfaces and separation of concerns

### Key Principles

1. **Explicit is better than implicit** - Use strict types, not `dict` or `Any`
2. **Interface segregation** - Clear input/output contracts for each function
3. **Dependency injection** - Pass repositories, don't create them inline
4. **Audit everything** - PMDA requires full traceability

---

## Type Safety with Pydantic

### Why Pydantic?

- **Automatic validation**: Catches bad data at runtime
- **IDE support**: IntelliSense autocomplete
- **Documentation**: Field descriptions become self-documenting
- **PMDA friendly**: Validation rules → compliance documentation

### Example: PERV Detection

#### ❌ Old Pattern (Weak Typing)

```python
def identify_perv_type(bam_file):
    # Returns dict with unknown structure
    return {
        'PERV-A': {'reads': 42, 'coverage': 0.85},
        'PERV-B': {'reads': 0, 'coverage': 0.0}
    }

# No type checking, easy to break
result = identify_perv_type("test.bam")
if result['PERV-A']['reads'] > 10:  # What if key doesn't exist?
    send_alert()
```

#### ✅ New Pattern (Type-Safe)

```python
from lib.models.pathogen import PERVTypingOutput, PERVDetectionResult

def identify_perv_type(bam_file: Path) -> PERVTypingOutput:
    """Identify PERV subtypes from aligned reads."""
    result = PERVTypingOutput(
        run_id="RUN-001",
        bam_file=bam_file,
        detections={
            PERVSubtype.PERV_A: PERVDetectionResult(
                subtype=PERVSubtype.PERV_A,
                reads_aligned=42,
                coverage=0.85,
                specific_motifs_found=['ATGGCAGCCACCACAGC'],
                mean_identity=96.5,
                confidence=PathogenConfidence.HIGH  # Auto-validated
            )
        }
    )

    # Type-safe access with IDE autocomplete
    if result.requires_sns_alert:  # Property, not dict key
        send_alert(result.to_audit_log())

    return result
```

### Available Models

| Model | Purpose | Location |
|-------|---------|----------|
| `PERVTypingOutput` | PERV detection results | `lib/models/pathogen.py` |
| `PMDA91PathogenResult` | 91 pathogen screening | `lib/models/pathogen.py` |
| `WorkflowExecution` | Pipeline workflow tracking | `lib/models/workflow.py` |
| `QCMetrics` | Quality control metrics | `lib/models/workflow.py` |
| `PathogenDetectionRecord` | Database record | `lib/models/database.py` |

### Creating New Models

```python
from pydantic import BaseModel, Field, field_validator

class MyNewModel(BaseModel):
    """Model description for documentation."""

    run_id: str = Field(
        ...,  # Required field
        pattern=r'^[A-Z0-9][A-Z0-9_-]{0,63}$',
        description="Unique run identifier"
    )

    reads: int = Field(
        ge=0,  # Greater than or equal to 0
        description="Number of reads"
    )

    coverage: float = Field(
        ge=0.0,
        le=1.0,  # Between 0 and 1
        description="Genome coverage (0-1)"
    )

    @field_validator('coverage')
    @classmethod
    def validate_coverage(cls, v: float, info) -> float:
        """Custom validation logic."""
        if v < 0.5 and info.data.get('reads', 0) > 100:
            raise ValueError("High read count but low coverage - suspicious")
        return v

    @property
    def is_valid(self) -> bool:
        """Computed property."""
        return self.reads > 10 and self.coverage > 0.8
```

---

## Repository Pattern

### Why Repository Pattern?

- **Testability**: Use SQLite in tests, RDS in production
- **Abstraction**: Business logic doesn't know about SQL
- **Swappability**: Easy to migrate databases
- **Mockability**: Unit tests don't need real databases

### Architecture

```
┌─────────────────────┐
│ Business Logic      │
│ (Lambda, Scripts)   │
└──────────┬──────────┘
           │ Depends on interface
           ▼
┌─────────────────────┐
│ Repository Protocol │ ← Interface (typing.Protocol)
│ (lib/repositories/  │
│  interfaces.py)     │
└──────────┬──────────┘
           │ Implemented by
    ┌──────┴──────┐
    ▼             ▼
┌─────────┐  ┌──────────┐
│   RDS   │  │  SQLite  │
│  Repo   │  │   Repo   │
│  (Prod) │  │  (Tests) │
└─────────┘  └──────────┘
```

### Example: Workflow Repository

#### Using in Lambda Handler

```python
from lib.models.workflow import WorkflowExecution, WorkflowStatus
from lib.repositories.rds_repository import RDSWorkflowRepository

def lambda_handler(event, context):
    # Initialize repository (production)
    repo = RDSWorkflowRepository(
        cluster_arn=os.environ['RDS_CLUSTER_ARN'],
        secret_arn=os.environ['RDS_SECRET_ARN']
    )

    # Create workflow
    workflow = WorkflowExecution(
        run_id=event['run_id'],
        sample_id=event['sample_id'],
        status=WorkflowStatus.INITIATED,
        operator_name=event['operator_name'],
        operator_email=event['operator_email']
    )

    repo.create(workflow)  # Type-safe, mockable

    # Update status
    repo.update_status(
        run_id=event['run_id'],
        status=WorkflowStatus.COMPLETED
    )

    # Query
    recent = repo.list_recent(limit=50, days=7)
    return {'recent_workflows': len(recent)}
```

#### Using in Tests

```python
from lib.repositories.sqlite_repository import SQLiteWorkflowRepository

def test_workflow_creation():
    # Use in-memory SQLite (no AWS needed!)
    repo = SQLiteWorkflowRepository(db_path=":memory:")

    workflow = WorkflowExecution(
        run_id="TEST-001",
        sample_id="SAMPLE-001",
        status=WorkflowStatus.INITIATED,
        operator_name="Test User",
        operator_email="test@example.com"
    )

    run_id = repo.create(workflow)

    retrieved = repo.get("TEST-001")
    assert retrieved.status == WorkflowStatus.INITIATED
```

### Available Repositories

| Repository | Interface | RDS Impl | SQLite Impl | Purpose |
|------------|-----------|----------|-------------|---------|
| Workflow | `WorkflowRepository` | ✅ | ✅ | Workflow tracking |
| Pathogen Detection | `PathogenDetectionRepository` | ✅ | ✅ | Detection results |
| Surveillance | `DynamoDBSurveillanceRepository` | N/A | N/A | Real-time alerts |

---

## Unified Logging

### Why AWS Lambda Powertools?

- **Structured JSON**: CloudWatch Logs Insights compatible
- **Correlation IDs**: Track requests across services
- **X-Ray Integration**: Distributed tracing
- **PMDA Audit**: Automatic audit trail

### Basic Usage

```python
from aws_lambda_powertools import Logger, Tracer, Metrics
from aws_lambda_powertools.logging import correlation_paths
from aws_lambda_powertools.metrics import MetricUnit

# Initialize (singleton pattern)
logger = Logger(service="pathogen-detection")
tracer = Tracer(service="pathogen-detection")
metrics = Metrics(namespace="MinION/Pipeline", service="pathogen-detection")

@logger.inject_lambda_context(correlation_id_path=correlation_paths.EVENT_BRIDGE)
@tracer.capture_lambda_handler
@metrics.log_metrics(capture_cold_start_metric=True)
def lambda_handler(event, context):
    run_id = event['run_id']

    logger.info(
        "Pathogen detection started",
        extra={
            "run_id": run_id,
            "databases": ["kraken2", "blast", "pmda"],
            "instance_type": "r5.4xlarge"
        }
    )

    try:
        result = detect_pathogens(run_id)

        metrics.add_metric(
            name="PathogensDetected",
            unit=MetricUnit.Count,
            value=len(result.detections)
        )

        logger.info(
            "Pathogen detection complete",
            extra={
                "run_id": run_id,
                "positive_count": len(result.positive_detections)
            }
        )

    except Exception as e:
        logger.exception(
            "Pathogen detection failed",
            extra={
                "run_id": run_id,
                "error_type": type(e).__name__
            }
        )
        raise

    return {"status": "success"}
```

### Audit Logging (PMDA Compliance)

```python
from lib.logging.logger import AuditLogger

audit_logger = AuditLogger(service="pathogen-detection")

# Log PERV detection (CRITICAL)
audit_logger.log_perv_detection(
    run_id="RUN-001",
    subtypes_detected=["PERV-A"],
    confidence_levels={"PERV-A": "HIGH"},
    operator_email="operator@example.com"
)

# Log workflow start
audit_logger.log_workflow_start(
    run_id="RUN-001",
    operator_email="operator@example.com",
    sample_id="SAMPLE-001"
)

# Log errors
audit_logger.log_error(
    run_id="RUN-001",
    phase="pathogen_detection",
    error_message="Kraken2 database not found",
    error_type="FileNotFoundError",
    operator_email="operator@example.com"
)
```

### Log Output (Structured JSON)

```json
{
  "timestamp": "2025-01-15T10:30:45.123Z",
  "level": "INFO",
  "service": "pathogen-detection",
  "message": "Pathogen detection complete",
  "run_id": "RUN-001",
  "positive_count": 2,
  "correlation_id": "abc-123-def-456",
  "cold_start": false,
  "function_name": "pathogen-detection",
  "function_memory_size": 1024,
  "xray_trace_id": "1-5e1c1234-5678abcd"
}
```

---

## CloudWatch Audit Queries

### Pre-Built Queries

Located in: `lib/audit/cloudwatch_queries.py`

#### Query 1: All PERV Detections

```python
from lib.audit.cloudwatch_queries import execute_query
import time

end_time = int(time.time())
start_time = end_time - (30 * 24 * 60 * 60)  # Last 30 days

results = execute_query(
    log_group_name="/aws/lambda/pathogen-detection",
    query_name="perv_detections",
    start_time=start_time,
    end_time=end_time
)

for result in results['results']:
    print(f"PERV detected in {result['run_id']} at {result['timestamp']}")
```

#### Query 2: Workflow Audit Trail

```python
# Get complete audit trail for specific run
results = execute_query(
    log_group_name="/aws/lambda/pathogen-detection",
    query_name="run_audit_trail",
    start_time=start_time,
    end_time=end_time,
    run_id="RUN-2025-001"  # Parameter substitution
)
```

### CloudWatch Dashboard

Deploy the dashboard:

```bash
aws cloudwatch put-dashboard \
  --dashboard-name MinION-PMDA-Compliance \
  --dashboard-body file://lib/audit/cloudwatch_dashboard.json \
  --region ap-northeast-1
```

Provides:
- PERV detection alerts
- 91 pathogen screening results
- Error analysis by phase
- Performance metrics
- Operator activity

---

## Migration Guide

### Phase 1: Adopt Models (Low Risk)

**Timeframe**: 1-2 weeks
**Effort**: Low
**Impact**: High (type safety)

1. **Update Phase 4 PERV detection**:

```python
# Old: scripts/phase4_pathogen/perv_typing.py
def identify_perv_type(bam_file) -> dict:
    return {'PERV-A': {...}}

# New:
from lib.models.pathogen import PERVTypingOutput

def identify_perv_type(bam_file: Path) -> PERVTypingOutput:
    return PERVTypingOutput(
        run_id=os.environ['RUN_ID'],
        bam_file=bam_file,
        detections={...}
    )
```

2. **Update all Phase 4 detection scripts**
3. **Run tests**: `pytest tests/integration/test_new_patterns.py`

### Phase 2: Add Logging (Medium Risk)

**Timeframe**: 3-5 days
**Effort**: Low
**Impact**: Critical (PMDA compliance)

1. **Update Lambda handlers** (start with one):

```python
# Add to lambda/phases/trigger_pathogen_detection.py
from aws_lambda_powertools import Logger, Tracer

logger = Logger(service="pathogen-detection")
tracer = Tracer(service="pathogen-detection")

@logger.inject_lambda_context()
@tracer.capture_lambda_handler
def lambda_handler(event, context):
    logger.info("Phase 4 started", extra={"run_id": event['run_id']})
    ...
```

2. **Deploy dashboard**: Upload `lib/audit/cloudwatch_dashboard.json`
3. **Test queries**: Run CloudWatch Insights queries
4. **Repeat for remaining Lambda handlers** (7 total)

### Phase 3: Repository Pattern (High Risk)

**Timeframe**: 1-2 weeks
**Effort**: Medium
**Impact**: Medium (testability)

1. **Add repository alongside existing code** (dual writes):

```python
# Keep existing SQL working
old_create_workflow(run_id, ...)

# Also write to repository
repo = RDSWorkflowRepository(...)
repo.create(workflow)
```

2. **Run integration tests** comparing outputs
3. **Monitor for 1-2 weeks**
4. **Remove old SQL once validated**

---

## Testing

### Unit Tests (Fast, No AWS)

```bash
# Test models only
pytest tests/unit/test_models.py -v

# Test with SQLite repositories (no RDS needed)
pytest tests/integration/test_new_patterns.py -v
```

### Integration Tests (Requires AWS)

```bash
# Full pipeline test with RDS
pytest tests/integration/test_rds_repositories.py -v \
  --rds-cluster-arn $RDS_CLUSTER_ARN \
  --rds-secret-arn $RDS_SECRET_ARN
```

### Coverage

```bash
pytest tests/ --cov=lib --cov-report=html
open htmlcov/index.html
```

---

## Quick Reference

### Import Cheat Sheet

```python
# Models
from lib.models.pathogen import (
    PERVTypingOutput,
    PMDA91PathogenResult,
    PathogenDetectionResult,
)
from lib.models.workflow import (
    WorkflowExecution,
    WorkflowStatus,
    QCMetrics,
)

# Repositories
from lib.repositories.rds_repository import (
    RDSWorkflowRepository,
    RDSPathogenDetectionRepository,
)
from lib.repositories.sqlite_repository import (
    SQLiteWorkflowRepository,  # For tests
)

# Logging
from lib.logging.logger import (
    get_logger,
    get_tracer,
    get_metrics,
    AuditLogger,
)

# Audit queries
from lib.audit.cloudwatch_queries import (
    get_query_by_name,
    execute_query,
)
```

### Common Patterns

| Task | Old Way | New Way |
|------|---------|---------|
| Return detection result | `return {'perv': {...}}` | `return PERVTypingOutput(...)` |
| Database insert | `rds.execute_statement(sql=...)` | `repo.create(workflow)` |
| Logging | `print(f"Starting {run_id}")` | `logger.info("Starting", extra={"run_id": run_id})` |
| Testing | Mock boto3 RDS | Use `SQLiteWorkflowRepository` |

---

## Support

- **Documentation**: This file + inline docstrings
- **Examples**: `lambda/phases/trigger_pathogen_detection_v2.py`
- **Tests**: `tests/integration/test_new_patterns.py`
- **Questions**: See CLAUDE.md for pipeline-specific guidance

---

**Last Updated**: 2025-01-15
**Version**: 2.0.0
**Author**: Claude Code Refactoring Team
