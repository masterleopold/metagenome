# Quick Reference Guide

**Version**: 2.0 (Updated 2025-01-15)
**New**: Added type-safe patterns, repository pattern, unified logging

> **v2.0 Update**: See [New Patterns](#new-patterns-v20) section for Pydantic models, repository pattern, and logging examples.

## Essential Commands

### Testing
```bash
pytest tests/                              # Run all tests
pytest tests/test_pmda_compliance.py -v   # Run specific test file
pytest tests/test_pmda_compliance.py::test_name -v  # Run single test
pytest --cov=scripts --cov=lambda tests/  # With coverage

# NEW v2.0: Test with new patterns (no AWS needed!)
pytest tests/integration/test_new_patterns.py -v  # Type-safe models + SQLite repos
pytest tests/ --cov=lib --cov-report=html  # Coverage for new lib/
```

### Code Quality
```bash
black scripts/ lambda/ tools/             # Format code
flake8 scripts/ lambda/                   # Lint check
mypy scripts/                             # Type checking
```

### Pipeline Operations
```bash
./tools/workflow_cli.py start --run-id RUN-2024-001 --bucket minion-data
./tools/workflow_cli.py status --run-id RUN-2024-001
./tools/workflow_cli.py logs --run-id RUN-2024-001 --phase 4
```

### Surveillance System Operations
```bash
# Start Streamlit Dashboard (auto-refresh every 30s)
streamlit run surveillance/dashboard/app.py --server.port 8501

# Start FastAPI REST API
cd surveillance/api && uvicorn main:app --host 0.0.0.0 --port 8000 --reload

# Manual test external collectors
python surveillance/external/estat_client.py --test
python surveillance/external/maff_scraper.py --test
python surveillance/external/academic_monitor.py --test

# Check surveillance detections (via API)
curl http://localhost:8000/api/v1/detections?limit=10
curl http://localhost:8000/api/v1/alerts/active

# View DynamoDB surveillance data
aws dynamodb scan --table-name surveillance-detections --region ap-northeast-1 --max-items 5
aws dynamodb scan --table-name surveillance-external-updates --region ap-northeast-1 --max-items 5
```

### AWS Operations
```bash
aws s3 ls s3://minion-data/runs/ --region ap-northeast-1
aws lambda invoke --function-name pipeline-orchestrator --payload '{}' response.json
aws ec2 describe-instances --filters "Name=tag:Pipeline,Values=MinION" --region ap-northeast-1
```

### Database Updates
```bash
# Update Kraken2 database
./scripts/database/update_kraken2.sh

# Update BLAST database
./scripts/database/update_blast_db.sh

# Update PERV reference database
./scripts/database/update_perv_db.sh
```

## Common Workflows

### 1. New Sample Processing
1. Upload FAST5 files to S3: `s3://minion-data/raw/{RUN_ID}/`
2. Start pipeline: `./tools/workflow_cli.py start --run-id {RUN_ID}`
3. Monitor progress: `./tools/workflow_cli.py status --run-id {RUN_ID}`
4. Review reports: `s3://minion-data/reports/{RUN_ID}/`

### 2. PERV Alert Response
1. Check alert details in SNS notification
2. Review typing results: `s3://minion-data/results/{RUN_ID}/phase4/perv_typing.json`
3. Generate detailed report: `./tools/perv_report.py --run-id {RUN_ID}`
4. Notify veterinary team immediately

### 3. Failed Run Recovery
1. Check logs: `./tools/workflow_cli.py logs --run-id {RUN_ID} --phase {PHASE}`
2. Identify failure point in Step Functions console
3. Fix issue (usually EC2 timeout or memory)
4. Restart from phase: `./tools/workflow_cli.py restart --run-id {RUN_ID} --from-phase {PHASE}`

### 4. Database Maintenance
1. Check current versions: `./tools/check_db_versions.sh`
2. Update if needed (monthly): `./scripts/database/update_all.sh`
3. Verify integrity: `./tools/verify_databases.py`
4. Update config: `templates/config/pmda_pathogens.json`

### 5. 4-Virus Surveillance Alert Response
1. Check alert details in dashboard or SNS notification
2. Review detection source: `surveillance/dashboard/app.py` (Active Alerts tab)
3. For CRITICAL alerts (Spumavirus >500 or ANY EEEV):
   - Review Phase 4 results: `s3://minion-data/results/{RUN_ID}/phase4/`
   - Verify external sources: Check MAFF/E-Stat/PubMed/J-STAGE reports
   - Generate detailed report via API: `GET /api/v1/detections?virus={virus_type}`
4. For external source alerts:
   - Review original publication/report in S3: `s3://surveillance-data/external/`
   - Cross-reference with internal detections
   - Update response protocols if needed

## Debug Commands

### Check EC2 Status
```bash
aws ec2 describe-instances \
  --filters "Name=instance-state-name,Values=running" \
  --query "Reservations[].Instances[].[InstanceId,Tags[?Key=='Name'].Value|[0],State.Name,LaunchTime]" \
  --output table
```

### View Lambda Logs
```bash
aws logs tail /aws/lambda/pipeline-orchestrator --follow
```

### S3 Data Verification
```bash
aws s3api head-object --bucket minion-data --key raw/{RUN_ID}/batch_001.fast5
```

### Step Functions Status
```bash
aws stepfunctions describe-execution \
  --execution-arn arn:aws:states:ap-northeast-1:ACCOUNT:execution:MinIONPipeline:{RUN_ID}
```

## Performance Optimization

### EC2 Instance Selection
- Phase 1 (Basecalling): `g4dn.xlarge` (GPU required)
- Phase 3 (Host Removal): `r5.4xlarge` (128GB RAM)
- Phase 4 (Pathogen): `c5.4xlarge` (4x parallel)
- Other phases: `t3.large` (sufficient)

### Cost Optimization
- Use Spot Instances: Add `--spot` flag to workflow_cli.py
- Schedule batch runs: Process multiple samples together
- Clean up S3: Remove intermediate files after 30 days

## Troubleshooting Quick Fixes

| Issue | Quick Fix |
|-------|-----------|
| EC2 timeout | Increase timeout in Lambda config |
| Out of memory | Use larger instance type |
| S3 access denied | Check IAM roles and bucket policies |
| Database not found | Verify EFS mount and paths |
| PERV false positive | Check contamination, rerun with clean sample |

## Important File Locations

### Configuration
- Pipeline config: `templates/config/pmda_pathogens.json`
- AWS config: `terraform/variables.tf`
- Database paths: `scripts/config/database_paths.json`

### Lambda Functions
- Orchestrator: `lambda/orchestration/pipeline_orchestrator.py`
- Phase handlers: `lambda/phase_*/handler.py`

### Analysis Scripts
- PERV typing: `scripts/phase4_pathogen/perv_typing.py`
- Quantification: `scripts/phase5_quantification/calculate_metrics.py`
- Report generation: `scripts/phase6_reports/generate_report.py`

### Databases (EFS Mount)
- Kraken2: `/mnt/efs/databases/kraken2/`
- BLAST: `/mnt/efs/databases/blast/`
- PERV: `/mnt/efs/databases/perv/`

### Surveillance System
- External collectors: `surveillance/external/`
- Internal listeners: `surveillance/internal/`
- Alerting system: `surveillance/alerting/`
- Dashboard: `surveillance/dashboard/app.py`
- REST API: `surveillance/api/main.py`
- Lambda functions: `surveillance/lambda/`
- Configuration: `surveillance/config/config.yaml`, `surveillance/config/severity_rules.yaml`
- Terraform: `infrastructure/surveillance/main.tf`
- DynamoDB tables: `surveillance-detections`, `surveillance-external-updates`, `surveillance-notifications`

## Emergency Contacts

- Pipeline Issues: DevOps team (via Slack #minion-pipeline)
- PERV Detection: Veterinary team (immediate notification required)
- AWS Issues: Cloud team (via PagerDuty)
- Database Updates: Bioinformatics team

---

## New Patterns (v2.0)

**Added**: 2025-01-15
**Full Documentation**: [docs/NEW_PATTERNS_GUIDE.md](NEW_PATTERNS_GUIDE.md)

### Quick Start: Using New Patterns

#### 1. Import Type-Safe Models

```python
# PERV detection
from lib.models.pathogen import PERVTypingOutput, PERVDetectionResult, PERVSubtype

# 91 pathogen screening
from lib.models.pathogen import PMDA91PathogenResult, PathogenDetectionResult

# Workflow tracking
from lib.models.workflow import WorkflowExecution, WorkflowStatus, QCMetrics
```

#### 2. Use Repository Pattern

```python
# Production (RDS)
from lib.repositories.rds_repository import RDSWorkflowRepository

repo = RDSWorkflowRepository(
    cluster_arn=os.environ['RDS_CLUSTER_ARN'],
    secret_arn=os.environ['RDS_SECRET_ARN']
)

workflow = WorkflowExecution(run_id="RUN-001", ...)
repo.create(workflow)

# Testing (SQLite - no AWS needed!)
from lib.repositories.sqlite_repository import SQLiteWorkflowRepository

test_repo = SQLiteWorkflowRepository(db_path=":memory:")
test_repo.create(workflow)  # Same interface!
```

#### 3. Add Unified Logging

```python
# Lambda handler
from lib.logging.logger import get_logger, AuditLogger
from aws_lambda_powertools import Logger, Tracer, Metrics

logger = get_logger("pathogen-detection")
tracer = Tracer(service="pathogen-detection")
metrics = Metrics(namespace="MinION/Pipeline")
audit = AuditLogger(service="pathogen-detection")

@logger.inject_lambda_context()
@tracer.capture_lambda_handler
def lambda_handler(event, context):
    logger.info("Started", extra={"run_id": event['run_id']})

    # PMDA audit logging
    audit.log_perv_detection(
        run_id=event['run_id'],
        subtypes_detected=["PERV-A"],
        confidence_levels={"PERV-A": "HIGH"},
        operator_email=event['operator_email']
    )

    metrics.add_metric(name="PERVDetected", unit=MetricUnit.Count, value=1)
```

#### 4. Query CloudWatch Audit Logs

```python
# Get all PERV detections
from lib.audit.cloudwatch_queries import execute_query
import time

end_time = int(time.time())
start_time = end_time - (30 * 24 * 60 * 60)  # 30 days

results = execute_query(
    log_group_name="/aws/lambda/pathogen-detection",
    query_name="perv_detections",
    start_time=start_time,
    end_time=end_time
)

# Get audit trail for specific run
audit_trail = execute_query(
    log_group_name="/aws/lambda/pathogen-detection",
    query_name="run_audit_trail",
    start_time=start_time,
    end_time=end_time,
    run_id="RUN-2025-001"
)
```

### Available CloudWatch Queries

Run these directly in CloudWatch Logs Insights:

1. `perv_detections` - All PERV alerts (CRITICAL for PMDA)
2. `failed_workflows` - Error analysis by phase
3. `workflow_durations` - Performance monitoring
4. `pathogen_screening` - 91 pathogen results
5. `run_audit_trail` - Complete run history
6. `qc_failures` - QC threshold failures
7. `operator_activity` - User action tracking
8. `high_priority_pathogens` - PERV/HantV/EEEV
9. `database_operations` - DB success rates
10. `phase_performance` - Duration by phase
11. `recent_workflows` - Last 7 days
12. `error_rate_by_phase` - SLA monitoring

### Deploy CloudWatch Dashboard

```bash
aws cloudwatch put-dashboard \
  --dashboard-name MinION-PMDA-Compliance \
  --dashboard-body file://lib/audit/cloudwatch_dashboard.json \
  --region ap-northeast-1
```

### Common Patterns Cheat Sheet

| Task | Old Way (v1.0) | New Way (v2.0) |
|------|----------------|----------------|
| Return detection result | `return {'perv': {...}}` | `return PERVTypingOutput(...)` |
| Database insert | `rds.execute_statement(sql=...)` | `repo.create(workflow)` |
| Logging | `print(f"Starting {run_id}")` | `logger.info("Starting", extra={"run_id": run_id})` |
| Testing database | Mock boto3 RDS | Use `SQLiteWorkflowRepository(":memory:")` |

### File Locations (v2.0 Components)

```
lib/
├── models/                # Pydantic models
│   ├── pathogen.py       # PERV, 91 pathogen models
│   ├── workflow.py       # Workflow, QC models
│   └── database.py       # DB record models
│
├── repositories/          # Repository pattern
│   ├── interfaces.py     # Protocol interfaces
│   ├── rds_repository.py # Production (RDS)
│   └── sqlite_repository.py # Testing (SQLite)
│
├── logging/               # Unified logging
│   ├── logger.py         # Lambda Powertools config
│   └── decorators.py     # Logging decorators
│
└── audit/                 # PMDA compliance
    ├── cloudwatch_queries.py # 12 pre-built queries
    └── cloudwatch_dashboard.json # Dashboard
```

### Example: Update Existing Script

```python
# BEFORE (v1.0)
def identify_perv_type(bam_file):
    """Identify PERV subtypes."""
    return {
        'PERV-A': {'reads': 42, 'coverage': 0.85},
        'PERV-B': {'reads': 0, 'coverage': 0.0}
    }

# AFTER (v2.0)
from lib.models.pathogen import PERVTypingOutput, PERVDetectionResult
from pathlib import Path

def identify_perv_type(bam_file: Path) -> PERVTypingOutput:
    """Identify PERV subtypes with type safety."""
    return PERVTypingOutput(
        run_id=os.environ['RUN_ID'],
        bam_file=bam_file,
        detections={
            PERVSubtype.PERV_A: PERVDetectionResult(
                subtype=PERVSubtype.PERV_A,
                reads_aligned=42,
                coverage=0.85,  # ✅ Auto-validated: 0.0 ≤ coverage ≤ 1.0
                specific_motifs_found=['ATGGCAGCCACCACAGC'],
                mean_identity=96.5,
                confidence=PathogenConfidence.HIGH  # ✅ Auto-calculated
            )
        }
    )

    # Type-safe property access
    if result.requires_sns_alert:
        send_alert(result.to_audit_log())
```

### Documentation Links

- **Full Guide**: [docs/NEW_PATTERNS_GUIDE.md](NEW_PATTERNS_GUIDE.md) (580 lines)
- **Implementation Summary**: [docs/REFACTORING_SUMMARY.md](REFACTORING_SUMMARY.md) (450 lines)
- **Architecture**: [docs/ARCHITECTURE.md#code-quality-architecture-v20](ARCHITECTURE.md#code-quality-architecture-v20)
- **Example Code**: [lambda/phases/trigger_pathogen_detection_v2.py](../lambda/phases/trigger_pathogen_detection_v2.py)
- **Tests**: [tests/integration/test_new_patterns.py](../tests/integration/test_new_patterns.py)