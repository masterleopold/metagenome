# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Project**: MinION Pathogen Screening Pipeline for PMDA xenotransplantation donor screening (91 pathogens)

## Critical Constraints

⚠️ **PERV Detection**: ANY PERV-A/B/C detection triggers immediate SNS alert (`scripts/phase4_pathogen/perv_typing.py`)
⚠️ **PMDA Compliance**: 91 pathogen coverage, PPA >95%, NPA >98% (`templates/config/pmda_pathogens.json`)
⚠️ **J-STAGE ToS**: 24h max retention - aggregated stats only (`surveillance/docs/JSTAGE_COMPLIANCE.md`)
⚠️ **AWS Region**: ALWAYS `ap-northeast-1` (Tokyo) - hardcoded for regulatory compliance
⚠️ **Security**: NO patient data in git, encrypted S3 only

## Commands

```bash
# Testing
pytest tests/                                    # All tests
pytest tests/integration/test_new_patterns.py -v # v2.0 patterns (Pydantic + Repository)
pytest tests/ --cov=scripts --cov=lambda        # With coverage

# Code quality
black scripts/ lambda/ tools/ lib/              # Format
flake8 scripts/ lambda/ lib/                    # Lint

# Pipeline execution
./tools/workflow_cli.py start --run-id RUN-001 --bucket minion-data --input-prefix runs/RUN-001/fast5/
./tools/workflow_cli.py status --run-id RUN-001 --watch
./tools/workflow_cli.py metrics --run-id RUN-001

# Surveillance dashboard
cd surveillance/dashboard && streamlit run app.py  # Port 8501
cd surveillance/api && uvicorn main:app           # Port 8000
```

## Architecture Overview

**7-Phase Containerless Pipeline**: Lambda orchestrates EC2 instances with custom AMIs (no Docker)

```
S3 Upload → Lambda Orchestrator → Step Functions → EC2 per phase (auto-terminate)
                                                      ↓
Phase 0: Sample Prep Routing (t3.small)
Phase 1: Basecalling (g4dn.xlarge, GPU Dorado duplex)
Phase 2: QC (t3.large, NanoPlot/PycoQC/RIN)
Phase 3: Host Removal (r5.4xlarge, 128GB RAM, Minimap2)
Phase 4: Pathogen Detection (4x parallel EC2, Kraken2/BLAST/PMDA)
Phase 5: Quantification (t3.large, copies/mL)
Phase 6: PMDA Reports (t3.large, PDF + JSON)
```

**Key AWS Services**:
- Lambda + Step Functions for orchestration
- EC2 Spot Instances (70% cost savings)
- EFS for databases (`/mnt/efs/databases/kraken2`, `/blast`, `/perv`)
- RDS Aurora Serverless v2 (PostgreSQL) for metadata
- S3 for data storage with lifecycle policies
- SNS for PERV alerts, SES for reports

## Code Architecture (v2.0)

**Type-Safe Data Models** (`lib/models/`):
- Pydantic v2 models with auto-validation
- `pathogen.py`: PERV detection, 91 pathogens, confidence levels
- `workflow.py`: Workflow execution, QC metrics, phase status
- `database.py`: Database records for repositories

**Repository Pattern** (`lib/repositories/`):
- `interfaces.py`: Protocol-based contracts (typing.Protocol)
- `rds_repository.py`: Production (RDS Data API, PostgreSQL)
- `sqlite_repository.py`: Testing (local SQLite)
- `dynamodb_repository.py`: Surveillance system state

**Unified Logging** (`lib/logging/`):
- AWS Lambda Powertools (Logger, Tracer, Metrics)
- Structured JSON logs for CloudWatch Insights
- 12 pre-built audit queries for PMDA compliance

## Critical Implementation Patterns

### File Validation (REQUIRED)
```python
from pathlib import Path

file = Path(file_path)
if not file.exists():
    raise FileNotFoundError(f"File not found: {file}")
if file.stat().st_size == 0:
    raise ValueError(f"Empty file: {file}")
```

### BAM File Indexing (REQUIRED)
```python
import pysam
from pathlib import Path

bam_path = Path(bam_file)
if not Path(f"{bam_path}.bai").exists():
    pysam.index(str(bam_path))
```

### AWS Client Initialization (REQUIRED)
```python
import boto3

# ALWAYS use Tokyo region for regulatory compliance
s3 = boto3.client('s3', region_name='ap-northeast-1')
ec2 = boto3.client('ec2', region_name='ap-northeast-1')
```

### Type-Safe Workflows (v2.0)
```python
from lib.models.workflow import WorkflowExecution, WorkflowStatus
from lib.repositories.rds_repository import RDSWorkflowRepository

# Create workflow with validation
workflow = WorkflowExecution(
    run_id="RUN-001",
    status=WorkflowStatus.PENDING,
    input_bucket="minion-data",
    input_prefix="runs/RUN-001/fast5/"
)

# Repository pattern with automatic backend selection
repo = RDSWorkflowRepository()  # Production: RDS PostgreSQL
# repo = SQLiteWorkflowRepository()  # Testing: Local SQLite

repo.create(workflow)
repo.update_status("RUN-001", WorkflowStatus.RUNNING)
```

### PERV Detection with Type Safety (v2.0)
```python
from lib.models.pathogen import PERVDetectionResult, PERVSubtype, PathogenConfidence

result = PERVDetectionResult(
    subtype=PERVSubtype.PERV_A,
    reads_aligned=15,
    coverage=0.92,
    specific_motifs_found=['ATGGCAGCCACCACAGC'],
    mean_identity=96.5,
    confidence=PathogenConfidence.HIGH
)

# Auto-validation: requires_alert property
if result.requires_alert:  # True for confidence >= MEDIUM
    send_sns_alert(result)
```

## Critical Files

| File | Purpose | DO NOT MODIFY WITHOUT |
|------|---------|----------------------|
| `scripts/phase4_pathogen/perv_typing.py` | PERV detection logic | PMDA validation |
| `templates/config/pmda_pathogens.json` | 91 pathogen definitions | Regulatory review |
| `lambda/orchestration/pipeline_orchestrator.py` | Main workflow orchestration | Integration testing |
| `lib/models/pathogen.py` | Type-safe PERV/pathogen models | Breaking change review |
| `lib/repositories/interfaces.py` | Repository contracts | All implementations |

## Surveillance System (4-Virus v2.3.0)

**Target Viruses**: Hantavirus, Polyomavirus, Spumavirus, EEEV

**Components**:
- `surveillance/external/`: MAFF scraper, E-Stat API, PubMed/J-STAGE monitor
- `surveillance/internal/`: Real-time Phase 4 result listener
- `surveillance/alerting/`: Severity engine (4 levels), Slack Bot API + Webhooks
- `surveillance/lambda/external_collector/`: Daily at 11:00 JST (cron: 0 2 * * ? *)

**Slack Integration**:
- Channels: #critical-alerts, #pathogen-alerts, #pathogen-monitoring
- Rich Block Kit formatting with action buttons
- Daily summary notifications

## Testing Strategy

**Integration Tests** (`tests/integration/test_new_patterns.py`):
- Validates Pydantic models with edge cases
- Tests repository pattern with SQLite backend
- Verifies type-safe workflows end-to-end

**Mocking AWS Services**:
```python
from moto import mock_s3, mock_lambda, mock_dynamodb

@mock_s3
@mock_dynamodb
def test_pipeline():
    # AWS services are mocked
    pass
```

## Documentation

- [ARCHITECTURE.md](docs/ARCHITECTURE.md) - System design, v2.0 improvements
- [QUICK_REFERENCE.md](docs/QUICK_REFERENCE.md) - Commands, troubleshooting
- [PATTERNS.md](docs/PATTERNS.md) - Code conventions
- [NEW_PATTERNS_GUIDE.md](docs/NEW_PATTERNS_GUIDE.md) - v2.0 patterns (Pydantic, Repository, Logging)
- [API_REFERENCE_V2.md](docs/API_REFERENCE_V2.md) - v2.0 API docs
- [RECENT_UPDATES.md](docs/RECENT_UPDATES.md) - Changelog
- [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md) - Detailed reference
- [docs/grants/](docs/grants/) - NVIDIA Academic Grant (DGX Spark ARM deployment)
