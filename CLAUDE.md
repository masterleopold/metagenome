# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Project**: MinION Pathogen Screening Pipeline - PMDA-compliant xenotransplantation donor pig screening (91 pathogens)

## Critical Constraints

**⚠️ PERV Detection**: ANY PERV-A/B/C detection → Immediate SNS alert (`scripts/phase4_pathogen/perv_typing.py:34`)
**⚠️ PMDA Compliance**: Must maintain 91 pathogen coverage, PPA >95%, NPA >98% (`templates/config/pmda_pathogens.json`)
**⚠️ Data Security**: No patient data in git, encrypted S3 only
**⚠️ AWS Region**: Always use `ap-northeast-1` (Tokyo)

## Essential Commands

```bash
# Testing
python -m pytest tests/                                    # Run all tests
python -m pytest tests/test_pmda_compliance.py::test_name -v  # Single test
python -m pytest tests/ --cov=scripts --cov=lambda         # With coverage

# Code Quality
black scripts/ lambda/ tools/                              # Format
flake8 scripts/ lambda/                                    # Lint

# Pipeline
./tools/workflow_cli.py start --run-id RUN-001 --bucket minion-data
./tools/workflow_cli.py status --run-id RUN-001 --watch    # Monitor
```

## Architecture

**7-Phase Pipeline**: phase0_sample_prep → phase1_basecalling → phase2_qc → phase3_host_removal → phase4_pathogen → phase5_quantification → phase6_reports

**AWS Pattern (Containerless)**: S3 upload → Lambda orchestrator → Step Functions → EC2 per phase (custom AMIs, auto-terminate)
- No Docker containers - uses custom AMIs with pre-installed tools
- Phase 1: g4dn.xlarge (GPU for Dorado basecalling)
- Phase 3: r5.4xlarge (128GB RAM for host removal)
- Phase 4: 4x parallel EC2 (Kraken2/BLAST pathogen detection)

## Code Structure

```
scripts/
├── phase0_sample_prep/        # Sample routing & workflow determination
├── phase1_basecalling/        # FAST5→FASTQ (Dorado duplex)
├── phase2_qc/                 # Quality metrics (NanoPlot/PycoQC)
├── phase3_host_removal/       # Host depletion (Minimap2)
├── phase4_pathogen/           # Multi-DB screening (Kraken2/BLAST)
│   ├── perv_typing.py         # CRITICAL: PERV-A/B/C subtype detection
│   ├── detect_4viruses.py     # NEW: Hantavirus/Polyomavirus/Spumavirus/EEEV
│   ├── detect_pmda_4viruses.py   # Polyoma/Hantavirus/EEEV/Spumavirus
│   └── detect_pmda_all_91_pathogens.py  # Full 91 pathogen coverage
├── phase5_quantification/     # Abundance calculation (copies/mL)
└── phase6_reports/            # PMDA-compliant report generation

lambda/
├── orchestration/             # Pipeline coordination (S3 trigger)
├── monitoring/                # EC2 lifecycle management
├── phases/                    # Phase-specific handlers
└── shared/                    # Common validation logic

surveillance/                  # NEW: 4-Virus Surveillance System
├── external/                  # External information collectors
│   ├── estat_client.py        # E-Stat API (government statistics)
│   ├── maff_scraper.py        # MAFF surveillance reports
│   └── academic_monitor.py    # PubMed + J-STAGE
├── internal/                  # Internal pipeline integration
│   └── pipeline_listener.py   # Phase 4 result monitoring
├── alerting/                  # Notification system
│   ├── severity_engine.py     # 4-level classification
│   └── notification_router.py # Multi-channel alerts
├── dashboard/                 # Streamlit UI (port 8501)
│   └── app.py
├── api/                       # FastAPI REST API (port 8000)
│   └── main.py
└── lambda/                    # Lambda functions
    ├── external_collector/    # Daily @ 11:00 JST
    └── pipeline_listener/     # Real-time S3 events

tools/
├── workflow_cli.py            # Main CLI (start/status/metrics)
├── database_setup.sh          # Reference database initialization
└── monitoring_dashboard.py    # Streamlit dashboard
```

## Key Patterns

### File Validation (Required)
```python
from pathlib import Path

def process_file(file_path: str):
    file = Path(file_path)
    if not file.exists():
        raise FileNotFoundError(f"File not found: {file}")
    if file.stat().st_size == 0:
        raise ValueError(f"Empty file: {file}")
```

### AWS Operations (Always ap-northeast-1)
```python
import boto3
s3 = boto3.client('s3', region_name='ap-northeast-1')
```

### PERV Detection Pattern
- Check `PERV_MARKERS` dict in `perv_typing.py:14-31`
- Any detection → immediate SNS alert
- Subtype-specific motifs for A/B/C classification

### BAM File Handling
```python
import pysam
bam = pysam.AlignmentFile(str(bam_file), "rb")
# Always validate index exists
if not Path(f"{bam_file}.bai").exists():
    pysam.index(str(bam_file))
```

## Key Configuration Files

- `templates/config/pmda_pathogens.json` - 91 pathogen definitions, Protocol 12 v2.1 specs
- `lambda/orchestration/pipeline_orchestrator.py` - Main orchestration logic
- `scripts/phase4_pathogen/perv_typing.py` - CRITICAL PERV detection
- `infrastructure/terraform/` - AWS infrastructure as code

## Database Locations (EFS Mount)

- Kraken2: `/mnt/efs/databases/kraken2/pmda_2024/`
- BLAST: `/mnt/efs/databases/blast/nt_2024/`
- PERV: `/mnt/efs/databases/perv/references/`
- Custom PMDA: `/mnt/efs/databases/pmda/2024.1/`

## Documentation

| Document | Purpose |
|----------|---------|
| [README](README.md) | Project overview, installation, usage |
| [Quick Reference](docs/QUICK_REFERENCE.md) | Commands, workflows, troubleshooting |
| [Architecture](docs/ARCHITECTURE.md) | Detailed system design |
| [Patterns](docs/PATTERNS.md) | Code conventions & best practices |
| [Protocols](docs/PROTOCOLS_GUIDE.md) | Lab protocols (Protocol 12 v2.1) |
| [Changelog](docs/CHANGELOG.md) | Version history & updates |

## Important Notes

- **Protocol 12 v2.1**: Includes Step 2.5 for circular/ssDNA virus support (PCV2, PCV3, TTV, PPV)
- **4-Virus Surveillance** (v2.2.0): Dual-source monitoring system (external: MAFF/E-Stat/PubMed/J-STAGE + internal: Phase 4 results)
  - Daily collection @ 11:00 JST, real-time internal monitoring
  - 4-level severity (CRITICAL/HIGH/MEDIUM/LOW), multi-channel alerts (SNS/SES/SMS/Slack/Dashboard)
  - **Slack Integration**: Bot API + Webhooks, severity-based channel routing (#critical-alerts, #pathogen-alerts, #pathogen-monitoring)
  - Cost: ~$7.40/month, DynamoDB + Lambda + S3
- **Test Pattern**: Use pytest with moto for AWS service mocking
- **Error Handling**: Always validate file existence/size before processing
- **EC2 Pattern**: Each phase runs on dedicated instance, auto-terminates on completion
- **Spot Instances**: 70% cost savings, handled by orchestrator
- **RDS Schema**: Aurora Serverless v2 for pipeline metadata

## Recent Updates (2025-11-15)

### Slack Notification Integration for 4-Virus Surveillance (v2.2.0)
- **Implementation**: `surveillance/alerting/slack_client.py` - Dual delivery (Bot API + Webhooks)
- **Features**: Rich Block Kit formatting, severity-based channel routing, action buttons for critical alerts
- **Channel Routing**: CRITICAL → #critical-alerts, HIGH → #pathogen-alerts, MEDIUM → #pathogen-monitoring
- **Setup**: `surveillance/docs/SLACK_SETUP.md` (complete guide), `surveillance/.env.template` (credentials)
- **Testing**: `surveillance/tests/test_slack_integration.py` - Connection, alerts, daily summaries
- **Deployment**: `surveillance/scripts/setup_lambda_env.sh` - Automated Lambda environment setup
- **Slack App**: App ID: A09TVLTGDSL, Scopes: `chat:write`, `chat:write.public`, `channels:read`
- **Security**: All credentials via environment variables (SLACK_BOT_TOKEN, SLACK_SIGNING_SECRET)
- **Cost**: Zero additional cost (Slack free tier)
- **Documentation**: Updated ARCHITECTURE.md, CHANGELOG.md, RECENT_UPDATES.md, surveillance README, portal page

### 4-Virus Surveillance System Implementation (v2.1.0) - 2025-11-14
- **Target viruses**: Hantavirus, Polyomavirus, Spumavirus (MHLW Special Management #5), EEEV
- **External collectors**: `surveillance/external/` (MAFF, E-Stat APP_ID: bae1f981a6d093a9676b03c8eea37324b8de421b, PubMed, J-STAGE web scraping)
- **Internal listener**: `surveillance/internal/pipeline_listener.py` monitors Phase 4 results
- **Alerting**: `surveillance/alerting/severity_engine.py` with YAML rules, notification_router.py for multi-channel
- **UIs**: Streamlit dashboard (8501), FastAPI (8000)
- **Infrastructure**: 3 DynamoDB tables, S3 data lake, Lambda functions, EventBridge schedule
- **Quality Assurance**: Ultra-thorough 7-layer audit - 4 bugs found and fixed, ZERO bugs remaining