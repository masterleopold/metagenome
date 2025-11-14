# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Project**: MinION Pathogen Screening Pipeline - PMDA xenotransplantation screening (91 pathogens)

## Critical Constraints

⚠️ **PERV Detection**: ANY PERV-A/B/C → Immediate SNS alert (`scripts/phase4_pathogen/perv_typing.py:34`)
⚠️ **PMDA Compliance**: 91 pathogen coverage, PPA >95%, NPA >98% (`templates/config/pmda_pathogens.json`)
⚠️ **Data Security**: No patient data in git, encrypted S3 only
⚠️ **AWS Region**: Always use `ap-northeast-1` (Tokyo)

## Quick Commands

```bash
# Testing
pytest tests/                           # Run all tests
pytest tests/ --cov=scripts --cov=lambda # With coverage

# Code Quality
black scripts/ lambda/ tools/           # Format
flake8 scripts/ lambda/                 # Lint

# Pipeline
./tools/workflow_cli.py start --run-id RUN-001 --bucket minion-data
./tools/workflow_cli.py status --run-id RUN-001 --watch
```

## Critical Files

| File | Purpose |
|------|---------|
| `scripts/phase4_pathogen/perv_typing.py` | PERV detection - CRITICAL |
| `templates/config/pmda_pathogens.json` | 91 pathogen definitions |
| `lambda/orchestration/pipeline_orchestrator.py` | Main orchestration |
| `surveillance/alerting/slack_client.py` | Slack notifications |
| `surveillance/config/severity_rules.yaml` | Alert severity rules |

## Key Patterns

```python
# ALWAYS validate files
from pathlib import Path
file = Path(file_path)
if not file.exists(): raise FileNotFoundError(f"File not found: {file}")
if file.stat().st_size == 0: raise ValueError(f"Empty file: {file}")

# ALWAYS use Tokyo region
s3 = boto3.client('s3', region_name='ap-northeast-1')

# ALWAYS check BAM index
if not Path(f"{bam_file}.bai").exists():
    pysam.index(str(bam_file))
```

## Database Paths (EFS)

- Kraken2: `/mnt/efs/databases/kraken2/pmda_2024/`
- BLAST: `/mnt/efs/databases/blast/nt_2024/`
- PERV: `/mnt/efs/databases/perv/references/`

## Documentation

| Doc | Location | Purpose |
|-----|----------|---------|
| Architecture | [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | System design |
| Commands | [docs/QUICK_REFERENCE.md](docs/QUICK_REFERENCE.md) | All commands |
| Updates | [docs/RECENT_UPDATES.md](docs/RECENT_UPDATES.md) | Latest changes |
| Patterns | [docs/PATTERNS.md](docs/PATTERNS.md) | Code conventions |
| Full Reference | [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md) | Detailed info |

## Current Features

- **7-Phase Pipeline**: Containerless AWS (Lambda → EC2 AMIs)
- **4-Virus Surveillance** (v2.2.0): Hantavirus/Polyomavirus/Spumavirus/EEEV
  - Slack alerts to #critical-alerts, #pathogen-alerts, #pathogen-monitoring
  - External: MAFF/E-Stat/PubMed/J-STAGE (daily 11:00 JST)
  - Internal: Real-time Phase 4 monitoring
- **Protocol 12 v2.1**: Circular/ssDNA virus support (PCV2, PCV3, TTV, PPV)

**See [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md) for detailed architecture, directory structure, and implementation notes.**