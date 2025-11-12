# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Project**: MinION Pathogen Screening Pipeline - PMDA-compliant xenotransplantation donor pig screening (91 pathogens)

## Critical Context

**⚠️ PERV Detection**: ANY PERV-A/B/C detection → Immediate SNS alert (`scripts/phase4_pathogen/perv_typing.py`)
**⚠️ PMDA Compliance**: Must maintain 91 pathogen coverage, PPA >95%, NPA >98%
**⚠️ Data Security**: No patient data in git, encrypted S3 only
**⚠️ AWS Region**: Always use `ap-northeast-1` (Tokyo)

## Essential Commands

```bash
pytest tests/                      # Run tests
pytest tests/test_pmda_compliance.py::test_name -v  # Single test
black scripts/ lambda/ tools/      # Format
./tools/workflow_cli.py start --run-id RUN-2024-001 --bucket minion-data
```

## Architecture

**7-Phase Pipeline**: phase0_sample_prep → phase1_basecalling → phase2_qc → phase3_host_removal → phase4_pathogen → phase5_quantification → phase6_reports

**AWS Pattern**: S3 upload → Lambda orchestrator → Step Functions → EC2 per phase (auto-terminates)

## Key Files

- Config: `templates/config/pmda_pathogens.json`
- Pipeline: `lambda/orchestration/pipeline_orchestrator.py`
- PERV: `scripts/phase4_pathogen/perv_typing.py`

## Documentation

| Quick Access | Details |
|-------------|---------|
| [Commands](docs/COMMANDS.md) | All commands & operations |
| [Patterns](docs/PATTERNS.md) | Code patterns & conventions |
| [Architecture](docs/ARCHITECTURE.md) | Detailed architecture |
| [Protocols](docs/PROTOCOLS_GUIDE.md) | Sample prep protocols |
| [Development](docs/DEVELOPMENT_GUIDE.md) | Dev guide |

## Recent Updates

**2025-11-13**: Protocol 12 v2.1 - Added Step 2.5 for circular/ssDNA virus support
- TRUE 91/91 pathogen coverage (87/91 → 91/91)
- Enables detection of PCV2, PCV3 (Special Management), TTV, PPV
- Cost: +¥5,000/sample (+3.2%), total ¥162,000
- Time: +2.5 hours (+19%), total 15.5 hours
- Bioinformatics: Reference duplication for circular genomes
- Documentation: 24 files created/modified, 8 Japanese translations

See [Session Log](docs/claude-sessions/2025-11-13-protocol-12-v2.1-circular-ssdna-support.md) for full details.