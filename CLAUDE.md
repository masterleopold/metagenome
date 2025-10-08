# CLAUDE.md - Quick Reference for Claude Code

This file provides essential guidance to Claude Code (claude.ai/code) when working with this repository.

## Project Overview

**MinION Pathogen Screening Pipeline** - Production-ready xenotransplantation donor pig screening system detecting 91 PMDA-designated pathogens using Oxford Nanopore MinION and AWS infrastructure.

**Tech Stack**: Python, AWS (Lambda, EC2, Step Functions, S3), Terraform, Next.js (docs portal)
**Target**: Yucatan miniature pigs with 69 genomic edits (3KO, 7TG, 59 PERV inactivations)
**Sample Type**: Plasma cfDNA/cfRNA
**Regulatory**: PMDA xenotransplantation guidelines compliance

## Quick Command Reference

```bash
# Setup
pip install -r requirements.txt
aws configure  # Set AWS_REGION=ap-northeast-1

# Testing
python -m pytest tests/                      # All tests
python -m pytest tests/ -k "not integration" # Unit tests only
python tests/test_pmda_compliance.py         # PMDA compliance

# Pipeline Operations
./tools/workflow_cli.py start --run-id RUN-2024-001 --bucket minion-data --input-prefix runs/RUN-2024-001/fast5/
./tools/workflow_cli.py status --run-id RUN-2024-001 --watch
./tools/workflow_cli.py metrics --run-id RUN-2024-001

# Code Quality
black scripts/ lambda/ lib/ tools/ tests/
flake8 scripts/ lambda/ lib/ tools/ tests/
mypy lib/ scripts/

# Documentation Portal
cd docs-portal
npm run dev   # http://localhost:3003
npm run build # Production build
```

## Critical Project Context

### PERV Detection (Highest Priority)
- **PERV-A/B/C** subtypes must be detected and quantified
- Immediate SNS alert on ANY PERV detection
- Scripts: `perv_typing.py`, `detect_recombinants.py`, `perv_phylogenetics.py`
- Must report in copies/mL plasma

### PMDA 91 Pathogens
Core regulatory requirement. Full list in `md/厚労省異種移植指針_91病原体リスト.md`
- Detection requirements: PPA >95%, NPA >98%, R² >0.90
- Custom database at `/mnt/efs/databases/pmda/2024.1/`

### Pipeline Architecture (6 Phases)
**Implementation**: Lambda functions orchestrate EC2 instances (custom AMIs) for each phase. **No Docker containers used**.

1. **Basecalling** - FAST5→FASTQ (Dorado, GPU g4dn.xlarge)
2. **QC** - Quality assessment (NanoPlot/PycoQC, t3.large)
3. **Host Removal** - Sus scrofa depletion (Minimap2, r5.4xlarge)
4. **Pathogen Detection** - Multi-database screening (4 parallel EC2 instances)
5. **Quantification** - Abundance calculation (t3.large)
6. **Reporting** - PMDA-compliant reports (t3.large)

**Key Architecture Decisions**:
- Lambda functions trigger EC2 instances with UserData scripts
- Custom AMIs pre-installed with analysis tools
- EFS for shared reference databases (Kraken2, BLAST, PERV)
- EC2 instances auto-terminate after completion
- Spot Instances for 70% cost savings

## Key Files & Directories

```
scripts/          # Pipeline phase scripts (1-6)
lambda/           # AWS Lambda functions
lib/              # Shared Python libraries
tools/            # CLI tools and utilities
tests/            # Test suite
infrastructure/   # Terraform IaC
md/               # Japanese strategic documents
docs-portal/      # Next.js documentation site
```

## Current Implementation Phase

**Phase 2** (Years 1-3): Analysis In-house Implementation
- In-house: Bioinformatics analysis
- Outsourced: Sequencing only
- Target: 12-24 analyses/year

## Essential Configuration

```yaml
# templates/config/default_pipeline.yaml
pipeline:
  phases:
    basecalling:
      min_qscore: 9  # Q30 accuracy target
    qc:
      min_reads: 10000000  # 10M minimum
      min_q30: 0.85        # 85% Q30
    pathogen:
      confidence_threshold: 0.1
      pmda_alert_enabled: true
```

## Documentation Index

### Detailed Guides
- [Development Workflow](docs/development/WORKFLOW_GUIDE.md) - Complete setup, commands, architecture
- [Coding Standards](docs/development/CODING_STANDARDS.md) - Style guides, commit conventions
- [API Documentation](docs/API_DOCUMENTATION.md) - REST endpoints, authentication
- [Deployment Guide](docs/DEPLOYMENT_GUIDE.md) - AWS infrastructure setup
- [Session History](docs/claude-sessions/README.md) - Development session logs

### Strategic Documents (Japanese)
- `異種移植用ドナーブタの病原体検査体制構築計画書.md` - Master construction plan
- `専門家向け詳細報告書：異種移植用ドナーブタにおける病原体メタゲノム解析体制の構築.md` - Expert report
- `NGS中心型病原体検査システム_最適化戦略v2.md` - NGS optimization strategy
- `MinION_Protocol_00_目次とマスタードキュメント.md` - Lab protocol master index (10 protocols + 3 appendices)

## Critical Warnings

⚠️ **PERV Detection** - Any detection triggers immediate alert
⚠️ **Data Security** - No patient data in git, use encrypted S3 only
⚠️ **PMDA Compliance** - All changes must maintain regulatory compliance
⚠️ **Host Depletion** - Must achieve >90% removal efficiency
⚠️ **Cost Threshold** - Alert if analysis exceeds $400

## Quick Troubleshooting

| Issue | Solution |
|-------|----------|
| Basecalling timeout | Check GPU, increase `max_runtime_hours`, use g4dn.2xlarge |
| High host contamination | Review extraction protocol, check Minimap2 params |
| No pathogen detection | Verify databases, check thresholds (0.1 default) |
| EC2 launch failure | Check spot availability, fallback to on-demand |

## Contact & Support

- **GitHub Issues**: https://github.com/masterleopold/metagenome/issues
- **Lead Developer**: Yoichiro Hara
- **AWS Region**: ap-northeast-1 (Tokyo)
- **Environment**: Production

## Recently Updated

- 2025-10-08: **MinION Protocol Appendices Complete** - Added 3 missing appendices to lab protocols
  - Appendix A: Reagents & Equipment List (¥21.5M initial investment)
  - Appendix B: Time & Cost Estimates (¥127k/sample, 4-year ROI)
  - Appendix C: Troubleshooting Guide (25 common issues)
- 2025-10-08: **Architecture Documentation Update** - Updated all core documentation to reflect Lambda + EC2 custom AMI architecture (containerless, no Docker) [`0a7eb04`](https://github.com/masterleopold/metagenome/commit/0a7eb049d1006595253573076bb01dd2d0979885)
- 2025-10-08: **CLAUDE.md optimized** - Reduced from 47KB to 5.3KB (88.5% reduction)
- 2025-10-08: Documentation portal with Linear-inspired design (#0089A7)
- 2025-10-08: Vercel deployment configuration with typed routes

For complete historical changes, see [Session History](docs/claude-sessions/README.md).