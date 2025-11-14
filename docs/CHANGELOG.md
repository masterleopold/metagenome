# Project Changelog

## Recent Updates

### 2025-11-14: 4-Virus Surveillance System Implementation
**Comprehensive dual-source monitoring system for 4 target viruses**
- **Target Viruses**: Hantavirus, Polyomavirus, Spumavirus (MHLW Special Management #5), EEEV
- **Architecture**: Standalone monitoring system with API integration to main pipeline
  - **External Sources**: Daily collection (11:00 JST) from MAFF, E-Stat (APP_ID: bae1f981a6d093a9676b03c8eea37324b8de421b), PubMed, J-STAGE
  - **Internal Pipeline**: Real-time Phase 4 result monitoring via S3 events
- **Components Implemented**:
  - External collectors: `surveillance/external/` (MAFF scraper, E-Stat API, Academic monitor with J-STAGE web scraping)
  - Internal listeners: `surveillance/internal/pipeline_listener.py`
  - Severity engine: 4-level classification (CRITICAL/HIGH/MEDIUM/LOW) with YAML rules
  - Notification router: Multi-channel (SNS/SES/SMS/Dashboard/API)
  - Streamlit dashboard: Real-time monitoring (port 8501, 30s auto-refresh)
  - FastAPI REST API: Programmatic access (port 8000)
- **AWS Infrastructure**:
  - 3 DynamoDB tables (detections, external-updates, notifications)
  - S3 data lake with lifecycle policies (365d external, 730d internal)
  - Lambda functions (external_collector, pipeline_listener, alert_processor)
  - EventBridge daily schedule (cron: 0 2 * * ? * = 11:00 JST)
  - SNS topics (critical-alerts, high-alerts, daily-summary)
- **Cost Efficiency**: ~$7.40/month operating cost (Lambda ~$2, DynamoDB ~$5, S3 ~$0.30, SNS/SES ~$0.10)
- **Documentation**:
  - Added comprehensive architecture section to `docs/ARCHITECTURE.md` (300+ lines)
  - Created docs-portal page: `docs-portal/src/app/surveillance/page.tsx`
  - Updated main `README.md`, sidebar navigation, homepage
  - Created `surveillance/README.md` user guide
  - Updated `docs/QUICK_REFERENCE.md`, `docs/DEVELOPMENT_GUIDE.md`, `docs/RECENT_UPDATES.md`
- **Session Log**: See [detailed log](claude-sessions/2025-11-14-4virus-surveillance-system.md)

### 2025-11-13: Protocol 12 v2.1 - Circular/ssDNA Virus Support
**Added Step 2.5 for complete pathogen coverage**
- **Coverage Achievement**: TRUE 91/91 pathogen coverage (87/91 → 91/91)
- **New Capability**: Enables detection of PCV2, PCV3 (Special Management), TTV, PPV
- **Cost Impact**: +¥5,000/sample (+3.2%), total ¥162,000
- **Time Impact**: +2.5 hours (+19%), total 15.5 hours
- **Technical Implementation**: Reference duplication for circular genomes in bioinformatics pipeline
- **Documentation**: 24 files created/modified, 8 Japanese translations completed
- **Session Log**: See [detailed log](claude-sessions/2025-11-13-protocol-12-v2.1-circular-ssdna-support.md)

### 2025-11-12: CLAUDE.md Optimization and Enhancement
**Optimized documentation structure for improved Claude Code performance**
- **CLAUDE.md Optimization**: Reduced from 181 lines (12KB) to 86 lines (4KB), then enhanced to 191 lines (8KB) with critical architecture
  - Removed 53% of redundant content (historical updates, detailed specs)
  - Created separate documentation files for better organization
  - Enhanced with essential architectural patterns and code conventions
- **New Documentation Files**:
  - `docs/CHANGELOG.md` - Historical updates and project timeline (this file)
  - `docs/TECHNICAL_DETAILS.md` - Complete architecture, phases, databases, performance requirements
  - `docs/TROUBLESHOOTING.md` - Common issues, solutions, and AWS debugging commands
- **CLAUDE.md Enhanced Sections**:
  - Architecture Overview: Lambda-triggered EC2 orchestration pattern (no containers)
  - Code Conventions: Python import patterns, error handling, test patterns, PERV markers
  - Key Technical Details: AWS infrastructure, data flow, important paths
  - Development Commands: Complete setup, testing, and deployment commands
- **Backup Created**: `CLAUDE.md.backup.20251112-221227` for reference
- **Technical Q&A**: Added `docs/technical-qa/MinION-DNA-Library-Preparation-QA.md`

### 2025-11-10: MinION Pipeline Technical Investigation Report
**Comprehensive technical analysis of the 91-pathogen screening system**
- Created `docs/minion-pipeline-technical-report.md` - 80+ page technical documentation covering:
  - Sample preprocessing strategy (PCR-free, universal approach, shotgun metagenomics)
  - 6-phase pipeline architecture with AWS Lambda + EC2 implementation details
  - Variant detection capabilities (currently disabled, implementation roadmap provided)
  - Data retention and reanalysis feasibility (5-year S3 lifecycle, FAST5/FASTQ/BAM reprocessing)
  - Plasma volume requirements (5-10 mL standard, LOD 50-100 copies/mL)
- Detailed investigation of capture probe usage (none - uses CpG methylation-based host depletion instead)
- PERV-specific detection with immediate SNS alerting (subtypes A/B/C, recombinant detection)
- Cost-benefit analysis for variant calling Phase 7 implementation (+$50-100/sample)

### 2025-11-10: Documentation Structure Reorganization
**Organized root markdown files into structured directories**
- Created `docs/audits/` with 10 audit reports (9 audits + summary, 37 bugs fixed)
- Created `docs/bug-fixes/` with detailed bug fix documentation
- Created `docs/sprints/` with sprint completion reports
- Added comprehensive README.md to each directory with timelines and summaries
- Improved documentation discoverability and maintainability

### 2025-10-09: NGS vs Traditional Methods Cost Analysis
**Created comprehensive cost-benefit analysis document**
- `NGS全量解析vs従来法ハイブリッド戦略_コスト・手間分析.md` - Compares 3 strategies for PMDA 91-pathogen screening
- Pattern A (NGS all 91): ¥162,574/sample, 20h hands-on, 3-5 days turnaround
- Pattern B (Hybrid NGS+Traditional): ¥449,574/sample, 72h hands-on, 7-10 days (2.8× cost increase)
- Pattern C (Traditional only): ¥315,000/sample, regulatory non-compliant
- NGS-only approach demonstrates superior cost efficiency, faster turnaround, and full regulatory compliance

### 2025-10-08: MinION Protocol Appendices Complete
**Added 3 missing appendices to lab protocols**
- Appendix A: Reagents & Equipment List (¥21.5M initial investment)
- Appendix B: Time & Cost Estimates (¥127k/sample, 4-year ROI)
- Appendix C: Troubleshooting Guide (25 common issues)

### 2025-10-08: Architecture Documentation Update
**Updated all core documentation to reflect Lambda + EC2 custom AMI architecture (containerless, no Docker)**

For complete historical changes and development session logs, see:
- [Session History](../docs/claude-sessions/README.md)
- [Sprint Reports](../docs/sprints/README.md)
- [Audit Reports](../docs/audits/README.md)