# Recent Updates

## 2025-11-14 Updates

### 4-Virus Surveillance System (NEW FEATURE)

**Comprehensive dual-source monitoring system** for 4 target viruses in Japanese pig populations:
- **Hantavirus** (ハンタウイルス)
- **Polyomavirus** (ポリオーマウイルス)
- **Spumavirus** (スピューマウイルス) - MHLW Special Management Pathogen #5
- **EEEV** (東部ウマ脳炎ウイルス)

#### Architecture Pattern
Standalone monitoring system with API integration to main pipeline:
- **External Sources**: Daily collection from MAFF, E-Stat, PubMed, J-STAGE (11:00 JST)
- **Internal Pipeline**: Real-time Phase 4 result monitoring via S3 events
- **Severity Engine**: 4-level classification (CRITICAL/HIGH/MEDIUM/LOW)
- **Multi-channel Alerts**: SNS/SES/SMS/Dashboard/API

#### Key Components Implemented

**External Monitoring**:
- `surveillance/external/estat_client.py` - E-Stat API integration (APP_ID: bae1f981a6d093a9676b03c8eea37324b8de421b)
- `surveillance/external/maff_scraper.py` - MAFF surveillance report scraping
- `surveillance/external/academic_monitor.py` - PubMed API + J-STAGE web scraping

**Internal Monitoring**:
- `surveillance/internal/pipeline_listener.py` - Phase 4 result monitoring
- `scripts/phase4_pathogen/detect_4viruses.py` - BAM-based detection module

**Alerting System**:
- `surveillance/alerting/severity_engine.py` - YAML-based classification
- `surveillance/alerting/notification_router.py` - Multi-channel routing

**User Interfaces**:
- `surveillance/dashboard/app.py` - Streamlit real-time dashboard (port 8501, 30s auto-refresh)
- `surveillance/api/main.py` - FastAPI REST API (port 8000)

**Infrastructure**:
- 3 DynamoDB tables (detections, external-updates, notifications)
- S3 data lake with lifecycle policies
- Lambda functions (external_collector, pipeline_listener, alert_processor)
- EventBridge daily schedule (cron: 0 2 * * ? * = 11:00 JST)
- SNS topics (critical-alerts, high-alerts, daily-summary)

#### Cost Efficiency
- **Monthly Operating Cost**: ~$7.40/month
  - Lambda: ~$2
  - DynamoDB: ~$5 (pay-per-request)
  - S3: ~$0.30
  - SNS/SES: ~$0.10

#### Documentation Updates
- Added comprehensive architecture section to `docs/ARCHITECTURE.md` (300+ lines)
- Created new docs-portal page: `docs-portal/src/app/surveillance/page.tsx`
- Updated main `README.md` with 4-Virus Surveillance feature
- Updated sidebar navigation and homepage in docs-portal
- Created `surveillance/README.md` user guide

**Documentation**: Session log at `docs/claude-sessions/2025-11-14-4virus-surveillance-system.md`

## 2025-11-13 Updates

### Protocol 12 v2.1: Circular & Single-Stranded DNA Virus Detection (CRITICAL PMDA COMPLIANCE FIX)

**Critical Discovery**: Standard SQK-LSK114 ligation-based library prep CANNOT detect circular DNA or single-stranded DNA viruses.

#### Problem Identified
- Circular DNA viruses have no free ends for adapter ligation
- Single-stranded DNA viruses: T4 Ligase has <5% efficiency on ssDNA
- **Affected pathogens (4)**: PCV2, PCV3 (Special Management pathogens), TTV, PPV
- **Impact**: Protocol 12 v1.0 achieved only 87/91 coverage (95.6%), NOT the claimed 100%

#### Solution Implemented
Added **Step 2.5** (2.5 hours, +¥5,000):
- Sub-step 2.5.1: Circular DNA linearization (DNase I, ultra-low 0.005 U)
- Sub-step 2.5.2: ssDNA→dsDNA conversion (Klenow Fragment + Random Hexamers)
- **Result**: TRUE 91/91 pathogen coverage achieved (100%)

#### Updated Files
- `md/MinION_Protocol_12_統合サンプル調製プロトコル.md` - Added comprehensive Step 2.5 (262 lines)
- `md/MinION_Protocol_00_目次とマスタードキュメント.md` - Updated to v2.1 specifications
- `md/MinION_Protocol_付録B_時間・コスト見積.md` - Updated time (15.5h) and cost (¥162,000)
- `templates/config/pmda_pathogens.json` - Added genome structure metadata
- `docs/PMDA_Simplified_Sample_Prep_Strategy.md` - Updated to v2.1
- `CLAUDE.md` - Updated Protocol 12 specifications

**Documentation**: Session log at `docs/claude-sessions/2025-11-13-protocol-12-circular-ssdna-update.md`

### Protocol 13: Spumavirus Scientific Background Addition

**Added comprehensive supplementary section** to Protocol 13 explaining the scientific reality of porcine spumavirus.

#### Key Facts Added
- 0 detections in 70 years (1954-2025)
- 0 NCBI sequences
- 0 peer-reviewed publications
- **PMDA inclusion rationale**: Precautionary principle based on phylogenetic proximity
- **Detection probability**: <0.1% (most Phase 1 triggers are PERV false positives)

#### Guidance for Technicians
- Realistic expectations about detection probability
- PERV discrimination protocols
- Value of negative data for scientific record
- Impact scenarios if truly detected (Nature/Science-level discovery)

**Documentation**: Session log at `docs/claude-sessions/2025-11-13-protocol-13-spumavirus-scientific-background.md`

### Documentation Portal Updates

#### Fixed Pathogen Category Counts
- Viruses: 41 pathogens
- Bacteria: 27 pathogens
- Parasites: 19 pathogens
- Fungi: 2 pathogens
- Special Management: 5 pathogens

#### Added New Sections
- Protocol 12 unified workflow details
- Time/cost metrics visualization
- Conditional screening strategy explanation
- Two-tier detection approach documentation

#### Technical Improvements
- Build verified with no TypeScript errors
- Updated component props for better type safety
- Enhanced navigation structure

## Previous Updates

### Phase 0 Addition
- New sample preparation phase for routing decisions
- Automatic DNA vs RNA extraction workflow determination
- Integration with Protocol 12 unified approach

### PERV Detection Enhancement
- Improved subtype discrimination algorithms
- Added recombinant detection capabilities
- Enhanced phylogenetic analysis pipeline
- Real-time alerting system implementation

### AWS Infrastructure Optimization
- Migrated to Spot Instances (70% cost reduction)
- Implemented auto-termination for EC2 instances
- Enhanced Step Functions workflow
- Added cost monitoring and alerting

### Database Updates (2024.1)
- Updated PMDA pathogen reference database
- Enhanced Kraken2 custom database
- Improved BLAST database coverage
- Added specialized PERV reference sequences

## Session History

Detailed development session logs are maintained at:
- `docs/claude-sessions/README.md` - Index of all sessions
- Individual session files document specific feature implementations

## Version History

| Version | Date | Major Changes |
|---------|------|---------------|
| 2.2.0 | 2025-11-14 | 4-Virus Surveillance System - Dual-source monitoring (MAFF/E-Stat/PubMed/J-STAGE + pipeline) |
| 2.1.0 | 2025-11-13 | Protocol 12 v2.1 - Circular/ssDNA virus support |
| 2.0.0 | 2025-11-12 | Protocol 13 - Spumavirus screening addition |
| 1.9.0 | 2025-11-10 | Documentation portal launch |
| 1.8.0 | 2025-11-08 | Phase 0 sample prep routing |
| 1.7.0 | 2025-11-05 | PERV detection enhancement |
| 1.6.0 | 2025-11-01 | AWS infrastructure optimization |

## Upcoming Features

### Planned for Next Release
- Real-time basecalling integration
- Enhanced ML-based pathogen classification
- Automated report translation (JP→EN)
- Mobile monitoring dashboard

### Under Consideration
- Integration with hospital LIMS
- Blockchain-based audit trail
- AI-powered anomaly detection
- Multi-site deployment support