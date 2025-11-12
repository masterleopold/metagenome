# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Project**: MinION Pathogen Screening Pipeline - PMDA-compliant xenotransplantation donor pig screening system for 91 designated pathogens.

## Development Commands

```bash
# Environment Setup
pip install -r requirements.txt
aws configure  # Set AWS_REGION=ap-northeast-1

# Testing
python -m pytest tests/                          # All tests
python -m pytest tests/ -k "not integration"     # Unit tests only
python -m pytest tests/test_pmda_compliance.py   # PMDA compliance
python -m pytest tests/test_pathogen_detection.py -v  # Specific test with verbose

# Pipeline Operations
./tools/workflow_cli.py start --run-id RUN-2024-001 --bucket minion-data --input-prefix runs/RUN-2024-001/fast5/
./tools/workflow_cli.py status --run-id RUN-2024-001 --watch
./tools/workflow_cli.py metrics --run-id RUN-2024-001

# Code Quality
black scripts/ lambda/ tools/ tests/    # Auto-format
flake8 scripts/ lambda/ tools/ tests/   # Linting
mypy scripts/                           # Type checking

# Documentation Portal (Next.js)
cd docs-portal && npm run dev    # http://localhost:3003
cd docs-portal && npm run build  # Production build
```

## Critical Context

### PERV Detection (HIGHEST PRIORITY)
- **PERV-A/B/C** subtypes must be detected and quantified
- Immediate SNS alert on ANY PERV detection
- Key scripts: `perv_typing.py`, `detect_recombinants.py`, `perv_phylogenetics.py`
- Must report in copies/mL plasma

### PMDA 91 Pathogens
- Full list: `md/厚労省異種移植指針_91病原体リスト.md`
- Requirements: PPA >95%, NPA >98%, R² >0.90
- Database: `/mnt/efs/databases/pmda/2024.1/`

### Sample Preparation Protocols (2025-11 Update)
**Recommended: Protocol 12 - Unified Workflow (covers all 91 pathogens)**
- **Workflows**: 2 universal (DNA + RNA) - replaces 3-4 complex variants
- **Time**: 13 hours hands-on (19% reduction)
- **Cost**: ¥157,000/sample (+¥5,000 over baseline)
- **Coverage**: 100% (91/91 pathogens)
- **LOD**: 100-500 copies/mL (screening-sufficient)
- **Key features**: Universal poly(A) selection, CpG methylation-based host depletion
- **Protocol**: `md/MinION_Protocol_12_統合サンプル調製プロトコル.md`

**Optional: Protocol 11 - High-Sensitivity Enhancement**
- Use ONLY when LOD <50 copies/mL required for specific viruses
- Target: Polyomavirus, Hantavirus, EEEV, Spumavirus
- **Protocol**: `md/MinION_Protocol_11_PMDA_4ウイルス高感度検出プロトコル.md`

**Conditional: Protocol 13 - Spumavirus-Specific Screening**
- Triggered only when Phase 1 detects retrovirus pol signatures (5-10% of samples)
- **Detection probability**: <0.1% (0 detections in 70 years)
- **Rationale**: PMDA precautionary testing despite zero historical detections
- **Method**: Nested PCR + phylogenetic analysis + PERV discrimination
- **Protocol**: `md/MinION_Protocol_13_スピューマウイルス専用検査プロトコル.md`
- **Context**: Section "科学的背景補足: ブタスピューマウイルスの実在性について" explains scientific background

## Architecture Overview

### Pipeline Orchestration Pattern
**Lambda-triggered EC2 pattern** (no containers):
1. S3 upload triggers Lambda orchestrator (`lambda/orchestration/pipeline_orchestrator.py`)
2. Lambda launches Step Functions workflow
3. Each phase runs on dedicated EC2 instance with custom AMI
4. EC2 instances use UserData scripts and auto-terminate
5. Phase completion triggers next Lambda → EC2 cycle

### Phase Structure (7-Phase Pipeline - Phase 0 NEW)
```
scripts/
├── phase0_sample_prep/     # NEW: Sample routing & workflow determination
│   ├── sample_router.py    # Determines DNA vs RNA extraction workflow
│   └── README.md           # Phase 0 usage documentation
├── phase1_basecalling/     # FAST5→FASTQ (GPU: g4dn.xlarge, Dorado)
├── phase2_qc/              # Quality metrics (t3.large, NanoPlot/PycoQC)
│                           # RNA integrity (RIN) for RNA viruses
├── phase3_host_removal/    # Host depletion (r5.4xlarge, Minimap2)
│                           # CpG methylation (DNA), rRNA depletion (RNA), poly(A) selection
├── phase4_pathogen/        # Multi-DB screening (4 parallel EC2, Kraken2/BLAST)
│   ├── perv_typing.py      # PERV-A/B/C subtype detection
│   ├── detect_recombinants.py
│   ├── perv_phylogenetics.py
│   ├── pmda_targeted_search.py
│   └── detect_pmda_4viruses.py  # NEW: Polyoma/Hantavirus/EEEV/Spumavirus
├── phase5_quantification/  # Abundance calc (t3.large, copies/mL)
└── phase6_reports/         # PMDA-compliant reports (t3.large)
```

### Lambda Functions Organization
```
lambda/
├── orchestration/          # Pipeline coordination
│   ├── pipeline_orchestrator.py  # Main S3 trigger handler
│   └── state_machine_handler.py  # Step Functions logic
├── monitoring/             # EC2 lifecycle and metrics
│   ├── check_phase_status.py
│   ├── collect_metrics.py
│   └── terminate_instances.py
├── api/                    # REST API endpoints
│   ├── get_workflow_status.py
│   └── list_workflows.py
└── shared/                 # Common utilities
    └── input_validation.py # Reusable validation logic
```

### Key Directories
```
tools/              # workflow_cli.py (start/status/metrics), monitoring_dashboard.py
tests/              # pytest suite (uses moto for AWS mocking)
infrastructure/     # Terraform IaC, CloudFormation templates
templates/          # Config YAML, Step Functions JSON, report HTML
docs/               # Technical docs (see TECHNICAL_DETAILS.md, TROUBLESHOOTING.md)
docs-portal/        # Next.js documentation site (port 3003)
md/                 # Japanese strategic documents and protocols
```

## Code Conventions

### Python Module Import Pattern
Phase scripts expect shared libraries in `/opt/minion/lib`:
```python
import sys
sys.path.append('/opt/minion/lib')
from workflow_manager import WorkflowManager
```

### Error Handling Pattern
Always validate files before opening:
```python
if not bam_file.exists():
    raise FileNotFoundError(f"BAM file not found: {bam_file}")
try:
    bam = pysam.AlignmentFile(str(bam_file), "rb")
except Exception as e:
    raise RuntimeError(f"Failed to open BAM file: {e}")
```

### Test Patterns
- Use `moto` for AWS service mocking
- Use absolute paths relative to test file:
```python
test_dir = Path(__file__).parent
repo_root = test_dir.parent
config_file = repo_root / 'templates' / 'config' / 'pmda_pathogens.json'
```

### PERV Detection Implementation
PERV markers hardcoded in `scripts/phase4_pathogen/perv_typing.py`:
```python
PERV_MARKERS = {
    'PERV-A': {'env_start': 5800, 'env_end': 7400,
               'specific_motifs': ['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC']},
    'PERV-B': {'env_start': 5800, 'env_end': 7400,
               'specific_motifs': ['ATGGCAACCACCGTAGC', 'TGGAAACCTGGAAAACC']},
    'PERV-C': {'env_start': 5800, 'env_end': 7400,
               'specific_motifs': ['ATGGCAGCCACCATAGG', 'TGGAGACCTGGAAGAAC']}
}
```

### Configuration Files
Main config: `templates/config/default_pipeline.yaml`
```yaml
pipeline:
  phases:
    basecalling:
      min_qscore: 9       # Q30 accuracy target
    qc:
      min_reads: 10000000 # 10M minimum
      min_q30: 0.85       # 85% Q30
    pathogen:
      confidence_threshold: 0.1
      pmda_alert_enabled: true
```
PMDA pathogen list: `templates/config/pmda_pathogens.json` (91 pathogens)

## Critical Warnings

⚠️ **PERV Detection** - Any detection triggers immediate alert
⚠️ **Data Security** - No patient data in git, use encrypted S3 only
⚠️ **PMDA Compliance** - All changes must maintain regulatory compliance
⚠️ **Host Depletion** - Must achieve >90% removal efficiency
⚠️ **Cost Threshold** - Alert if analysis exceeds $400

## Key Technical Details

### AWS Infrastructure
- **Region**: ap-northeast-1 (Tokyo)
- **Database**: RDS PostgreSQL (Aurora Serverless v2)
- **Shared Storage**: EFS for reference databases (Kraken2, BLAST, PERV)
- **Cost Optimization**: Spot Instances (70% savings), auto-terminating EC2
- **Monitoring**: CloudWatch + SNS alerts, Prometheus metrics

### Data Flow
S3 Upload (FAST5) → Lambda Orchestrator → Step Functions → EC2 Phase Execution → Results to S3 → SNS Notifications

### Important Paths
- Reference DBs: `/mnt/efs/databases/pmda/2024.1/`
- Shared libs: `/opt/minion/lib`
- Config templates: `templates/config/`
- Step Functions: `templates/stepfunctions/workflow_definition.json`

## Documentation

### Technical Documentation
- [Technical Details](docs/TECHNICAL_DETAILS.md) - Full architecture, phases, databases
- [MinION Pipeline Technical Report](docs/minion-pipeline-technical-report.md) - Comprehensive 80+ page technical investigation
- [Troubleshooting](docs/TROUBLESHOOTING.md) - Common issues and AWS debugging
- [Changelog](docs/CHANGELOG.md) - Recent updates
- [Development Workflow](docs/development/WORKFLOW_GUIDE.md) - Complete setup guide
- [API Documentation](docs/API_DOCUMENTATION.md) - REST API endpoints

### PMDA Compliance Documentation
- [PMDA Simplified Sample Prep Strategy](docs/PMDA_Simplified_Sample_Prep_Strategy.md) - Protocol 12 technical justification
- [PMDA Simplified Workflow Flowchart](docs/PMDA_Simplified_Workflow_Flowchart.md) - Visual workflows and decision trees
- [Session History](docs/claude-sessions/README.md) - Development session logs and decisions

## Recent Updates (2025-11-13)

### Protocol 13: Spumavirus Scientific Background Addition
- **Added comprehensive supplementary section** to Protocol 13 explaining the scientific reality that porcine spumavirus has never been detected
- **Key facts**: 0 detections in 70 years (1954-2025), 0 NCBI sequences, 0 peer-reviewed publications
- **PMDA inclusion rationale**: Precautionary principle based on phylogenetic proximity to foamy virus-positive species
- **Detection probability**: <0.1% (most Phase 1 triggers are PERV false positives)
- **Guidance for technicians**: Realistic expectations, PERV discrimination protocols, value of negative data
- **Impact scenarios**: What would happen if truly detected (Nature/Science-level discovery, immediate PMDA report)
- **Documentation**: Session log at `docs/claude-sessions/2025-11-13-protocol-13-spumavirus-scientific-background.md`

### Documentation Portal Updates
- **Fixed pathogen category counts**: Viruses (41), Bacteria (27), Parasites (19), Fungi (2), Special Management (5)
- **Added Protocol 12 section**: Unified workflow details, time/cost metrics, key features
- **Added Conditional Screening Strategy section**: Two-tier detection approach, universal vs conditional testing
- **Build verified**: Documentation portal builds successfully with no TypeScript errors