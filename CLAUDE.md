# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository implements a production-ready MinION pathogen screening pipeline for xenotransplantation donor pig screening, with comprehensive strategic planning documentation. The system detects 91 PMDA-designated pathogens using Oxford Nanopore MinION long-read sequencing and AWS cloud infrastructure.

**Architecture**: Serverless AWS pipeline with Lambda orchestration, EC2 compute, and Step Functions workflow management.

**Key Components**:
- Production Python codebase (`scripts/`, `lambda/`, `lib/`, `tools/`)
- Terraform infrastructure-as-code (`infrastructure/terraform/`)
- Strategic planning documents (`md/` - Japanese language protocols and reports)
- Testing framework (`tests/`)
- Documentation portal (`docs-portal/`)

**Key Documents:**
- `異種移植用ドナーブタの病原体検査体制構築計画書.md` - Detailed construction plan for pathogen testing system
- `専門家向け詳細報告書：異種移植用ドナーブタにおける病原体メタゲノム解析体制の構築.md` - Expert-level detailed report on system construction
- `NGS中心型病原体検査システム_最適化戦略v2.md` - NGS-centric optimization strategy (v2.0) with PCR-minimal approach
- `MinION_Protocol_00_目次とマスタードキュメント.md` - Master protocol index for MinION-based workflow (10 protocol documents)
- `MinION_AWS_クラウド解析パイプライン詳細設計.md` - AWS cloud-based analysis pipeline architecture for MinION data
- `MinION単独臨床試験システム_PMDA適合性評価.md` - PMDA compliance evaluation for MinION standalone system
- `MinION検証研究計画書_MiSeq外注版.md` - 12-month validation research plan comparing MinION in-house vs MiSeq outsourcing
- `NGSプラットフォーム完全比較_Nanopore含む.md` - Complete NGS platform comparison including Nanopore technology
- `臨床試験フェーズ向けNGSプラットフォーム比較_マルチプレックス最適化.md` - NGS platform comparison optimized for clinical trial phases with multiplexing strategies
- `厚労省異種移植指針_91病原体リスト.md` - MHLW xenotransplantation guidelines and 91 pathogen list

## Project Context

### Research Background
- **Target Animal**: Yucatan miniature pig (genetically modified with 69 genomic edits: 3KO, 7TG, 59 PERV inactivations)
- **Sample Type**: Plasma samples for cfDNA/cfRNA extraction
- **Regulatory Framework**: PMDA guidelines for xenotransplantation safety
- **Primary Goal**: Comprehensive detection of both known (91 PMDA-designated pathogens) and unknown pathogens

### Technical Approach
The project employs a three-phase workflow with multiple platform options:

1. **Lab Work Phase**:
   - Blood plasma separation and preservation
   - cfDNA/cfRNA extraction using magnetic bead methods (Zymo Quick-cfDNA/cfRNA™ Serum & Plasma Kit)
   - Host DNA depletion using endonuclease-based methods (NEBNext Microbiome DNA Enrichment Kit)

2. **Sequencing Phase** (Platform Options):
   - **Illumina Short-Read Platforms**: NovaSeq 6000, NextSeq 2000, or MiSeq (Paired-end 150bp)
   - **Oxford Nanopore MinION**: Long-read sequencing with Duplex mode for higher accuracy
     - DNA libraries: SQK-LSK114 kit
     - RNA libraries: Direct RNA or cDNA workflows
     - Real-time analysis capability
     - Lower initial investment (~$4,950 starter pack for MinION Mk1D includes device + 5 flow cells vs $50,000-$1M for Illumina)

3. **Bioinformatics Phase**:
   - Quality control (FastQC, MultiQC, PycoQC for MinION)
   - Host genome removal (BWA/Bowtie2 for short-read, Minimap2 for long-read)
   - Pathogen detection via multi-database search (RVDB, Kraken2, BLAST, Diamond)
   - Quantitative analysis and reporting
   - Optional: AWS cloud-based pipeline for scalable MinION data processing

## Phased Implementation Strategy

The NGS analysis workflow follows a 4-phase implementation approach over 5+ years:

### Phase 1 (0-1 year): Minimal In-house Implementation
- **In-house**: Nucleic acid extraction, quality evaluation
- **Outsourced**: Sequencing and analysis
- **Key Strategy**: Acquire complete analysis program code from contractor
- **Initial Investment**: ¥13M
- **Annual Cost**: ¥15.9M (6 analyses/year)

### Phase 2 (1-3 years): Analysis In-house Implementation
- **In-house**: Add bioinformatics analysis after receiving FASTQ data
- **Outsourced**: Sequencing only
- **Investment**: ¥30M (computing infrastructure)
- **Annual Cost**: ¥9.8M (12 analyses/year)
- **ROI**: 2.2 year payback period

### Phase 3 (3-5 years): Complete Sequencing In-house
- **In-house**: All processes including NGS sequencing
- **Investment**: ¥20M (NGS equipment)
- **Annual Cost**: ¥11M (24 analyses/year)
- **Trigger**: 24+ analyses per year

### Phase 4 (5+ years): Full R&D Capability
- **Focus**: Technology innovation, AI/ML integration
- **Investment**: ¥20M (R&D facilities)
- **Goal**: Competitive advantage through proprietary methods

## Key Technical Components

### Database Architecture
Multi-layered database strategy for comprehensive pathogen coverage:

- **RVDB v30.0**: Manually curated viral sequences (primary virus database)
- **Kraken2 Standard DB**: RefSeq archaea, bacteria, viruses, plasmids
- **NCBI RefSeq**: High-quality reference genomes for bacteria/viruses
- **GISAID**: Pandemic-relevant pathogen sequences
- **PATRIC**: Zoonotic bacteria genomics
- **PMDA Custom DB**: 91 designated pathogen sequences specific to Japanese regulations
- **Host Genome DB**: Sus scrofa + modified pig genome sequences

### Bioinformatics Pipeline

```
Raw FASTQ → QC (FastQC/Trimmomatic) → Host Removal (BWA/Bowtie2) →
Parallel Pathogen Detection:
  ├─ Kraken2/Bracken (rapid screening)
  ├─ RVDB/BLAST (virus-specific)
  ├─ Diamond (protein-level)
  └─ De novo assembly (unknown pathogens)
→ Result Integration → Quantification → Report Generation
```

### Quality Control Framework

**Analytical Validation Requirements:**
- Limit of Detection (LOD) determination using spiked samples
- Reproducibility across operators, instruments, and dates
- Specificity verification against host genome cross-reactivity
- Data quality metrics (Q30 >85%, >20M reads/sample, host depletion >95%)

**ALCOA+ Principles Implementation:**
- Attributable, Legible, Contemporaneous, Original, Accurate
- Complete, Consistent, Enduring, Available

## Regulatory Considerations

### PMDA Guidelines Compliance
- Detection of both "known and unknown infectious diseases"
- Donor animals from closed, documented environments
- Screening during minimum 3-week quarantine period
- Detection of pathogens during latency periods

### Validation Strategy
- Reference FDA guidance for NGS-based IVDs
- BloodPAC analytical validation protocols for cfDNA/cfRNA assays
- ISO 15189 medical laboratory accreditation
- GLP (Good Laboratory Practice) compliance

## Important Technical Specifications

### Sample Processing
- **Sample Volume**: 5-10mL plasma with EDTA anticoagulant
- **Storage**: -80°C for long-term preservation
- **Extraction Kit**: Zymo Quick-cfDNA/cfRNA™ Serum & Plasma Kit (or equivalent)
- **Expected Yield**: 1-100ng/mL cfNA (avg ~30ng/mL)
- **Host Depletion**: NEBNext Microbiome DNA Enrichment Kit (90-99% removal)

### Sequencing Parameters
- **Read Type**: Paired-end 150bp
- **Coverage Target**: Minimum 10M reads/sample
- **Sequencing Depth**: 50× average
- **Quality Threshold**: Q30 >85%
- **PhiX Control**: 1% spike-in

### Computational Requirements
- **CPU**: 64 cores × 8 nodes
- **Memory**: 1TB RAM per node
- **Storage**: 500TB NAS + 100TB SSD (RAID6)
- **Workflow Management**: Nextflow or Snakemake
- **Containerization**: Docker/Singularity

## Cost Structure Summary

### Total Investment by Phase
- Phase 1: ¥13M (nucleic acid extraction equipment, analysis program development)
- Phase 2: ¥43M cumulative (computing infrastructure, bioinformatics talent)
- Phase 3: ¥63M cumulative (NGS equipment, sequencing facilities)
- Phase 4: ¥83M cumulative (R&D facilities, research development)

### Annual Operating Costs (varies by phase and analysis frequency)
- Phase 1: ¥9.9M-47.6M (depending on 6-24 analyses/year, heavily outsourced)
- Phase 2: ¥9.8M-19.2M (analysis in-house, sequencing outsourced)
- Phase 3: ¥11M-18M (fully in-house)
- Phase 4: ¥14M-22M (includes R&D operations)

## Decision Criteria for Phase Transitions

### Phase 1 → Phase 2
- Annual analysis frequency: 12+ times/year
- Cumulative cost reduction: ¥5M+
- Technical proficiency: 80%+
- Quality metrics achievement: 95%+

### Phase 2 → Phase 3
- Annual analysis frequency: 24+ times/year
- Cumulative cost reduction: ¥20M+
- Technical proficiency: 90%+
- Quality metrics achievement: 98%+

### Phase 3 → Phase 4
- Annual analysis frequency: 36+ times/year
- Business expansion requirements
- Competitive advantage needs
- Technical proficiency: 95%+

## Risk Management

### High-Priority Risks
1. **Technology Learning Delays** (Medium probability, High impact)
   - Mitigation: External expert support, comprehensive training programs

2. **Equipment Failure** (Low probability, High impact)
   - Mitigation: Maintenance contracts, backup outsourcing arrangements

3. **Talent Attrition** (Medium probability, Medium impact)
   - Mitigation: Documentation, cross-training multiple staff members

4. **Regulatory Changes** (Medium probability, Medium impact)
   - Mitigation: Continuous regulatory monitoring, flexible response systems

## Working with This Repository

### When Adding New Documentation
- Use Japanese for primary documents (matches existing content)
- Include technical diagrams in Mermaid format where appropriate
- Reference PMDA guidelines and international standards (FDA, ISO)
- Document cost implications and ROI analysis

### When Updating Strategic Plans
- Maintain phase-based structure (Phases 1-4)
- Update cost projections with current market prices
- Include timeline adjustments based on project progress
- Cross-reference with regulatory requirement changes

### When Developing Technical Protocols
- Follow the three-phase workflow structure (Lab → Sequencing → Bioinformatics)
- Include quality control checkpoints at each stage
- Reference specific equipment models and software versions
- Document ALCOA+ compliance measures

## Key Technical Terms and Abbreviations

**Biological/Medical:**
- **cfDNA/cfRNA**: Cell-free DNA/RNA circulating in plasma
- **PERV**: Porcine Endogenous Retrovirus
- **3KO**: 3 Knockout genes (GGTA1, CMAH, B4GALNT2/B4GALNT2L)
- **7TG**: 7 Transgene (CD46, CD55, THBD, PROCR, CD47, TNFAIP3, HMOX1)

**Sequencing/Analysis:**
- **NGS**: Next Generation Sequencing
- **MinION**: Oxford Nanopore's portable sequencer
- **Duplex mode**: Nanopore sequencing mode that sequences both DNA strands for higher accuracy
- **FAST5**: Raw signal data format from Nanopore sequencers
- **FASTQ**: Text-based format for nucleotide sequences with quality scores
- **Basecalling**: Converting raw electrical signals to nucleotide sequences
- **TPM/RPM**: Transcripts/Reads Per Million (normalization methods)
- **UMI**: Unique Molecular Identifier for absolute quantification

**Databases:**
- **RVDB**: Reference Viral Database
- **NCBI RefSeq**: NCBI Reference Sequence Database
- **PMDA Custom DB**: 91 PMDA-designated pathogen sequences

**Regulatory/Quality:**
- **PMDA**: Pharmaceuticals and Medical Devices Agency (Japan)
- **MHLW**: Ministry of Health, Labour and Welfare (Japan)
- **LOD**: Limit of Detection
- **ALCOA+**: Data integrity principles (Attributable, Legible, Contemporaneous, Original, Accurate, Complete, Consistent, Enduring, Available)
- **GLP**: Good Laboratory Practice
- **IVD**: In Vitro Diagnostic

**Cloud/Infrastructure:**
- **AWS**: Amazon Web Services
- **S3**: Simple Storage Service (AWS object storage)
- **EC2**: Elastic Compute Cloud (AWS virtual servers)
- **TAT**: Turnaround Time

## External Resources

### Regulatory Guidance
- PMDA xenotransplantation guidelines (referenced throughout documents)
- FDA NGS-based IVD guidance
- BloodPAC analytical validation protocols
- ISO 15189 medical laboratory standards

### Technical Databases
- RVDB: https://rvdb.dbi.udel.edu/
- NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
- GISAID: https://www.gisaid.org/
- PATRIC: https://www.patricbrc.org/

### Key Technologies

**Sequencing Platforms:**
- Illumina: NovaSeq 6000, NextSeq 2000, MiSeq (short-read)
- Oxford Nanopore: MinION Mk1D, GridION (long-read)

**Workflow Management:**
- Nextflow, Snakemake (local pipelines)
- AWS Step Functions, AWS Batch (cloud-based workflows)

**Analysis Tools:**
- Classification: Kraken2, Bracken, KrakenUniq
- Alignment: BWA, Bowtie2 (short-read), Minimap2 (long-read)
- Assembly: MEGAHIT, metaSPAdes
- Annotation: BLAST, Diamond
- Basecalling: Guppy, Dorado (for MinION)

**QC Tools:**
- FastQC, MultiQC, Trimmomatic (short-read)
- PycoQC, NanoPlot (long-read/MinION)

**Cloud Infrastructure:**
- AWS S3 (data storage), EC2 (compute), RDS (databases)
- Docker/Singularity for containerization

## Document Organization

The repository contains three main protocol series:

### 1. MinION Protocol Series (10 documents)
Detailed wet-lab protocols for MinION-based pathogen detection:
- `MinION_Protocol_00_目次とマスタードキュメント.md` - Master index
- `MinION_Protocol_01_*.md` through `MinION_Protocol_10_*.md` - Step-by-step protocols
- Covers: Sample collection → Extraction → Library prep → Sequencing → QC
- Total lab time: ~16 hours + 6-48 hours sequencing

### 2. Strategic Planning Documents
High-level strategy and system design:
- NGS platform comparisons and optimization strategies (including multiplexing for clinical trials)
- Cost-benefit analysis for different implementation phases
- PMDA compliance evaluation
- AWS cloud architecture designs
- Validation research plan (MinION vs MiSeq outsourcing comparison)

### 3. Regulatory and Reference Documents
- PMDA/MHLW guidelines and pathogen lists
- Quality management frameworks
- Validation criteria (PPA >95%, NPA >98%, R² >0.90)

## Validation Research Framework

A critical component of the project is the 12-month validation study (`MinION検証研究計画書_MiSeq外注版.md`):

**Research Design:** Prospective comparative validation (MinION in-house vs MiSeq outsourcing)

**Three-Phase Approach:**
- **Phase 1 (Months 1-7):** Parallel validation - MinION internal analysis vs MiSeq outsourced analysis
- **Phase 2 (Months 8-10):** Independent verification - Optimized MinION protocol vs blinded MiSeq comparison
- **Phase 3 (Months 11-12):** Complete internalization - MinION standalone operation with quarterly MiSeq QC

**Exit Criteria for MiSeq Outsourcing Termination:**
1. Species-level identification concordance >95%
2. Quantification correlation R² >0.90
3. Positive Percent Agreement (PPA) >95%
4. Negative Percent Agreement (NPA) >98%
5. All criteria met in MinION independent verification

## Development Workflow

### Common Commands

**Setup and Installation:**
```bash
# Install dependencies
pip install -r requirements.txt

# Set up AWS credentials
aws configure
export AWS_REGION=ap-northeast-1
export ENVIRONMENT=production

# Deploy infrastructure
cd infrastructure/terraform
terraform init
terraform plan
terraform apply

# Set up reference databases
./tools/database_setup.sh --all

# Deploy pipeline components
./tools/deployment_script.sh deploy
```

**Running the Pipeline:**
```bash
# Start workflow via CLI
./tools/workflow_cli.py start \
  --run-id RUN-2024-001 \
  --bucket minion-data \
  --input-prefix runs/RUN-2024-001/fast5/

# Monitor workflow
./tools/workflow_cli.py status --run-id RUN-2024-001 --watch

# View metrics
./tools/workflow_cli.py metrics --run-id RUN-2024-001

# Launch monitoring dashboard
streamlit run tools/monitoring_dashboard.py
```

**Testing:**
```bash
# Run all tests
python -m pytest tests/

# Run unit tests only
python -m pytest tests/ -k "not integration"

# Run integration tests
python tests/test_pipeline_integration.py

# Run PMDA compliance tests
python tests/test_pmda_compliance.py

# Run specific phase tests
python tests/test_basecalling_module.py
python tests/test_pathogen_detection.py

# Run with coverage
python -m pytest tests/ --cov=lib --cov=scripts --cov-report=term-missing
```

**Configuration Management:**
```bash
# Create custom config
./tools/workflow_cli.py config --create-default

# Validate configuration
./tools/workflow_cli.py config --validate custom.yaml

# Default config location
vi templates/config/default_pipeline.yaml
```

**Development Tools:**
```bash
# Code formatting
black scripts/ lambda/ lib/ tools/ tests/

# Linting
flake8 scripts/ lambda/ lib/ tools/ tests/

# Type checking
mypy lib/ scripts/

# Pre-commit hooks
pre-commit install
pre-commit run --all-files
```

### Code Architecture

**High-Level System Flow:**
```
MinION Sequencer → S3 Upload → Lambda Trigger → Step Functions State Machine
                                                          ↓
                                            [6 Sequential EC2 Phases]
                                                          ↓
                                              Reports + Notifications
```

**Phase Organization (scripts/):**
The pipeline executes in 6 sequential phases, each with dedicated scripts:

1. **Phase 1: Basecalling** (`scripts/phase1_basecalling/`)
   - Converts FAST5/POD5 signal data to FASTQ sequences
   - Uses GPU-accelerated Dorado with duplex mode for Q30 accuracy
   - Key: `generate_summary_from_fastq.py` - Post-basecalling metrics

2. **Phase 2: QC** (`scripts/phase2_qc/`)
   - Quality assessment using NanoPlot/PycoQC
   - Key: `qc_check.py` - Validates quality thresholds

3. **Phase 3: Host Removal** (`scripts/phase3_host_removal/`)
   - Aligns reads to Sus scrofa genome using Minimap2
   - Key: `calculate_depletion_rate.py` - Host depletion metrics

4. **Phase 4: Pathogen Detection** (`scripts/phase4_pathogen/`)
   - Multi-database pathogen screening
   - **Critical PERV scripts:**
     - `perv_typing.py` - PERV-A/B/C subtype identification
     - `detect_recombinants.py` - PERV recombination detection
     - `perv_phylogenetics.py` - Phylogenetic analysis

5. **Phase 5: Quantification** (`scripts/phase5_quantification/`)
   - `kraken_quantify.py` - Kraken2-based abundance
   - `blast_quantify.py` - BLAST-based quantification
   - `spike_in_normalization.py` - PhiX174 spike-in normalization
   - `absolute_copy_number.py` - Copies/mL calculation

6. **Phase 6: Reporting** (`scripts/phase6_reports/`)
   - `generate_pmda_report.py` - PDF/JSON/HTML reports
   - `generate_pmda_checklist.py` - 91-pathogen compliance checklist

**Lambda Functions (lambda/):**
- `orchestration/pipeline_orchestrator.py` - Main orchestration logic
- `orchestration/state_machine_handler.py` - Step Functions integration
- `phases/` - Per-phase Lambda handlers that launch EC2 instances
- `api/` - REST API endpoints for workflow management
- `monitoring/` - CloudWatch metrics and alerting

**Core Libraries (lib/):**
Shared Python modules used across scripts and Lambda functions:
- `workflow_manager.py` - Workflow lifecycle management
- `config_manager.py` - Configuration validation and loading
- `s3_utils.py` - S3 data operations
- `ec2_manager.py` - EC2 instance lifecycle (launch, monitor, terminate)
- `database_client.py` - RDS Aurora interactions
- `reference_manager.py` - Reference database management
- `monitoring_client.py` - CloudWatch metrics
- `notification_client.py` - SNS/email alerting
- `report_generator.py` - Multi-format report generation

**Infrastructure (infrastructure/terraform/):**
Modular Terraform configuration:
- `main.tf` - Provider and backend configuration
- `variables.tf` - Input variables
- `outputs.tf` - Export values
- `s3.tf` - Data buckets and lifecycle policies
- `lambda.tf` - Lambda functions and layers
- `ec2_ami.tf` - Custom AMIs with pre-installed tools
- `efs.tf` - Reference database storage
- `rds.tf` - Aurora Serverless metadata database
- `iam.tf` - IAM roles and policies
- `cloudwatch.tf` - Monitoring dashboards and alarms
- `sns.tf` - Notification topics
- `eventbridge.tf` - Event-driven triggers

**Tools (tools/):**
- `workflow_cli.py` - CLI for workflow management (start, status, metrics)
- `monitoring_dashboard.py` - Streamlit web dashboard
- `database_setup.sh` - Download and configure reference databases
- `deployment_script.sh` - Deploy pipeline components to AWS

**Configuration System:**
- Default config: `templates/config/default_pipeline.yaml`
- Custom configs override defaults
- Validation via `ConfigManager` class
- Per-phase settings (basecalling, QC, pathogen detection, etc.)
- Resource allocation (EC2 instance types, timeouts)
- PMDA compliance settings (critical pathogens, alert thresholds)

**Data Flow:**
1. Raw data uploaded to S3 bucket (`s3://minion-data/{run_id}/fast5/`)
2. Lambda orchestrator triggered, starts Step Functions
3. Each phase:
   - Lambda launches EC2 instance with phase-specific AMI
   - EC2 pulls input from S3, executes scripts, pushes output to S3
   - Updates workflow status in DynamoDB
   - Auto-terminates on completion
4. Final reports written to `s3://minion-data/{run_id}/reports/`
5. SNS notifications sent for critical events (PERV detection, failures)

**Testing Strategy:**
- Unit tests for individual functions/modules
- Integration tests for end-to-end pipeline phases
- PMDA compliance tests verify 91-pathogen coverage
- Mock AWS services using `moto` library
- Test data in `tests/data/` (synthetic samples, no real patient data)

### Working with PERV Detection (Critical)

PERV (Porcine Endogenous Retrovirus) detection is the **highest priority** for PMDA compliance:

**Key Scripts:**
- `scripts/phase4_pathogen/perv_typing.py` - Identifies PERV-A, PERV-B, PERV-C subtypes using envelope gene markers
- `scripts/phase4_pathogen/detect_recombinants.py` - Detects PERV-A/C recombinants
- `scripts/phase4_pathogen/perv_phylogenetics.py` - Phylogenetic placement

**PERV Markers:**
```python
PERV_MARKERS = {
    'PERV-A': {
        'env_start': 5800,
        'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC']
    },
    'PERV-B': {...},
    'PERV-C': {...}
}
```

**Detection Requirements:**
- Immediate SNS alert on any PERV detection
- Separate reporting section for PERV-specific results
- Quantification in copies/mL plasma
- Phylogenetic analysis for strain identification

### Database Management

**Reference Databases (stored on EFS):**
- **Kraken2 Standard DB**: RefSeq bacteria, viruses, archaea (~50GB)
- **RVDB v30.0**: Curated viral sequences (~2GB)
- **PMDA Custom DB**: 91 designated pathogen sequences
- **Sus scrofa genome**: Host reference for depletion (~3GB)

**Setup:**
```bash
./tools/database_setup.sh --all              # Download all databases
./tools/database_setup.sh --kraken2          # Single database
./tools/database_setup.sh --local            # Local testing setup
```

**Database Locations on EFS:**
```
/mnt/efs/databases/
├── kraken2/standard/
├── rvdb/v30.0/
├── pmda/2024.1/
└── host_genomes/sus_scrofa_11.1/
```

**Version Management:**
- Kraken2: `standard_20231009` (auto_update: false)
- RVDB: `v30.0` (auto_update: false)
- PMDA: `2024.1` (auto_update: true) - Critical for regulatory compliance

### API Integration

**Base URL:** `https://api.minion-pipeline.com/{environment}/`

**Authentication:** All requests require `x-api-key` header

**Key Endpoints:**
- `POST /workflows` - Start new workflow
- `GET /workflows/{workflow_id}` - Get workflow status
- `GET /workflows/{workflow_id}/metrics` - Get performance metrics
- `GET /workflows/{workflow_id}/results` - Get analysis results
- `DELETE /workflows/{workflow_id}` - Cancel workflow

**Example API Call:**
```bash
curl -X POST https://api.minion-pipeline.com/production/workflows \
  -H "x-api-key: YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "run_id": "RUN-2024-001",
    "bucket": "minion-data",
    "input_prefix": "runs/RUN-2024-001/fast5/"
  }'
```

**Response:**
```json
{
  "run_id": "RUN-2024-001",
  "workflow_id": "123456",
  "execution_arn": "arn:aws:states:...",
  "status": "STARTED",
  "estimated_duration_hours": 8.5
}
```

Full API documentation: `docs/API_DOCUMENTATION.md`

### Monitoring and Debugging

**CloudWatch Dashboard:**
```
https://console.aws.amazon.com/cloudwatch/home?region=ap-northeast-1#dashboards:name=minion-pipeline-production
```

**Key Metrics Tracked:**
- Workflow execution time (per phase and total)
- Phase completion rates and success/failure ratios
- Pathogen detection counts (especially PERV detections)
- EC2 resource utilization (CPU, memory, GPU)
- Cost per analysis run
- Error rates and retry counts

**Log Access:**
```bash
# Lambda logs
aws logs tail /aws/lambda/minion-pipeline-production --follow

# EC2 phase logs (via SSM)
aws ssm start-session --target INSTANCE_ID

# Step Functions execution history
aws stepfunctions get-execution-history --execution-arn ARN
```

**Alert Configuration:**
- **CRITICAL**: PERV detection (immediate SNS notification)
- **CRITICAL**: Workflow failure (immediate notification)
- **WARNING**: QC threshold violations
- **WARNING**: High host contamination (>95%)
- **WARNING**: Cost threshold exceeded ($400 per run)

**Troubleshooting Common Issues:**

1. **Basecalling Timeout:**
   - Check GPU availability on EC2 instance
   - Verify CUDA drivers installed
   - Increase `max_runtime_hours` in config
   - Consider larger instance type (g4dn.2xlarge)

2. **High Host Contamination:**
   - Review nucleic acid extraction protocol
   - Check host depletion efficiency (<90%)
   - Verify Minimap2 alignment parameters
   - Validate Sus scrofa reference genome version

3. **Low/No Pathogen Detection:**
   - Verify database integrity and versions
   - Check confidence thresholds (default 0.1)
   - Review sample quality metrics
   - Ensure sufficient read depth (>10M reads)

4. **EC2 Launch Failures:**
   - Check spot instance availability
   - Fallback to on-demand if spot unavailable
   - Verify IAM role permissions
   - Check VPC subnet capacity

## Notes for Future Development

When working with this codebase:

1. **Maintain consistency** with the phased implementation approach (Phases 1-4 over 5+ years)
2. **Ensure PMDA compliance** - all technical specifications must align with PMDA regulatory requirements
3. **Update cost projections** to reflect current market conditions (especially for reagents and cloud services)
4. **Consider scalability** when proposing technical solutions (current focus: 24 samples/year for Phase I clinical trial)
5. **Document validation strategies** for any new methods - must include LOD, reproducibility, specificity, PPA/NPA metrics
6. **Maintain bilingual capability** (Japanese primary, English technical terms)
7. **Cross-reference documents** - when updating protocols, check if strategic documents need corresponding updates
8. **Version control** - major protocol changes should be tracked with version numbers (e.g., v2.0 for strategy shift from PCR-hybrid to NGS-only)
9. **Reference the validation framework** - all protocol modifications should consider impact on the MinION vs MiSeq validation study

### Coding Standards

**Python Code Requirements:**
- Use Google-style docstrings for all public functions and classes
- Type hints required for function signatures
- Follow PEP 8 style guide (enforced by black + flake8)
- Maximum line length: 88 characters (black default)
- Use `logging` module, never `print()` statements in production code

**Example Function:**
```python
def detect_pathogens(
    self,
    fastq_path: str,
    pmda_only: bool = False
) -> Dict[str, float]:
    """
    Detect pathogens in FASTQ file.

    Args:
        fastq_path: Path to input FASTQ file
        pmda_only: If True, only detect PMDA pathogens

    Returns:
        Dictionary mapping pathogen codes to confidence scores

    Raises:
        FileNotFoundError: If FASTQ file does not exist
        ValueError: If FASTQ file is malformed

    Example:
        >>> detector = PathogenDetector(config)
        >>> results = detector.detect_pathogens("sample.fastq")
        >>> print(results["PERV-A"])
        0.98
    """
    logger.info(f"Starting pathogen detection on {fastq_path}")
    # Implementation
```

**Shell Scripts:**
- Use `#!/usr/bin/env bash` shebang
- Set strict mode: `set -euo pipefail`
- Include header with description, author, version
- Implement error handling with trap
- Use readonly for constants
- Quote all variable expansions

**Terraform:**
- Use descriptive resource names: `${var.project_name}-{resource}-${var.environment}`
- Include comprehensive tags (Name, Environment, Purpose, ManagedBy, CostCenter, Owner)
- Enable monitoring and IMDSv2 for EC2 instances
- Document dependencies in comments
- Use `lifecycle` blocks to prevent unwanted resource replacement

**Git Commit Messages:**
Follow conventional commits specification:
```bash
feat(pathogen): add PERV-C detection algorithm
fix(basecalling): correct duplex mode parameters
docs(api): update endpoint documentation
perf(analysis): optimize Kraken2 memory usage
test(compliance): add PMDA validation tests
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`, `revert`

**Pull Request Requirements:**
- Minimum 2 approvals required
- PMDA-related changes need compliance team review
- All CI/CD checks must pass
- Code coverage must not decrease
- No new security vulnerabilities
- Signed commits (GPG) required

For complete contribution guidelines, see `CONTRIBUTING.md`
