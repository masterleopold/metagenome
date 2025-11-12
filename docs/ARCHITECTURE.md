# Architecture Documentation

## Pipeline Orchestration Pattern

**Lambda-triggered EC2 pattern** (no containers):
1. S3 upload triggers Lambda orchestrator (`lambda/orchestration/pipeline_orchestrator.py`)
2. Lambda launches Step Functions workflow
3. Each phase runs on dedicated EC2 instance with custom AMI
4. EC2 instances use UserData scripts and auto-terminate
5. Phase completion triggers next Lambda → EC2 cycle

## Phase Structure (7-Phase Pipeline)

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

## Lambda Functions Organization

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

## Key Directories

```
tools/              # workflow_cli.py (start/status/metrics), monitoring_dashboard.py
tests/              # pytest suite (uses moto for AWS mocking)
infrastructure/     # Terraform IaC, CloudFormation templates
templates/          # Config YAML, Step Functions JSON, report HTML
docs/               # Technical docs (see TECHNICAL_DETAILS.md, TROUBLESHOOTING.md)
docs-portal/        # Next.js documentation site (port 3003)
md/                 # Japanese strategic documents and protocols
```

## AWS Infrastructure

- **Region**: ap-northeast-1 (Tokyo)
- **Database**: RDS PostgreSQL (Aurora Serverless v2)
- **Shared Storage**: EFS for reference databases (Kraken2, BLAST, PERV)
- **Cost Optimization**: Spot Instances (70% savings), auto-terminating EC2
- **Monitoring**: CloudWatch + SNS alerts, Prometheus metrics

## Data Flow

```
S3 Upload (FAST5) → Lambda Orchestrator → Step Functions → EC2 Phase Execution → Results to S3 → SNS Notifications
```

## EC2 Instance Types by Phase

| Phase | Instance Type | Purpose | Key Features |
|-------|--------------|---------|--------------|
| Phase 1 | g4dn.xlarge | GPU Basecalling | NVIDIA T4 GPU for Dorado |
| Phase 2 | t3.large | QC Analysis | CPU-based quality metrics |
| Phase 3 | r5.4xlarge | Host Removal | High memory for alignment |
| Phase 4 | Multiple t3.xlarge | Pathogen Detection | Parallel processing |
| Phase 5 | t3.large | Quantification | Calculation intensive |
| Phase 6 | t3.large | Report Generation | Document creation |

## Database Structure

### Reference Databases (EFS Mount)
- **Kraken2 DB**: `/mnt/efs/databases/kraken2/pmda_2024/`
- **BLAST DB**: `/mnt/efs/databases/blast/nt_2024/`
- **PERV DB**: `/mnt/efs/databases/perv/references/`
- **Custom PMDA DB**: `/mnt/efs/databases/pmda/2024.1/`

### RDS PostgreSQL Schema
- Pipeline runs metadata
- Sample information
- Pathogen detection results
- Quality metrics
- Report generation logs

## Security Architecture

1. **VPC Configuration**: Private subnets for EC2, public for NAT Gateway
2. **IAM Roles**: Least privilege per Lambda/EC2 function
3. **S3 Encryption**: AES-256 for all data at rest
4. **Secrets Manager**: Database credentials and API keys
5. **CloudTrail**: Audit logging for all API calls

## Cost Optimization Strategies

1. **Spot Instances**: 70% cost reduction for non-critical phases
2. **Auto-termination**: EC2 instances terminate after task completion
3. **EFS Intelligent Tiering**: Automatic cost optimization for reference DBs
4. **Lambda Concurrency Limits**: Prevent runaway costs
5. **S3 Lifecycle Policies**: Archive old run data to Glacier

## Monitoring and Alerting

### CloudWatch Metrics
- Pipeline execution time per phase
- EC2 instance utilization
- Failed run counts
- Cost per analysis

### SNS Alert Topics
- PERV detection (CRITICAL)
- Pipeline failures
- Cost threshold exceeded ($400)
- Quality metrics below threshold

## Scalability Considerations

- **Horizontal Scaling**: Phase 4 runs 4 parallel EC2 instances
- **Queue-based Processing**: SQS for batch sample processing
- **Database Scaling**: Aurora Serverless v2 auto-scales
- **EFS Performance Mode**: Max I/O for high throughput
- **Step Functions**: Handles complex workflow orchestration