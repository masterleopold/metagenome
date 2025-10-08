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

