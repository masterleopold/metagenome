# Commands Reference

Complete command reference for the MinION Pathogen Screening Pipeline.

## Environment Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Configure AWS
aws configure
# Set AWS_REGION=ap-northeast-1 (Tokyo)
export AWS_REGION=ap-northeast-1
export AWS_DEFAULT_REGION=ap-northeast-1
```

## Testing Commands

```bash
# Run all tests
python -m pytest tests/

# Run unit tests only (exclude integration tests)
python -m pytest tests/ -k "not integration"

# Run specific test file
python -m pytest tests/test_pmda_compliance.py

# Run specific test function with verbose output
python -m pytest tests/test_pmda_compliance.py::test_function_name -v

# Run tests with coverage
python -m pytest tests/ --cov=scripts --cov=lambda

# Run PMDA compliance tests specifically
python -m pytest tests/test_pmda_compliance.py -v

# Run pathogen detection tests
python -m pytest tests/test_pathogen_detection.py -v
```

## Pipeline Operations

```bash
# Start a new pipeline run
./tools/workflow_cli.py start \
  --run-id RUN-2024-001 \
  --bucket minion-data \
  --input-prefix runs/RUN-2024-001/fast5/

# Check pipeline status
./tools/workflow_cli.py status --run-id RUN-2024-001

# Watch status in real-time
./tools/workflow_cli.py status --run-id RUN-2024-001 --watch

# Get pipeline metrics
./tools/workflow_cli.py metrics --run-id RUN-2024-001

# List all workflows
./tools/workflow_cli.py list

# Cancel a running pipeline
./tools/workflow_cli.py cancel --run-id RUN-2024-001
```

## Code Quality

```bash
# Format code with Black
black scripts/ lambda/ tools/ tests/

# Check specific file
black scripts/phase4_pathogen/perv_typing.py --check

# Lint with flake8
flake8 scripts/ lambda/ tools/ tests/

# Lint specific directory
flake8 scripts/phase4_pathogen/

# Type checking with mypy
mypy scripts/
mypy lambda/
mypy tools/

# Run all quality checks
black scripts/ lambda/ tools/ tests/ && \
flake8 scripts/ lambda/ tools/ tests/ && \
mypy scripts/
```

## Documentation Portal

```bash
# Start development server
cd docs-portal
npm install  # First time only
npm run dev  # Runs on http://localhost:3003

# Build for production
cd docs-portal
npm run build

# Run tests
cd docs-portal
npm test

# Check for TypeScript errors
cd docs-portal
npm run type-check
```

## Database Management

```bash
# Setup all reference databases
./tools/database_setup.sh --all

# Update PMDA pathogen database
./tools/database_setup.sh --pmda

# Update Kraken2 database
./tools/database_setup.sh --kraken2

# Update BLAST database
./tools/database_setup.sh --blast

# Check database integrity
./tools/database_check.py --verify
```

## AWS Operations

```bash
# Deploy infrastructure with Terraform
cd infrastructure/terraform
terraform init
terraform plan
terraform apply

# Update Lambda functions
./tools/deployment_script.sh deploy-lambda

# Update Step Functions
./tools/deployment_script.sh update-stepfunctions

# Check EC2 instances
aws ec2 describe-instances --filters "Name=tag:Project,Values=minion-pipeline"

# View CloudWatch logs
aws logs tail /aws/lambda/minion-pipeline-production --follow

# Check S3 buckets
aws s3 ls s3://minion-data/
aws s3 ls s3://minion-results/
```

## Monitoring

```bash
# Launch monitoring dashboard
streamlit run tools/monitoring_dashboard.py

# Get cost report
./tools/cost_analyzer.py --month 2024-11

# Check alerts
./tools/check_alerts.py --last 24h

# View metrics
./tools/metrics_viewer.py --phase all --run-id RUN-2024-001
```

## Development Tools

```bash
# Create new phase script
./tools/create_phase_script.py --phase 7 --name "validation"

# Generate test data
./tools/generate_test_data.py --type fastq --reads 10000

# Validate configuration
./tools/config_validator.py templates/config/pmda_pathogens.json

# Run local simulation
./tools/local_simulator.py --phases 1,2,3 --test-data sample.fastq
```

## Troubleshooting

```bash
# Check pipeline logs
./tools/get_logs.py --run-id RUN-2024-001 --phase all

# Debug failed phase
./tools/debug_phase.py --run-id RUN-2024-001 --phase 4

# Retry failed phase
./tools/retry_phase.py --run-id RUN-2024-001 --phase 4

# Export debug info
./tools/export_debug.py --run-id RUN-2024-001 --output debug.tar.gz
```

## Quick Scripts

```bash
# Run PERV detection on specific file
python scripts/phase4_pathogen/perv_typing.py \
  --bam aligned.bam \
  --output perv_results.json

# Generate PMDA report
python scripts/phase6_reports/generate_pmda_report.py \
  --results results.json \
  --output report.pdf

# Check QC metrics
python scripts/phase2_qc/check_qc_metrics.py \
  --fastq reads.fastq \
  --min-q30 0.85

# Calculate host depletion rate
python scripts/phase3_host_removal/calculate_depletion_rate.py \
  --before raw.fastq \
  --after filtered.fastq
```