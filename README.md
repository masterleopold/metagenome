# MinION Pathogen Screening Pipeline

PMDA-compliant metagenomic analysis pipeline for xenotransplantation donor pig screening using Oxford Nanopore MinION Mk1D.

## Overview

This pipeline provides comprehensive pathogen detection for 91 PMDA-designated pathogens in donor pigs intended for xenotransplantation. The system uses Oxford Nanopore long-read sequencing technology combined with AWS cloud infrastructure for scalable, cost-effective analysis.

## Key Features

- **PMDA Compliance**: Full coverage of 91 designated pathogens
- **PERV Detection**: Critical detection of Porcine Endogenous Retroviruses (PERV-A, B, C)
- **Real-time Analysis**: Streaming analysis capability with MinION
- **Cloud-Native**: Serverless architecture on AWS (Lambda + EC2 on-demand)
- **Automated Workflow**: End-to-end automation from basecalling to reporting
- **Quality Assured**: Q30+ accuracy with duplex basecalling
- **Cost Optimized**: Spot instances and on-demand scaling

## System Architecture

```
MinION Sequencer → S3 Upload → Lambda Orchestration → EC2 Processing → Reports
                                        ↓
                              Step Functions Workflow
                                        ↓
                    [Basecalling → QC → Host Removal → Pathogen Detection →
                     Quantification → Report Generation]
```

## Requirements

### Hardware
- Oxford Nanopore MinION Mk1D
- GPU-enabled EC2 instances (g4dn.xlarge for basecalling)
- High-memory EC2 instances (r5.4xlarge for pathogen detection)

### Software
- Python 3.9+
- Terraform 1.0+
- AWS CLI configured
- Docker (optional)

### AWS Services
- S3 for data storage
- Lambda for orchestration
- EC2 for compute
- Step Functions for workflow
- RDS Aurora Serverless for metadata
- EFS for reference databases
- SNS for notifications

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/your-org/minion-pipeline.git
cd minion-pipeline
```

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

### 3. Configure AWS

```bash
aws configure
export AWS_REGION=ap-northeast-1
export ENVIRONMENT=production
```

### 4. Deploy Infrastructure

```bash
cd infrastructure/terraform
terraform init
terraform plan
terraform apply
```

### 5. Setup Databases

```bash
./tools/database_setup.sh --all
```

### 6. Deploy Pipeline

```bash
./tools/deployment_script.sh deploy
```

## Usage

### Starting a Workflow

```bash
# Using CLI
./tools/workflow_cli.py start \
  --run-id RUN-2024-001 \
  --bucket minion-data \
  --input-prefix runs/RUN-2024-001/fast5/

# Using API
curl -X POST https://api.your-domain.com/workflows \
  -H "x-api-key: YOUR_API_KEY" \
  -d '{
    "run_id": "RUN-2024-001",
    "bucket": "minion-data",
    "input_prefix": "runs/RUN-2024-001/fast5/"
  }'
```

### Monitoring Progress

```bash
# Check status
./tools/workflow_cli.py status --run-id RUN-2024-001 --watch

# View metrics
./tools/workflow_cli.py metrics --run-id RUN-2024-001

# Launch dashboard
streamlit run tools/monitoring_dashboard.py
```

### Accessing Reports

Reports are automatically generated in multiple formats:
- **PDF**: Comprehensive report with visualizations
- **JSON**: Machine-readable PMDA checklist
- **HTML**: Interactive web report

## Pipeline Phases

### 1. Basecalling
- Converts raw signal (FAST5/POD5) to sequences (FASTQ)
- Uses Dorado with duplex mode for Q30 accuracy
- GPU-accelerated on g4dn instances

### 2. Quality Control
- Assesses read quality metrics
- Filters low-quality reads
- Generates QC reports with NanoPlot

### 3. Host Genome Removal
- Aligns reads to Sus scrofa reference genome
- Removes host DNA contamination
- Retains unmapped (potentially pathogenic) reads

### 4. Pathogen Detection
- **Kraken2**: Rapid taxonomic classification
- **RVDB**: Viral database search
- **BLAST**: PMDA custom database alignment
- **PERV-specific**: Targeted PERV detection and typing

### 5. Quantification
- Spike-in normalization (PhiX174)
- Absolute copy number calculation (copies/mL)
- Confidence interval estimation

### 6. Report Generation
- PMDA compliance checklist (91 pathogens)
- Detection summary with risk assessment
- Quality metrics and validation

## Configuration

### Default Configuration

See `templates/config/default_pipeline.yaml` for full configuration options.

```yaml
phases:
  basecalling:
    skip_duplex: false  # Use duplex for Q30
    min_quality: 9

  pathogen_detection:
    databases:
      - kraken2
      - rvdb
      - pmda
    confidence_threshold: 0.1

  reporting:
    formats:
      - pdf
      - json
    pmda_checklist: true
```

### Custom Configuration

```bash
# Create custom config
./tools/workflow_cli.py config --create-default

# Edit config
vi /etc/minion-pipeline/custom.yaml

# Validate
./tools/workflow_cli.py config --validate custom.yaml
```

## PMDA Compliance

### 91 Pathogen Coverage

The pipeline screens for all 91 PMDA-designated pathogens:
- **Viruses**: 50 pathogens including PERV variants
- **Bacteria**: 35 pathogens
- **Parasites**: 5 pathogens
- **Fungi**: 5 pathogens
- **Prions**: 1 pathogen

### Critical Pathogens

Immediate alerts are triggered for:
- PERV-A, PERV-B, PERV-C (xenotransplantation critical)
- ASFV, CSFV, FMDV (disease outbreak risks)
- Prion detection

### Reporting Requirements

All reports include:
- Complete 91 pathogen checklist
- PERV-specific analysis section
- Quantitative results (copies/mL)
- Detection confidence levels
- Quality control metrics

## Testing

```bash
# Run unit tests
python -m pytest tests/

# Run integration tests
python tests/test_pipeline_integration.py

# Run PMDA compliance tests
python tests/test_pmda_compliance.py
```

## Monitoring

### CloudWatch Dashboard

Access the dashboard at:
```
https://console.aws.amazon.com/cloudwatch/home?region=ap-northeast-1#dashboards:name=minion-pipeline-production
```

### Metrics

Key metrics tracked:
- Workflow execution time
- Phase completion rates
- Pathogen detection counts
- Resource utilization
- Cost per analysis

### Alerts

Configured alerts:
- PERV detection (CRITICAL)
- High error rate
- Long-running workflows
- Cost threshold exceeded

## Cost Optimization

### Spot Instances
- Basecalling uses GPU spot instances (70% cost reduction)
- Automatic fallback to on-demand if spot unavailable

### Auto-scaling
- EC2 instances auto-terminate after phase completion
- Lambda functions scale automatically

### Data Lifecycle
- Raw data: 30-day retention
- Processed data: 90-day retention
- Reports: 365-day retention

## Troubleshooting

### Common Issues

1. **Basecalling Timeout**
   - Check GPU availability
   - Verify CUDA drivers
   - Increase instance size

2. **High Host Contamination**
   - Review extraction protocol
   - Check depletion efficiency
   - Adjust alignment parameters

3. **Low Pathogen Detection**
   - Verify database integrity
   - Check confidence thresholds
   - Review sample quality

### Logs

Access logs via:
```bash
# CloudWatch logs
aws logs tail /aws/lambda/minion-pipeline-production --follow

# EC2 instance logs
aws ssm start-session --target INSTANCE_ID
```

## Support

- **Documentation**: See `/docs` directory
- **Issues**: GitHub Issues
- **Contact**: support@your-org.com

## License

Proprietary - All rights reserved

## Citations

If using this pipeline, please cite:
- Oxford Nanopore Technologies
- Kraken2 (Wood et al., 2019)
- RVDB (Goodacre et al., 2018)
- PMDA Xenotransplantation Guidelines (2024)

## Version

Current Version: 1.0.0

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.