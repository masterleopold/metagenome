# Quick Reference Guide

## Essential Commands

### Testing
```bash
pytest tests/                              # Run all tests
pytest tests/test_pmda_compliance.py -v   # Run specific test file
pytest tests/test_pmda_compliance.py::test_name -v  # Run single test
pytest --cov=scripts --cov=lambda tests/  # With coverage
```

### Code Quality
```bash
black scripts/ lambda/ tools/             # Format code
flake8 scripts/ lambda/                   # Lint check
mypy scripts/                             # Type checking
```

### Pipeline Operations
```bash
./tools/workflow_cli.py start --run-id RUN-2024-001 --bucket minion-data
./tools/workflow_cli.py status --run-id RUN-2024-001
./tools/workflow_cli.py logs --run-id RUN-2024-001 --phase 4
```

### AWS Operations
```bash
aws s3 ls s3://minion-data/runs/ --region ap-northeast-1
aws lambda invoke --function-name pipeline-orchestrator --payload '{}' response.json
aws ec2 describe-instances --filters "Name=tag:Pipeline,Values=MinION" --region ap-northeast-1
```

### Database Updates
```bash
# Update Kraken2 database
./scripts/database/update_kraken2.sh

# Update BLAST database
./scripts/database/update_blast_db.sh

# Update PERV reference database
./scripts/database/update_perv_db.sh
```

## Common Workflows

### 1. New Sample Processing
1. Upload FAST5 files to S3: `s3://minion-data/raw/{RUN_ID}/`
2. Start pipeline: `./tools/workflow_cli.py start --run-id {RUN_ID}`
3. Monitor progress: `./tools/workflow_cli.py status --run-id {RUN_ID}`
4. Review reports: `s3://minion-data/reports/{RUN_ID}/`

### 2. PERV Alert Response
1. Check alert details in SNS notification
2. Review typing results: `s3://minion-data/results/{RUN_ID}/phase4/perv_typing.json`
3. Generate detailed report: `./tools/perv_report.py --run-id {RUN_ID}`
4. Notify veterinary team immediately

### 3. Failed Run Recovery
1. Check logs: `./tools/workflow_cli.py logs --run-id {RUN_ID} --phase {PHASE}`
2. Identify failure point in Step Functions console
3. Fix issue (usually EC2 timeout or memory)
4. Restart from phase: `./tools/workflow_cli.py restart --run-id {RUN_ID} --from-phase {PHASE}`

### 4. Database Maintenance
1. Check current versions: `./tools/check_db_versions.sh`
2. Update if needed (monthly): `./scripts/database/update_all.sh`
3. Verify integrity: `./tools/verify_databases.py`
4. Update config: `templates/config/pmda_pathogens.json`

## Debug Commands

### Check EC2 Status
```bash
aws ec2 describe-instances \
  --filters "Name=instance-state-name,Values=running" \
  --query "Reservations[].Instances[].[InstanceId,Tags[?Key=='Name'].Value|[0],State.Name,LaunchTime]" \
  --output table
```

### View Lambda Logs
```bash
aws logs tail /aws/lambda/pipeline-orchestrator --follow
```

### S3 Data Verification
```bash
aws s3api head-object --bucket minion-data --key raw/{RUN_ID}/batch_001.fast5
```

### Step Functions Status
```bash
aws stepfunctions describe-execution \
  --execution-arn arn:aws:states:ap-northeast-1:ACCOUNT:execution:MinIONPipeline:{RUN_ID}
```

## Performance Optimization

### EC2 Instance Selection
- Phase 1 (Basecalling): `g4dn.xlarge` (GPU required)
- Phase 3 (Host Removal): `r5.4xlarge` (128GB RAM)
- Phase 4 (Pathogen): `c5.4xlarge` (4x parallel)
- Other phases: `t3.large` (sufficient)

### Cost Optimization
- Use Spot Instances: Add `--spot` flag to workflow_cli.py
- Schedule batch runs: Process multiple samples together
- Clean up S3: Remove intermediate files after 30 days

## Troubleshooting Quick Fixes

| Issue | Quick Fix |
|-------|-----------|
| EC2 timeout | Increase timeout in Lambda config |
| Out of memory | Use larger instance type |
| S3 access denied | Check IAM roles and bucket policies |
| Database not found | Verify EFS mount and paths |
| PERV false positive | Check contamination, rerun with clean sample |

## Important File Locations

### Configuration
- Pipeline config: `templates/config/pmda_pathogens.json`
- AWS config: `terraform/variables.tf`
- Database paths: `scripts/config/database_paths.json`

### Lambda Functions
- Orchestrator: `lambda/orchestration/pipeline_orchestrator.py`
- Phase handlers: `lambda/phase_*/handler.py`

### Analysis Scripts
- PERV typing: `scripts/phase4_pathogen/perv_typing.py`
- Quantification: `scripts/phase5_quantification/calculate_metrics.py`
- Report generation: `scripts/phase6_reports/generate_report.py`

### Databases (EFS Mount)
- Kraken2: `/mnt/efs/databases/kraken2/`
- BLAST: `/mnt/efs/databases/blast/`
- PERV: `/mnt/efs/databases/perv/`

## Emergency Contacts

- Pipeline Issues: DevOps team (via Slack #minion-pipeline)
- PERV Detection: Veterinary team (immediate notification required)
- AWS Issues: Cloud team (via PagerDuty)
- Database Updates: Bioinformatics team