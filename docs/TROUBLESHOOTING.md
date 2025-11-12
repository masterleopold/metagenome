# Troubleshooting Guide

## Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| **Basecalling timeout** | Check GPU availability, increase `max_runtime_hours` in config, consider using g4dn.2xlarge for larger datasets |
| **High host contamination** | Review DNA extraction protocol, check Minimap2 parameters, verify CpG methylation-based depletion is working |
| **No pathogen detection** | Verify reference databases are mounted correctly, check confidence thresholds (default 0.1), ensure sufficient sequencing depth |
| **EC2 launch failure** | Check spot instance availability in region, enable fallback to on-demand instances, verify AMI permissions |
| **Memory errors in QC** | Reduce batch size in NanoPlot, use streaming mode for large files, consider r5.xlarge instance |
| **PERV alert not triggering** | Verify SNS topic configuration, check IAM permissions, test with known PERV-positive control |
| **Slow pathogen detection** | Ensure 4 parallel EC2 instances are launching, check EFS throughput limits, optimize database indexing |
| **Report generation failure** | Check template files exist, verify all required fields in metadata, ensure S3 write permissions |
| **Cost exceeding threshold** | Review spot instance usage, optimize EC2 instance types, check for stuck instances not auto-terminating |
| **Database update failures** | Verify EFS mount points, check available storage space, ensure update scripts have correct permissions |

## Advanced Troubleshooting

### Debugging Lambda Functions
```bash
aws logs tail /aws/lambda/minion-orchestrator --follow
aws lambda invoke --function-name minion-orchestrator response.json --log-type Tail
```

### Monitoring EC2 Instances
```bash
aws ec2 describe-instances --filters "Name=tag:Pipeline,Values=MinION" --query 'Reservations[*].Instances[*].[InstanceId,State.Name,LaunchTime]'
```

### Checking Step Functions Execution
```bash
aws stepfunctions describe-execution --execution-arn <execution-arn>
```

## Contact Support

For issues not covered here, see:
- [MinION Protocol Troubleshooting](../md/MinION_Protocol_00_目次とマスタードキュメント.md) - Appendix C
- [GitHub Issues](https://github.com/masterleopold/metagenome/issues)
- AWS Support Console for infrastructure issues