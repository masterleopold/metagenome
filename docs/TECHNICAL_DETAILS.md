# Technical Details

## Pipeline Architecture (6 Phases)

**Implementation**: Lambda functions orchestrate EC2 instances (custom AMIs) for each phase. **No Docker containers used**.

### Phase Breakdown

1. **Basecalling** - FAST5→FASTQ (Dorado, GPU g4dn.xlarge)
2. **QC** - Quality assessment (NanoPlot/PycoQC, t3.large)
3. **Host Removal** - Sus scrofa depletion (Minimap2, r5.4xlarge)
4. **Pathogen Detection** - Multi-database screening (4 parallel EC2 instances)
5. **Quantification** - Abundance calculation (t3.large)
6. **Reporting** - PMDA-compliant reports (t3.large)

### Key Architecture Decisions

- Lambda functions trigger EC2 instances with UserData scripts
- Custom AMIs pre-installed with analysis tools
- EFS for shared reference databases (Kraken2, BLAST, PERV)
- EC2 instances auto-terminate after completion
- Spot Instances for 70% cost savings

## Implementation Phases

### Current: Phase 2 (Years 1-3) - Analysis In-house Implementation
- **In-house**: Bioinformatics analysis
- **Outsourced**: Sequencing only
- **Target**: 12-24 analyses/year

### Phase 1 (Year 0) - Initial Setup
- Infrastructure deployment
- Reference database creation
- Pipeline validation

### Phase 3 (Years 3-5) - Scale-up
- Increase to 50+ analyses/year
- Full automation
- Cost optimization

## Database Architecture

### PMDA Pathogen Database
- Location: `/mnt/efs/databases/pmda/2024.1/`
- Coverage: 91 PMDA-designated pathogens
- Update cycle: Quarterly

### PERV Detection Database
- Custom curated PERV-A/B/C sequences
- Recombinant detection algorithms
- Real-time SNS alerting

### Host Reference
- Sus scrofa genome (Sscrofa11.1)
- CpG methylation patterns for depletion

## Performance Requirements

- **Basecalling**: Q30+ accuracy (duplex mode)
- **Host Removal**: >90% efficiency
- **Detection**: PPA >95%, NPA >98%, R² >0.90
- **Turnaround**: 48-72 hours total
- **Cost Cap**: $400/sample alert threshold

## Scalability Considerations

- Auto-scaling EC2 fleet (1-100 concurrent analyses)
- S3 lifecycle policies (5-year retention)
- CloudWatch monitoring and alerting
- Step Functions for workflow orchestration