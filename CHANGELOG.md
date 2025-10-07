# Changelog

All notable changes to the MinION Pathogen Screening Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-01-15

### Added
- Initial release of MinION Pathogen Screening Pipeline
- PMDA-compliant screening for 91 designated pathogens
- PERV-specific detection module with critical alerts
- AWS serverless architecture (Lambda + EC2 on-demand)
- Terraform infrastructure as code for AWS deployment
- Step Functions workflow orchestration
- Real-time basecalling with Dorado duplex mode (Q30 accuracy)
- Multi-database pathogen detection (Kraken2, RVDB, BLAST)
- Spike-in normalization for absolute quantification
- Comprehensive reporting in PDF, JSON, and HTML formats
- CloudWatch monitoring and alerting
- Cost optimization with spot instances and auto-scaling
- S3 event-driven Lambda triggers
- RDS Aurora Serverless for metadata storage
- EFS for reference database storage
- API Gateway for RESTful API access
- Complete test suite with PMDA compliance validation

### Infrastructure Components
- **VPC Configuration**: Multi-AZ deployment with public/private subnets
- **S3 Buckets**: Data storage with lifecycle policies
- **Lambda Functions**: 16 serverless functions for orchestration
- **EC2 Instances**: GPU-enabled for basecalling, high-memory for analysis
- **Step Functions**: 6-phase analysis workflow
- **RDS Aurora**: Serverless PostgreSQL for metadata
- **EFS**: Shared storage for reference databases
- **CloudWatch**: Comprehensive monitoring and logging
- **SNS**: Alert notifications for critical findings
- **API Gateway**: RESTful API with authentication

### Analysis Pipeline Phases
1. **Basecalling**: FAST5/POD5 to FASTQ conversion with Q30 accuracy
2. **Quality Control**: Read quality assessment and filtering
3. **Host Removal**: Sus scrofa genome alignment and removal
4. **Pathogen Detection**: Multi-database screening
5. **Quantification**: Absolute copy number calculation
6. **Report Generation**: PMDA-compliant reporting

### Database Components
- Kraken2 Standard Database
- RVDB v30.0 (viral sequences)
- PMDA Custom Database (91 pathogens)
- Sus scrofa reference genome (Sscrofa11.1)
- PERV-specific sequence database

### Documentation
- Comprehensive deployment guide
- API documentation with examples
- PMDA compliance testing framework
- Development plan with 68 files across 7 phases

## [0.9.0] - 2024-01-01 (Pre-release)

### Added
- Beta version for internal testing
- Core basecalling functionality
- Basic pathogen detection
- Initial AWS infrastructure

### Changed
- Migrated from Docker to serverless architecture
- Optimized for Lambda + EC2 on-demand approach

### Fixed
- GPU instance availability issues
- Memory limitations in Lambda functions

## [0.8.0] - 2023-12-15 (Alpha)

### Added
- Proof of concept implementation
- Manual workflow execution
- Basic Kraken2 integration

### Known Issues
- No PERV-specific detection
- Limited to single database search
- Manual result integration required

## Future Releases

### [1.1.0] - Planned Q2 2024
- [ ] Machine learning model for pathogen prediction
- [ ] Real-time streaming analysis
- [ ] Mobile app for monitoring
- [ ] Integration with LIMS systems
- [ ] Support for GridION/PromethION

### [1.2.0] - Planned Q3 2024
- [ ] AI-powered unknown pathogen detection
- [ ] Automated reanalysis with updated databases
- [ ] Multi-language report generation (Japanese/English)
- [ ] Integration with national pathogen surveillance systems

### [2.0.0] - Planned 2025
- [ ] Complete redesign for clinical trial phase
- [ ] GLP compliance features
- [ ] 21 CFR Part 11 compliance
- [ ] Integration with electronic health records
- [ ] Support for human clinical samples

## Versioning Strategy

- **Major versions (X.0.0)**: Breaking changes, new regulatory compliance
- **Minor versions (0.X.0)**: New features, backwards compatible
- **Patch versions (0.0.X)**: Bug fixes, security updates

## Support Policy

- **Current version (1.0.0)**: Full support
- **Previous minor version**: Security updates only
- **Older versions**: No support

## Migration Notes

### From 0.9.0 to 1.0.0
1. Update Terraform to v1.0+
2. Migrate RDS database schema
3. Update Lambda function permissions
4. Rebuild AMIs with latest Dorado
5. Update API endpoints

### From Docker to Serverless
1. Export data from Docker volumes
2. Deploy new infrastructure with Terraform
3. Migrate data to S3 and RDS
4. Update workflow configurations
5. Validate with test dataset

## Contributors

- Architecture Team: Infrastructure design and implementation
- Bioinformatics Team: Pipeline development and validation
- DevOps Team: Deployment and monitoring
- QA Team: Testing and compliance validation

## License

Proprietary - All rights reserved

## Contact

- Issues: GitHub Issues
- Support: support@your-org.com
- Security: security@your-org.com