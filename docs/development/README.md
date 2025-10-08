# Development Documentation

This directory contains detailed technical documentation for developers working on the MinION pathogen screening pipeline.

## Available Guides

### [WORKFLOW_GUIDE.md](./WORKFLOW_GUIDE.md)
Complete development workflow including:
- Setup and installation procedures
- Running the pipeline
- Testing strategies
- Configuration management
- Development tools and commands
- Code architecture overview
- Database management
- API integration
- Monitoring and debugging
- Working with PERV detection

### [CODING_STANDARDS.md](./CODING_STANDARDS.md)
Project coding standards and conventions:
- Python code requirements (Google-style docstrings, PEP 8)
- Shell script standards
- Terraform best practices
- Git commit message format (conventional commits)
- Pull request requirements
- Code review guidelines

## Quick Links

- [Main CLAUDE.md](../../CLAUDE.md) - Quick reference guide
- [API Documentation](../API_DOCUMENTATION.md) - REST API endpoints
- [Deployment Guide](../DEPLOYMENT_GUIDE.md) - Infrastructure setup
- [Session History](../claude-sessions/README.md) - Development session logs

## Project Structure

```
metagenome/
├── scripts/              # Pipeline phase scripts
│   ├── phase1_basecalling/
│   ├── phase2_qc/
│   ├── phase3_host_removal/
│   ├── phase4_pathogen/  # CRITICAL: PERV detection scripts here
│   ├── phase5_quantification/
│   └── phase6_reports/
├── lambda/              # AWS Lambda functions
├── lib/                 # Shared Python libraries
├── infrastructure/      # Terraform IaC
├── tests/              # Test suite
└── docs-portal/        # Next.js documentation site
```

## Key Development Principles

1. **PMDA Compliance First** - All changes must maintain regulatory compliance
2. **PERV Detection Priority** - Any PERV-related code requires extra review
3. **Test Coverage** - Maintain >80% code coverage
4. **Documentation** - Update docs with any API or workflow changes
5. **Security** - No patient data in repositories, use encrypted storage

## Getting Help

- Review [Session History](../claude-sessions/README.md) for past problem solutions
- Check [WORKFLOW_GUIDE.md](./WORKFLOW_GUIDE.md) for common operations
- Submit issues to [GitHub](https://github.com/masterleopold/metagenome/issues)