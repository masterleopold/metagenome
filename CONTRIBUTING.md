# Contributing to MinION Pathogen Screening Pipeline

Thank you for your interest in contributing to the MinION Pathogen Screening Pipeline project. This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Contribution Process](#contribution-process)
- [Coding Standards](#coding-standards)
- [Testing Requirements](#testing-requirements)
- [Documentation Standards](#documentation-standards)
- [Security Guidelines](#security-guidelines)
- [PMDA Compliance](#pmda-compliance)

## Code of Conduct

### Our Standards

- Use welcoming and inclusive language
- Be respectful of differing viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the project and compliance requirements
- Show empathy towards other contributors

### Unacceptable Behavior

- Harassment, discrimination, or offensive comments
- Publishing private information without consent
- Unprofessional conduct that could affect regulatory compliance
- Any behavior that violates company policies

## Getting Started

### Prerequisites

Before contributing, ensure you have:

1. **Technical Requirements**
   - Python 3.9+
   - AWS CLI configured
   - Terraform 1.0+
   - Git with GPG signing configured
   - Access to AWS development environment

2. **Knowledge Requirements**
   - Understanding of NGS data analysis
   - Familiarity with AWS services
   - Knowledge of PMDA xenotransplantation guidelines
   - Experience with bioinformatics pipelines

3. **Access Requirements**
   - GitHub repository access
   - AWS IAM credentials
   - Access to reference databases
   - Slack/Teams for communication

## Development Setup

### 1. Clone the Repository

```bash
git clone https://github.com/your-org/minion-pipeline.git
cd minion-pipeline
git checkout -b feature/your-feature-name
```

### 2. Set Up Virtual Environment

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

### 3. Configure Pre-commit Hooks

```bash
pre-commit install
pre-commit run --all-files  # Test hooks
```

### 4. Set Up Local Testing Environment

```bash
# Download test data
./scripts/download_test_data.sh

# Set up local databases
./tools/database_setup.sh --local

# Configure environment variables
cp .env.example .env
# Edit .env with your settings
```

## Contribution Process

### 1. Issue Creation

Before starting work:

```markdown
## Issue Template

**Type**: Bug/Feature/Enhancement/Documentation

**Priority**: Critical/High/Medium/Low

**PMDA Impact**: Yes/No

**Description**:
[Clear description of the issue or feature]

**Acceptance Criteria**:
- [ ] Criterion 1
- [ ] Criterion 2

**Testing Requirements**:
- [ ] Unit tests
- [ ] Integration tests
- [ ] PMDA compliance validation
```

### 2. Branch Strategy

```bash
# Feature branches
git checkout -b feature/pmda-pathogen-update

# Bugfix branches
git checkout -b bugfix/perv-detection-threshold

# Hotfix branches (critical production issues)
git checkout -b hotfix/critical-alert-failure
```

### 3. Commit Guidelines

Follow conventional commits specification:

```bash
# Format: <type>(<scope>): <subject>

# Examples:
git commit -m "feat(pathogen): add PERV-C detection algorithm"
git commit -m "fix(basecalling): correct duplex mode parameters"
git commit -m "docs(api): update endpoint documentation"
git commit -m "perf(analysis): optimize Kraken2 memory usage"
git commit -m "test(compliance): add PMDA validation tests"
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `test`: Test additions or fixes
- `chore`: Maintenance tasks
- `revert`: Reverting changes

### 4. Pull Request Process

#### PR Template

```markdown
## Description
[Brief description of changes]

## Related Issue
Fixes #123

## Type of Change
- [ ] Bug fix (non-breaking change)
- [ ] New feature (non-breaking change)
- [ ] Breaking change
- [ ] Documentation update
- [ ] Performance improvement

## PMDA Compliance Impact
- [ ] No impact on PMDA compliance
- [ ] Updates PMDA pathogen list
- [ ] Modifies detection algorithms
- [ ] Changes reporting format

## Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] PMDA compliance tests pass
- [ ] Manual testing completed

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Comments added for complex logic
- [ ] Documentation updated
- [ ] No new warnings generated
- [ ] Dependent changes merged
```

#### Review Process

1. **Automated Checks**
   - CI/CD pipeline must pass
   - Code coverage must not decrease
   - No security vulnerabilities

2. **Peer Review**
   - Minimum 2 approvals required
   - PMDA-related changes need compliance team review
   - Architecture changes need tech lead approval

3. **Merge Requirements**
   - All conversations resolved
   - Branch up to date with main
   - Signed commits (GPG)

## Coding Standards

### Python Style Guide

```python
"""
Module docstring following Google style.

This module handles pathogen detection using multiple databases.
"""

from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


class PathogenDetector:
    """
    Detects pathogens in metagenomic data.

    Attributes:
        config: Configuration dictionary
        databases: List of database paths
    """

    def __init__(self, config: Dict[str, any]) -> None:
        """
        Initialize the PathogenDetector.

        Args:
            config: Configuration dictionary containing:
                - database_path: Path to pathogen database
                - confidence_threshold: Minimum confidence score

        Raises:
            ValueError: If config is missing required fields
        """
        self.config = self._validate_config(config)
        self.databases = self._load_databases()

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

        Example:
            >>> detector = PathogenDetector(config)
            >>> results = detector.detect_pathogens("sample.fastq")
            >>> print(results["PERV-A"])
            0.98
        """
        # Implementation
        pass
```

### Shell Script Standards

```bash
#!/usr/bin/env bash
#
# Script: basecall_duplex.sh
# Description: Run Dorado duplex basecalling for Q30 accuracy
# Author: Your Name
# Date: 2024-01-15
# Version: 1.0.0

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# Configuration
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly LOG_FILE="/var/log/minion/basecalling.log"

# Function: Log message with timestamp
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

# Function: Handle errors
error_handler() {
    local line_num=$1
    log "ERROR: Script failed at line ${line_num}"
    exit 1
}

trap 'error_handler ${LINENO}' ERR

# Main execution
main() {
    log "Starting basecalling process"

    # Validate inputs
    if [[ $# -lt 2 ]]; then
        log "ERROR: Usage: $0 <input_dir> <output_dir>"
        exit 1
    fi

    local input_dir="$1"
    local output_dir="$2"

    # Process files
    # ...
}

# Run main function
main "$@"
```

### Terraform Standards

```hcl
# Resource: EC2 instance for basecalling
# Purpose: GPU-enabled instance for Dorado basecalling
# Dependencies: VPC, Security Group, IAM Role

resource "aws_instance" "basecalling" {
  # Use consistent naming convention
  tags = {
    Name        = "${var.project_name}-basecalling-${var.environment}"
    Environment = var.environment
    Purpose     = "Basecalling"
    ManagedBy   = "Terraform"
    CostCenter  = var.cost_center
    Owner       = var.owner_email
  }

  # Instance configuration
  ami           = data.aws_ami.gpu_optimized.id
  instance_type = var.basecalling_instance_type  # g4dn.xlarge

  # Best practices
  monitoring                  = true
  ebs_optimized              = true
  associate_public_ip_address = false

  metadata_options {
    http_endpoint = "enabled"
    http_tokens   = "required"  # IMDSv2
  }

  # User data for initialization
  user_data = templatefile("${path.module}/templates/basecalling_init.sh", {
    region     = var.aws_region
    bucket     = aws_s3_bucket.data.id
    efs_dns    = aws_efs_file_system.databases.dns_name
  })

  lifecycle {
    ignore_changes = [ami]  # Prevent replacement on AMI updates
  }
}
```

## Testing Requirements

### Test Categories

1. **Unit Tests** (Required for all functions)
   ```python
   python -m pytest tests/unit/ --cov=src --cov-report=term-missing
   ```

2. **Integration Tests** (Required for pipelines)
   ```python
   python -m pytest tests/integration/ -v
   ```

3. **PMDA Compliance Tests** (Required for pathogen detection)
   ```python
   python tests/test_pmda_compliance.py
   ```

4. **Performance Tests** (For optimization PRs)
   ```python
   python tests/performance/test_throughput.py
   ```

### Test Coverage Requirements

- Minimum overall coverage: 80%
- Critical paths (PERV detection): 95%
- New features: 90%
- Bug fixes: Include regression tests

### Test Data

- Use provided test datasets in `tests/data/`
- Never use real patient/animal data in tests
- Mock external API calls
- Use fixtures for database connections

## Documentation Standards

### Code Documentation

- All public functions must have docstrings
- Complex algorithms need inline comments
- Configuration options must be documented
- API endpoints need OpenAPI specification

### User Documentation

- Update README.md for user-facing changes
- Update API_DOCUMENTATION.md for API changes
- Update DEPLOYMENT_GUIDE.md for infrastructure changes
- Add migration guides for breaking changes

### Technical Documentation

Location: `docs/technical/`

- Architecture decisions
- Database schemas
- Algorithm explanations
- Performance benchmarks

## Security Guidelines

### Sensitive Data

- Never commit credentials, keys, or tokens
- Use AWS Secrets Manager for secrets
- Encrypt sensitive configuration
- Follow principle of least privilege

### Code Security

```python
# DO NOT hardcode credentials
# BAD:
api_key = "sk-abc123def456"  # Never do this

# GOOD:
import boto3
secrets = boto3.client('secretsmanager')
api_key = secrets.get_secret_value(SecretId='api-key')['SecretString']
```

### Dependency Management

- Review dependencies for vulnerabilities
- Keep dependencies up to date
- Use specific version pins
- Document security advisories

## PMDA Compliance

### Critical Requirements

1. **Pathogen Coverage**
   - All 91 PMDA pathogens must be detectable
   - PERV detection is highest priority
   - Cannot remove pathogens without approval

2. **Detection Accuracy**
   - Maintain >95% sensitivity for critical pathogens
   - Document any changes to detection algorithms
   - Validate changes with reference samples

3. **Reporting Compliance**
   - Reports must include all required fields
   - PMDA checklist format cannot change
   - Maintain audit trail for all results

### Validation Process

Before merging PMDA-related changes:

1. Run full compliance test suite
2. Generate test reports
3. Review with compliance team
4. Document validation results
5. Update validation records

## Review Checklist

### For Contributors

- [ ] Code follows style guidelines
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] No hardcoded values
- [ ] Error handling implemented
- [ ] Logging added for debugging
- [ ] Performance impact considered
- [ ] Security review completed
- [ ] PMDA compliance verified

### For Reviewers

- [ ] Code quality acceptable
- [ ] Tests comprehensive
- [ ] Documentation clear
- [ ] No security issues
- [ ] Performance acceptable
- [ ] Error handling appropriate
- [ ] PMDA requirements met
- [ ] Breaking changes documented

## Getting Help

### Resources

- Project Wiki: Internal documentation
- Slack Channel: #minion-pipeline
- Office Hours: Thursdays 2-3 PM JST
- Email: minion-dev@your-org.com

### Key Contacts

- Technical Lead: tech-lead@your-org.com
- PMDA Compliance: compliance@your-org.com
- Security Team: security@your-org.com
- DevOps Team: devops@your-org.com

## Recognition

Contributors will be recognized in:
- CHANGELOG.md for significant contributions
- Annual contributor report
- Team presentations
- Internal awards program

Thank you for contributing to the MinION Pathogen Screening Pipeline!