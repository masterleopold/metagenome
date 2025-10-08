## Notes for Future Development

When working with this codebase:

1. **Maintain consistency** with the phased implementation approach (Phases 1-4 over 5+ years)
2. **Ensure PMDA compliance** - all technical specifications must align with PMDA regulatory requirements
3. **Update cost projections** to reflect current market conditions (especially for reagents and cloud services)
4. **Consider scalability** when proposing technical solutions (current focus: 24 samples/year for Phase I clinical trial)
5. **Document validation strategies** for any new methods - must include LOD, reproducibility, specificity, PPA/NPA metrics
6. **Maintain bilingual capability** (Japanese primary, English technical terms)
7. **Cross-reference documents** - when updating protocols, check if strategic documents need corresponding updates
8. **Version control** - major protocol changes should be tracked with version numbers (e.g., v2.0 for strategy shift from PCR-hybrid to NGS-only)
9. **Reference the validation framework** - all protocol modifications should consider impact on the MinION vs MiSeq validation study

### Coding Standards

**Python Code Requirements:**
- Use Google-style docstrings for all public functions and classes
- Type hints required for function signatures
- Follow PEP 8 style guide (enforced by black + flake8)
- Maximum line length: 88 characters (black default)
- Use `logging` module, never `print()` statements in production code

**Example Function:**
```python
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

    Raises:
        FileNotFoundError: If FASTQ file does not exist
        ValueError: If FASTQ file is malformed

    Example:
        >>> detector = PathogenDetector(config)
        >>> results = detector.detect_pathogens("sample.fastq")
        >>> print(results["PERV-A"])
        0.98
    """
    logger.info(f"Starting pathogen detection on {fastq_path}")
    # Implementation
```

**Shell Scripts:**
- Use `#!/usr/bin/env bash` shebang
- Set strict mode: `set -euo pipefail`
- Include header with description, author, version
- Implement error handling with trap
- Use readonly for constants
- Quote all variable expansions

**Terraform:**
- Use descriptive resource names: `${var.project_name}-{resource}-${var.environment}`
- Include comprehensive tags (Name, Environment, Purpose, ManagedBy, CostCenter, Owner)
- Enable monitoring and IMDSv2 for EC2 instances
- Document dependencies in comments
- Use `lifecycle` blocks to prevent unwanted resource replacement

**Git Commit Messages:**
Follow conventional commits specification:
```bash
feat(pathogen): add PERV-C detection algorithm
fix(basecalling): correct duplex mode parameters
docs(api): update endpoint documentation
perf(analysis): optimize Kraken2 memory usage
test(compliance): add PMDA validation tests
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`, `revert`

**Pull Request Requirements:**
- Minimum 2 approvals required
- PMDA-related changes need compliance team review
- All CI/CD checks must pass
- Code coverage must not decrease
- No new security vulnerabilities
- Signed commits (GPG) required

For complete contribution guidelines, see `CONTRIBUTING.md`
