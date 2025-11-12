# Development Guide

This guide contains detailed development commands, code conventions, and patterns for the MinION Pathogen Screening Pipeline.

## Development Commands

### Environment Setup
```bash
pip install -r requirements.txt
aws configure  # Set AWS_REGION=ap-northeast-1
```

### Testing
```bash
python -m pytest tests/                          # All tests
python -m pytest tests/ -k "not integration"     # Unit tests only
python -m pytest tests/test_pmda_compliance.py   # PMDA compliance
python -m pytest tests/test_pathogen_detection.py -v  # Specific test with verbose
```

### Pipeline Operations
```bash
./tools/workflow_cli.py start --run-id RUN-2024-001 --bucket minion-data --input-prefix runs/RUN-2024-001/fast5/
./tools/workflow_cli.py status --run-id RUN-2024-001 --watch
./tools/workflow_cli.py metrics --run-id RUN-2024-001
```

### Code Quality
```bash
black scripts/ lambda/ tools/ tests/    # Auto-format
flake8 scripts/ lambda/ tools/ tests/   # Linting
mypy scripts/                           # Type checking
```

### Documentation Portal (Next.js)
```bash
cd docs-portal && npm run dev    # http://localhost:3003
cd docs-portal && npm run build  # Production build
```

## Code Conventions

### Python Module Import Pattern
Phase scripts expect shared libraries in `/opt/minion/lib`:
```python
import sys
sys.path.append('/opt/minion/lib')
from workflow_manager import WorkflowManager
```

### Error Handling Pattern
Always validate files before opening:
```python
if not bam_file.exists():
    raise FileNotFoundError(f"BAM file not found: {bam_file}")
try:
    bam = pysam.AlignmentFile(str(bam_file), "rb")
except Exception as e:
    raise RuntimeError(f"Failed to open BAM file: {e}")
```

### Test Patterns
- Use `moto` for AWS service mocking
- Use absolute paths relative to test file:
```python
test_dir = Path(__file__).parent
repo_root = test_dir.parent
config_file = repo_root / 'templates' / 'config' / 'pmda_pathogens.json'
```

### PERV Detection Implementation
PERV markers hardcoded in `scripts/phase4_pathogen/perv_typing.py`:
```python
PERV_MARKERS = {
    'PERV-A': {'env_start': 5800, 'env_end': 7400,
               'specific_motifs': ['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC']},
    'PERV-B': {'env_start': 5800, 'env_end': 7400,
               'specific_motifs': ['ATGGCAACCACCGTAGC', 'TGGAAACCTGGAAAACC']},
    'PERV-C': {'env_start': 5800, 'env_end': 7400,
               'specific_motifs': ['ATGGCAGCCACCATAGG', 'TGGAGACCTGGAAGAAC']}
}
```

## Configuration Files

### Main Pipeline Configuration
Main config: `templates/config/default_pipeline.yaml`
```yaml
pipeline:
  phases:
    basecalling:
      min_qscore: 9       # Q30 accuracy target
    qc:
      min_reads: 10000000 # 10M minimum
      min_q30: 0.85       # 85% Q30
    pathogen:
      confidence_threshold: 0.1
      pmda_alert_enabled: true
```

### PMDA Pathogen List
PMDA pathogen list: `templates/config/pmda_pathogens.json` (91 pathogens)

## Important Paths
- Reference DBs: `/mnt/efs/databases/pmda/2024.1/`
- Shared libs: `/opt/minion/lib`
- Config templates: `templates/config/`
- Step Functions: `templates/stepfunctions/workflow_definition.json`

## Best Practices
1. Always run tests before committing
2. Use type hints for better code documentation
3. Follow PEP 8 style guide (enforced by black)
4. Document complex algorithms with clear comments
5. Use meaningful variable names
6. Keep functions focused on single responsibilities