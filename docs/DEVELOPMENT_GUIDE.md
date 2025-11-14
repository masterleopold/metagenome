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

### Surveillance System
```bash
# Dashboard (Streamlit)
streamlit run surveillance/dashboard/app.py --server.port 8501

# REST API (FastAPI)
cd surveillance/api && uvicorn main:app --reload --port 8000

# Test external collectors
python -m pytest surveillance/tests/test_external_collectors.py
python surveillance/external/estat_client.py --test
python surveillance/external/maff_scraper.py --test
python surveillance/external/academic_monitor.py --test

# Deploy Lambda functions
cd infrastructure/surveillance && terraform apply

# Manual Lambda trigger (for testing)
aws lambda invoke --function-name surveillance-external-collector --payload '{}' response.json
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

### 4-Virus Surveillance Implementation
Target virus configuration in `surveillance/internal/pipeline_listener.py`:
```python
TARGET_VIRUSES = {
    'hantavirus': {
        'taxa_names': ['Hantaan virus', 'Hantavirus', 'Seoul virus'],
        'ncbi_taxids': [1980519, 11594, 11596, 1980416]
    },
    'polyomavirus': {
        'taxa_names': ['Sus scrofa polyomavirus', 'Polyomavirus'],
        'ncbi_taxids': [1891763]
    },
    'spumavirus': {
        'taxa_names': ['Porcine type-C oncovirus', 'Spumavirus'],
        'ncbi_taxids': [35268]
    },
    'eeev': {
        'taxa_names': ['Eastern equine encephalitis virus', 'EEEV'],
        'ncbi_taxids': [11021, 11019]
    }
}
```

Severity classification in `surveillance/config/severity_rules.yaml`:
```yaml
virus_rules:
  spumavirus:
    severity_thresholds:
      - condition:
          source: internal_pipeline
          copies_per_ml: "> 500"
        severity: critical
        reason: "High viral load - xenotransplantation risk"
        notification_channels: ["sns", "email", "sms", "dashboard"]
```

External source keyword detection pattern:
```python
VIRUS_KEYWORDS = {
    'hantavirus': ['ハンタウイルス', 'ハンタ', 'hantavirus', 'HANTAV'],
    'polyomavirus': ['ポリオーマウイルス', 'ポリオーマウィルス', 'ポリオーマ'],
    'spumavirus': ['スピューマウイルス', 'スピューマウィルス', 'スピューマ'],
    'eeev': ['東部ウマ脳炎', 'EEEV', 'Eastern equine encephalitis']
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
- Surveillance config: `surveillance/config/config.yaml`, `surveillance/config/severity_rules.yaml`
- Surveillance data (S3): `s3://surveillance-data/external/`, `s3://surveillance-data/internal/`
- DynamoDB tables: `surveillance-detections`, `surveillance-external-updates`, `surveillance-notifications`

## Best Practices
1. Always run tests before committing
2. Use type hints for better code documentation
3. Follow PEP 8 style guide (enforced by black)
4. Document complex algorithms with clear comments
5. Use meaningful variable names
6. Keep functions focused on single responsibilities