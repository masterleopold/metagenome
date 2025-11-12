# Protocols Guide

## PERV Detection (HIGHEST PRIORITY)

- **PERV-A/B/C** subtypes must be detected and quantified
- Immediate SNS alert on ANY PERV detection
- Key scripts: `perv_typing.py`, `detect_recombinants.py`, `perv_phylogenetics.py`
- Must report in copies/mL plasma

### PERV Markers
```python
PERV_MARKERS = {
    'PERV-A': {
        'env_start': 5800, 'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCACAGC', 'TGGAGACCTGGAAGACC']
    },
    'PERV-B': {
        'env_start': 5800, 'env_end': 7400,
        'specific_motifs': ['ATGGCAACCACCGTAGC', 'TGGAAACCTGGAAAACC']
    },
    'PERV-C': {
        'env_start': 5800, 'env_end': 7400,
        'specific_motifs': ['ATGGCAGCCACCATAGG', 'TGGAGACCTGGAAGAAC']
    }
}
```

## PMDA 91 Pathogens

- Full list: `md/厚労省異種移植指針_91病原体リスト.md`
- Requirements: PPA >95%, NPA >98%, R² >0.90
- Database: `/mnt/efs/databases/pmda/2024.1/`
- Categories:
  - Viruses: 41 pathogens
  - Bacteria: 27 pathogens
  - Parasites: 19 pathogens
  - Fungi: 2 pathogens
  - Special Management: 5 pathogens (PCV2, PCV3, PERV-A/B/C)

## Sample Preparation Protocols

### Protocol 12 v2.1 - Unified Workflow (RECOMMENDED)

**TRUE 91/91 pathogen coverage - Primary protocol for all samples**

- **Version**: 2.1 - Now includes circular ssDNA virus support (PCV2/PCV3/TTV/PPV)
- **Workflows**: 2 universal (DNA + RNA) - replaces 3-4 complex variants
- **Time**: 15.5 hours hands-on (includes Step 2.5 for circular/ssDNA handling)
- **Cost**: ¥162,000/sample (+¥5,000 for circular/ssDNA reagents)
- **Coverage**: TRUE 100% (91/91 pathogens) - no exceptions
- **LOD**: 100-500 copies/mL (screening-sufficient, all pathogens)

#### Key Features
- Universal poly(A) selection
- CpG methylation-based host depletion
- **NEW Step 2.5**: Critical for circular/ssDNA viruses
  - Sub-step 2.5.1: Circular DNA linearization (DNase I, ultra-low 0.005 U)
  - Sub-step 2.5.2: ssDNA→dsDNA conversion (Klenow Fragment + Random Hexamers)
- Ensures detection of PCV2/PCV3 (Special Management pathogens), TTV, PPV

**Protocol Document**: `docs/protocols/MinION_Protocol_12_統合サンプル調製プロトコル.md`

### Protocol 11 - High-Sensitivity Enhancement (OPTIONAL)

**Use ONLY when LOD <50 copies/mL required for specific viruses**

- **Target viruses**: Polyomavirus, Hantavirus, EEEV, Spumavirus
- **Method**: Targeted amplification and enrichment
- **Time**: Additional 8 hours
- **Cost**: +¥45,000/sample
- **When to use**: Only after Protocol 12 if ultra-high sensitivity needed

**Protocol Document**: `docs/protocols/MinION_Protocol_11_PMDA_4ウイルス高感度検出プロトコル.md`

### Protocol 13 - Spumavirus-Specific Screening (CONDITIONAL)

**Triggered only when Phase 1 detects retrovirus pol signatures**

- **Trigger rate**: 5-10% of samples (mostly PERV false positives)
- **Detection probability**: <0.1% (0 detections in 70 years globally)
- **Rationale**: PMDA precautionary testing despite zero historical detections
- **Method**: Nested PCR + phylogenetic analysis + PERV discrimination
- **Time**: 12 hours
- **Cost**: ¥35,000/sample

#### Scientific Context
- **Historical data**: 0 detections in pigs (1954-2025)
- **NCBI sequences**: 0 porcine spumavirus sequences
- **Publications**: 0 peer-reviewed papers confirming detection
- **PMDA inclusion reason**: Phylogenetic proximity to other species with foamy viruses
- **Impact if detected**: Nature/Science-level discovery, immediate PMDA notification

**Protocol Document**: `docs/protocols/MinION_Protocol_13_スピューマウイルス専用検査プロトコル.md`

## Protocol Selection Flowchart

```
Start → Run Protocol 12 v2.1 (Universal)
         ↓
    [Phase 1-4 Analysis]
         ↓
    Retrovirus pol detected?
         ├─ No (90-95%) → Continue standard pipeline
         └─ Yes (5-10%) → Trigger Protocol 13
                           ↓
                      [PERV discrimination]
                           ↓
                      Is it PERV?
                      ├─ Yes (>99%) → Report PERV detection
                      └─ No (<1%) → Potential spumavirus
                                    (Scientific discovery)
```

## Critical Circular/ssDNA Virus Information

### Problem Discovered (2025-11-13)
Standard SQK-LSK114 ligation-based library prep CANNOT detect:
1. **Circular DNA viruses**: No free ends for adapter ligation
2. **Single-stranded DNA viruses**: T4 Ligase has <5% efficiency on ssDNA

### Affected Pathogens
- **PCV2** (Porcine Circovirus 2) - Special Management pathogen
- **PCV3** (Porcine Circovirus 3) - Special Management pathogen
- **TTV** (Torque Teno Virus)
- **PPV** (Porcine Parvovirus)

### Solution (Protocol 12 v2.1 Step 2.5)
1. **Linearization**: DNase I at ultra-low concentration (0.005 U) creates single nick in circular DNA
2. **Conversion**: Klenow Fragment converts ssDNA to dsDNA using random hexamers
3. **Result**: All DNA now compatible with standard library prep

## Quality Control Requirements

### Minimum Standards
- **Reads**: ≥10M total reads
- **Q30 percentage**: ≥85%
- **Host depletion**: >90% efficiency
- **RIN score** (RNA samples): ≥7.0

### Critical Thresholds
- **PERV detection**: ANY detection triggers immediate alert
- **Special Management pathogens**: Require confirmatory testing
- **Cost alert**: Analysis exceeding $400

## Documentation References

- **Japanese protocols**: `docs/protocols/` directory (18 protocol files)
- **PMDA compliance**: `docs/PMDA_Simplified_Sample_Prep_Strategy.md`
- **Workflow diagrams**: `docs/PMDA_Simplified_Workflow_Flowchart.md`
- **Technical details**: `docs/TECHNICAL_DETAILS.md`
- **Session logs**: `docs/claude-sessions/` - Development session documentation