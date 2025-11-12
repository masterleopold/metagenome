# Phase 0: Sample Preparation

Sample preparation router for PMDA 4-virus high-sensitivity detection.

## Overview

Phase 0 determines the optimal nucleic acid extraction and library preparation workflow based on target viruses:

- **Polyomavirus** (dsDNA) → cfDNA extraction → CpG depletion → DNA library
- **Hantavirus** (ssRNA-, no poly(A)) → cfRNA extraction → rRNA depletion → Amplicon RT-PCR
- **EEEV** (ssRNA+, poly(A)+) → cfRNA extraction → Poly(A) selection → Direct RNA
- **Spumavirus** (proviral DNA) → PBMC genomic DNA → Nested PCR

## Usage

### Basic Usage

```bash
# Single virus
python sample_router.py \
  --targets polyomavirus \
  --metadata sample_001.json \
  --output workflows/sample_001/

# Multiple viruses
python sample_router.py \
  --targets polyomavirus hantavirus eeev \
  --metadata sample_001.json \
  --output workflows/sample_001/
```

### Validation Only

```bash
# Check if sample meets requirements without generating workflows
python sample_router.py \
  --targets polyomavirus hantavirus eeev spumavirus \
  --metadata sample_001.json \
  --output workflows/sample_001/ \
  --validate-only
```

## Sample Metadata Format

See `example_sample_metadata.json` for a complete example.

Required fields:
- `sample_id`: Unique sample identifier
- `blood_volume_ml`: Blood volume collected
- `rnase_inhibitor_added`: Boolean, required for RNA viruses
- `plasma_volume_ml`: Plasma volume after separation
- `pbmc_count`: PBMC count if spumavirus screening required

## Output Files

The script generates:

1. **workflow_summary.json** - Overall workflow with all processing paths
2. **{virus}_extraction_protocol.json** - Detailed protocol for each virus

Example output structure:
```
workflows/sample_001/
├── workflow_summary.json
├── polyomavirus_extraction_protocol.json
├── hantavirus_extraction_protocol.json
├── eeev_extraction_protocol.json
└── spumavirus_extraction_protocol.json
```

## Decision Logic

The router automatically determines:

1. **Blood collection requirements**
   - Volume: 5 mL (DNA only) or 10 mL (RNA viruses)
   - RNase inhibitor: Required for RNA viruses

2. **Sample separation**
   - Plasma: Required for polyomavirus, hantavirus, EEEV
   - PBMCs: Required for spumavirus

3. **Extraction method**
   - Dual DNA/RNA: If both DNA and RNA viruses targeted
   - cfDNA only: DNA viruses only
   - cfRNA only: RNA viruses only
   - PBMC genomic DNA: Spumavirus provirus

4. **Host depletion method**
   - CpG methylation: DNA viruses (polyomavirus)
   - rRNA depletion: RNA viruses without poly(A) (hantavirus)
   - Poly(A) selection: RNA viruses with poly(A) (EEEV)
   - None: Proviral DNA (spumavirus)

5. **Library preparation**
   - DNA ligation: Polyomavirus
   - Direct RNA/cDNA: EEEV
   - Amplicon RT-PCR: Hantavirus (highest sensitivity)
   - Nested PCR: Spumavirus

## Validation Checks

The script validates:
- ✓ Sufficient blood volume
- ✓ RNase inhibitor added for RNA viruses
- ✓ Sufficient plasma volume
- ✓ Sufficient PBMC count (≥1×10⁶)

## Integration with Pipeline

Phase 0 outputs feed into:
- **Phase 1**: Basecalling (after MinION sequencing)
- **Phase 2**: QC (RNA integrity, read quality)
- **Phase 3**: Host removal (method determined by Phase 0)
- **Phase 4**: Pathogen detection (virus-specific algorithms)

## Expected LOD

Based on Protocol 11 validation:
- Polyomavirus: 50-500 copies/mL (with probe capture)
- Hantavirus: 100-500 copies/mL (amplicon approach)
- EEEV: 50-500 copies/mL (poly(A) selection + capture)
- Spumavirus: 1-10 copies/10⁵ PBMCs (nested PCR)

## References

- Protocol 11: PMDA 4ウイルス高感度検出プロトコル
- Appendix D: RNAウイルス検出技術詳解
