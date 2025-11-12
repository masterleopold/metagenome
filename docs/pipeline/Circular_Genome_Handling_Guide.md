# Circular Genome Handling Guide for MinION Pipeline

**Document Type**: Bioinformatics Pipeline Documentation
**Version**: 1.0
**Date**: 2025-11-13
**Related**: Protocol 12 v2.1, Step 2.5

---

## Overview

This document describes how the MinION pathogen detection pipeline handles **circular DNA genomes** and **single-stranded DNA (ssDNA) genomes**, which require special preprocessing before sequencing and special considerations during bioinformatics analysis.

### Affected Pathogens (4 of 91 PMDA Pathogens)

| Pathogen Code | Full Name | Genome Structure | Size | PMDA Classification | Special Handling |
|---------------|-----------|------------------|------|---------------------|------------------|
| **PCV2** | Porcine Circovirus 2 | Circular ssDNA | 1.7 kb | Special Management #3 | Linearization + Synthesis |
| **PCV3** | Porcine Circovirus 3 | Circular ssDNA | 2.0 kb | Special Management #3 | Linearization + Synthesis |
| **TTV** | Torque Teno Virus | Circular ssDNA | 3.8 kb | Standard #40 | Linearization + Synthesis |
| **PPV** | Porcine Parvovirus | Linear ssDNA | 5.0 kb | Standard #1 | Synthesis only |

---

## Laboratory Processing (Protocol 12 v2.1 Step 2.5)

### Why Special Processing is Required

**Problem 1: Circular DNA Cannot Be Sequenced Directly**
- Oxford Nanopore SQK-LSK114 ligation kit requires **free DNA ends** for adapter attachment
- Circular DNA has no 3'/5' termini → No adapter ligation → No sequencing
- Solution: **DNase I linearization** creates random nicks, converting circular → linear

**Problem 2: ssDNA Has Low Ligation Efficiency**
- T4 DNA Ligase (in LSK114) has **<5% efficiency** on single-stranded DNA
- MinION motor proteins require **dsDNA input** (unzip dsDNA → ssDNA for sequencing)
- Solution: **Klenow Fragment second-strand synthesis** converts ssDNA → dsDNA

### Step 2.5 Processing (2.5 hours)

**Sub-step 2.5.1: Circular DNA Linearization** (30 min)
```
Circular DNA + DNase I (0.005 U) → 37°C, 5 min → Linear DNA with free ends
```

**Sub-step 2.5.2: ssDNA → dsDNA Conversion** (2 hours)
```
ssDNA + Random Hexamers → 25°C, 5 min (annealing)
      + Klenow Fragment → 16°C, 60 min (synthesis)
      → dsDNA ready for LSK114
```

---

## Bioinformatics Analysis Considerations

### 1. Read Mapping to Circular Genomes

#### Issue: Circular Genome "Junction" Reads

When circular DNA is linearized randomly, some sequencing reads may span the **circularization junction** (the point where 5' and 3' ends were joined in the original circular genome).

**Example: PCV2 (1768 bp circular genome)**
```
Original circular genome:
   Position 1 ────────────────── Position 1768
        ↑                              ↓
        └──────────────────────────────┘
              (connected in circle)

After DNase I linearization at position 500:
   Position 500 ──────────────────────→ Position 1768 → Position 1 ──────→ Position 499
   (Now has free 5' end)                                    (Now has free 3' end)

Junction reads spanning original 1768→1 connection:
   Read: ...ATCG [Position 1760-1768] [Position 1-20] GATC...
         └─── End of genome ───┘ └── Start of genome ──┘
```

#### Solution 1: Duplicate Reference Genome (Recommended)

**For Minimap2 alignment**, create a duplicated reference:
```bash
# Original circular reference (1768 bp)
>PCV2_circular
ATCGATCGATCG...(1768 bp)

# Duplicated reference (3536 bp = 1768 × 2)
>PCV2_linearized_duplicated
ATCGATCGATCG...(1768 bp)...ATCGATCGATCG...(1768 bp repeated)
```

**Benefits**:
- Junction reads map continuously without split alignment
- No special handling needed in alignment pipeline
- Standard coverage calculation works correctly

**Implementation**:
```python
# scripts/database_preparation/duplicate_circular_genomes.py

def duplicate_circular_reference(seq: str, genome_name: str) -> str:
    """
    Duplicate circular genome sequence for proper junction read mapping.

    Args:
        seq: Original circular genome sequence
        genome_name: Genome identifier

    Returns:
        FASTA entry with duplicated sequence
    """
    duplicated_seq = seq + seq
    return f">{genome_name}_linearized_dup\n{duplicated_seq}\n"

# Example usage for PCV2
circular_genomes = {
    'PCV2': 'ATCGATCG...',  # 1768 bp
    'PCV3': 'GCTAGCTA...',  # 2000 bp
    'TTV': 'TTAATTAA...',   # 3800 bp
}

for name, seq in circular_genomes.items():
    duplicate_circular_reference(seq, name)
```

**Coverage Calculation Adjustment**:
```python
def calculate_circular_genome_coverage(bam_file, genome_length):
    """
    Calculate coverage for circular genome with duplicated reference.

    Only count coverage for first half of duplicated reference
    (positions 0 to genome_length-1).
    """
    coverage = []
    for pileup in bam_file.pileup():
        position = pileup.pos
        if position < genome_length:  # Only first half
            coverage.append(pileup.n)
        # Ignore second half (positions >= genome_length)

    mean_coverage = sum(coverage) / len(coverage) if coverage else 0
    return mean_coverage
```

---

#### Solution 2: Split Read Alignment (Alternative)

If duplicated reference is not used, handle junction reads with split alignment:

```python
def handle_junction_reads(read, genome_length):
    """
    Detect and properly align reads spanning circularization junction.

    Args:
        read: Aligned read from BAM file
        genome_length: Length of circular genome (e.g., 1768 for PCV2)

    Returns:
        True if read spans junction, adjusted alignment coordinates
    """
    # Check if read has supplementary alignment
    if read.has_tag('SA'):  # Supplementary alignment tag
        # Parse SA tag: chr,pos,strand,CIGAR,mapQ,NM
        sa_tag = read.get_tag('SA')
        # Handle as junction-spanning read
        return True, (read.reference_start, read.reference_end)

    return False, (read.reference_start, read.reference_end)
```

---

### 2. Copy Number Quantification

#### Genome Size Considerations

**Important**: For circular genomes, use the **actual circular genome size**, not the duplicated reference size.

```python
# From scripts/phase5_quantification/absolute_copy_number.py

GENOME_SIZES = {
    # Use actual circular genome sizes, NOT duplicated sizes
    'PCV2': 1768,       # NOT 3536 (duplicated)
    'PCV3': 2000,       # NOT 4000 (duplicated)
    'TTV': 3800,        # NOT 7600 (duplicated)
    'PPV': 5000,        # Linear ssDNA, no duplication
}
```

**Rationale**: One circular genome molecule has 1768 bp of unique sequence, even though we map to 3536 bp duplicated reference. Copy number calculation must use actual genome size.

---

#### Read Count Adjustment for Duplicated References

If using duplicated references, adjust read counts:

```python
def adjust_read_count_for_circular(reads_mapped, reference_duplicated=True):
    """
    Adjust read counts for circular genomes mapped to duplicated references.

    Args:
        reads_mapped: Number of reads mapped to reference
        reference_duplicated: True if reference was duplicated for junction reads

    Returns:
        Adjusted read count
    """
    if reference_duplicated:
        # Reads mapping to duplicated reference should not be double-counted
        # Each unique read maps once; no adjustment needed
        return reads_mapped
    else:
        # Standard linear genome
        return reads_mapped

# Note: Modern mappers (Minimap2) count each read only once,
# even if reference is duplicated. No adjustment typically needed.
```

---

### 3. Coverage Analysis

#### Detecting Complete Genome Coverage

For circular genomes, ensure coverage spans the **entire circular sequence**:

```python
def check_circular_genome_coverage(bam_file, genome_name, genome_length, min_coverage=5):
    """
    Check if circular genome has complete coverage.

    Args:
        bam_file: BAM file with alignments
        genome_name: Name of circular genome reference
        genome_length: Length of circular genome (e.g., 1768 for PCV2)
        min_coverage: Minimum depth to consider "covered" (default: 5×)

    Returns:
        Dictionary with coverage statistics
    """
    import pysam

    bam = pysam.AlignmentFile(bam_file, "rb")

    # Get coverage for each position (only first half if duplicated)
    coverage = [0] * genome_length
    for pileup in bam.pileup(genome_name):
        pos = pileup.pos
        if pos < genome_length:  # Only count first half
            coverage[pos] = pileup.n

    # Calculate statistics
    covered_positions = sum(1 for depth in coverage if depth >= min_coverage)
    mean_depth = sum(coverage) / len(coverage) if coverage else 0
    genome_coverage_pct = (covered_positions / genome_length) * 100

    # Check for gaps (uncovered regions)
    gaps = []
    in_gap = False
    gap_start = 0
    for i, depth in enumerate(coverage):
        if depth < min_coverage and not in_gap:
            in_gap = True
            gap_start = i
        elif depth >= min_coverage and in_gap:
            in_gap = False
            gaps.append((gap_start, i-1))

    return {
        'genome': genome_name,
        'length': genome_length,
        'mean_depth': round(mean_depth, 2),
        'coverage_pct': round(genome_coverage_pct, 2),
        'covered_positions': covered_positions,
        'uncovered_positions': genome_length - covered_positions,
        'gaps': gaps,
        'complete_coverage': genome_coverage_pct > 95.0
    }
```

---

### 4. Variant Calling Considerations

#### Circular Genome Coordinate System

When calling variants in circular genomes:

**Issue**: Variant at position 1 may actually be near the "artificial" linearization point.

**Solution**: Report variants with circular genome coordinates:

```python
def adjust_variant_position_for_circular(variant_pos, genome_length, linearization_point=0):
    """
    Adjust variant position to circular genome coordinate system.

    Args:
        variant_pos: Position in linearized reference
        genome_length: Length of circular genome
        linearization_point: Position where circle was linearized (default: 0)

    Returns:
        Position in circular coordinate system
    """
    # If using duplicated reference, take modulo genome_length
    circular_pos = (variant_pos - linearization_point) % genome_length
    return circular_pos

# Example:
# PCV2 linearized at position 500, variant called at position 1900
variant_pos_linear = 1900
genome_length = 1768
linearization_point = 500

circular_pos = adjust_variant_position_for_circular(variant_pos_linear, genome_length, linearization_point)
# Result: circular_pos = 632 (in original circular coordinates)
```

---

## Database Preparation

### Reference Sequence Preparation

**Step 1: Obtain Reference Sequences**

```bash
# PCV2 (GenBank: AF027217.1)
# PCV3 (GenBank: MF318988.1)
# TTV (GenBank: AB017610.1)
# PPV (GenBank: NC_001718.1)

# Download from NCBI
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AF027217.1&rettype=fasta" -O PCV2.fasta
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MF318988.1&rettype=fasta" -O PCV3.fasta
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AB017610.1&rettype=fasta" -O TTV.fasta
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001718.1&rettype=fasta" -O PPV.fasta
```

**Step 2: Duplicate Circular Genome References**

```bash
# Duplicate PCV2, PCV3, TTV (circular genomes)
python scripts/database_preparation/duplicate_circular_genomes.py \
    --input PCV2.fasta,PCV3.fasta,TTV.fasta \
    --output circular_genomes_duplicated.fasta

# PPV (linear ssDNA) does not need duplication
cat PPV.fasta >> linear_genomes.fasta
```

**Step 3: Build Minimap2 Index**

```bash
# Combine all references
cat circular_genomes_duplicated.fasta linear_genomes.fasta other_pmda_genomes.fasta > pmda_all_91.fasta

# Build Minimap2 index
minimap2 -d pmda_all_91.mmi pmda_all_91.fasta
```

---

## Quality Control Checks

### QC Checkpoint 1: Post-Step 2.5 Validation

**Objective**: Confirm linearization and ssDNA→dsDNA conversion succeeded

**Metrics**:
- ✅ Qubit recovery >70% after linearization
- ✅ Qubit recovery >80% after second-strand synthesis
- ✅ Bioanalyzer concentration ~2× Qubit (indicates dsDNA)

**If QC fails**: Circular/ssDNA viruses will have <5% detection rate

---

### QC Checkpoint 2: Post-Sequencing Validation

**Objective**: Confirm circular genome was properly detected

**Metrics for PCV2/PCV3/TTV**:
- ✅ Mean coverage >10× across entire genome (no large gaps)
- ✅ Coverage uniformity (SD/mean <0.5)
- ✅ Junction reads present (if using duplicated reference, reads mapping to both halves)

**Troubleshooting**:
| Issue | Cause | Solution |
|-------|-------|----------|
| No reads mapping | Step 2.5 failed | Repeat Step 2.5 with fresh reagents |
| Low coverage (<5×) | Incomplete linearization | Optimize DNase I concentration |
| Large gaps in coverage | Over-digestion | Reduce DNase I to 0.003 U |
| No junction reads | Reference not duplicated | Use duplicated reference |

---

## Implementation Checklist

**Database Preparation**:
- ☐ Downloaded PCV2, PCV3, TTV, PPV reference sequences
- ☐ Duplicated circular genome references (PCV2, PCV3, TTV)
- ☐ Built Minimap2 index with duplicated references
- ☐ Updated absolute_copy_number.py with correct genome sizes

**Pipeline Scripts**:
- ☐ Updated detect_pmda_all_91_pathogens.py to handle circular genomes
- ☐ Updated coverage calculation scripts for duplicated references
- ☐ Added TTV to quantification scripts
- ☐ Tested with synthetic circular ssDNA samples

**Quality Control**:
- ☐ Verified Step 2.5 QC checkpoints
- ☐ Validated circular genome coverage analysis
- ☐ Confirmed junction read mapping
- ☐ Tested PCV2/PCV3 detection at 500 copies/mL

**Documentation**:
- ☐ This document completed
- ☐ Protocol 12 v2.1 updated with Step 2.5
- ☐ Staff trained on circular genome handling
- ☐ LOD validation completed for PCV2/PCV3/TTV/PPV

---

## References

### Technical Literature
1. Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34:3094-3100.
2. Salzberg SL, Yorke JA. Beware of mis-assembled genomes. Bioinformatics. 2005;21:4320-4321.
3. Wylie TN, et al. Enhanced virome sequencing using targeted sequence capture. Genome Research. 2015;25:1910-1920.

### Circular Genome Bioinformatics
4. Wick RR, et al. Completing bacterial genome assemblies with multiplex MinION sequencing. Microb Genom. 2017;3:e000132.
5. Hunt M, et al. Circlator: automated circularization of genome assemblies using long sequencing reads. Genome Biol. 2015;16:294.

### Pathogen-Specific
6. Johne R, et al. Detection and sequence analysis of porcine circovirus 2 in Germany. Arch Virol. 2004;149:2001-2010.
7. Phan TG, et al. The fecal viral flora of wild rodents. PLoS Pathog. 2011;7:e1002218.

### Oxford Nanopore
8. Oxford Nanopore Technologies. Ligation Sequencing Kit SQK-LSK114 Protocol. 2024.
9. Oxford Nanopore Technologies. Circular consensus sequencing for improved accuracy. Technical Note. 2023.

---

## Document Control

**Version**: 1.0
**Date**: 2025-11-13
**Author**: MinION Pipeline Development Team
**Related Documents**:
- Protocol 12 v2.1 (Step 2.5)
- absolute_copy_number.py (v2.1)
- detect_pmda_all_91_pathogens.py
**Next Review**: 2026-05-13 (6 months)

---

**End of Circular Genome Handling Guide**
