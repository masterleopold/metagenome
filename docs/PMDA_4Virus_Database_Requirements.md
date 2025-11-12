# PMDA 4-Virus Database Requirements

**Version:** 1.0
**Date:** 2025-11-12
**Purpose:** Reference genome database requirements for high-sensitivity detection of 4 PMDA-required viruses

---

## Overview

This document specifies the reference genome databases required for detecting:
1. Polyomavirus (dsDNA)
2. Hantavirus (ssRNA-, trisegmented)
3. Eastern Equine Encephalitis Virus - EEEV (ssRNA+)
4. Porcine Spumavirus (retrovirus, proviral DNA)

All databases should be stored in: `/mnt/efs/databases/pmda/2024.1/`

---

## 1. Polyomavirus Database

### 1.1 Directory Structure
```
/mnt/efs/databases/pmda/2024.1/polyomavirus/
├── polyoma_all.fasta           # All polyomavirus genomes
├── polyoma_all.mmi             # Minimap2 index
├── polyoma_kraken2/            # Kraken2 database
│   ├── hash.k2d
│   ├── opts.k2d
│   └── taxo.k2d
├── polyoma_blast/              # BLAST database
│   ├── polyoma.nhr
│   ├── polyoma.nin
│   └── polyoma.nsq
└── metadata.tsv                # Species, accessions, lineages
```

### 1.2 Required Reference Genomes

| Species | Accession | Length (bp) | Notes |
|---------|-----------|-------------|-------|
| BK polyomavirus | NC_001538 | 5,153 | Human polyomavirus, reference |
| JC polyomavirus | NC_001699 | 5,130 | Human polyomavirus |
| SV40 | NC_001669 | 5,243 | Simian virus 40 |
| Sus scrofa polyomavirus 2 | MH381769 | 5,185 | **Critical: Porcine-specific** |
| Merkel cell polyomavirus | NC_010277 | 5,387 | Cross-reference |
| Karolinska Institutet polyomavirus | NC_018102 | 5,024 | Cross-reference |

### 1.3 Database Construction

```bash
# Download sequences from NCBI
esearch -db nucleotide -query "Polyomaviridae[Organism] AND complete genome" | \
  efetch -format fasta > polyoma_all.fasta

# Build Minimap2 index
minimap2 -d polyoma_all.mmi polyoma_all.fasta

# Build Kraken2 database
kraken2-build --download-taxonomy --db polyoma_kraken2/
kraken2-build --add-to-library polyoma_all.fasta --db polyoma_kraken2/
kraken2-build --build --db polyoma_kraken2/

# Build BLAST database
makeblastdb -in polyoma_all.fasta -dbtype nucl -out polyoma_blast/polyoma
```

### 1.4 Detection Parameters

- **Mapping tool:** Minimap2 (`-ax map-ont`)
- **Minimum read length:** 100 bp
- **Mapping quality threshold:** Q10
- **Confidence threshold:** ≥100 reads AND ≥10× coverage
- **Cross-reactivity check:** Distinguish from PERV, porcine parvovirus

---

## 2. Hantavirus Database

### 2.1 Directory Structure
```
/mnt/efs/databases/pmda/2024.1/hantavirus/
├── hantavirus_L_segment.fasta
├── hantavirus_M_segment.fasta
├── hantavirus_S_segment.fasta
├── hantavirus_all.fasta        # Concatenated L+M+S
├── hantavirus_all.mmi
├── hantavirus_kraken2/
├── amplicon_primers/
│   ├── hantavirus_pool_A.bed   # Primer Pool A coordinates
│   ├── hantavirus_pool_B.bed   # Primer Pool B coordinates
│   └── primer_sequences.tsv
└── metadata.tsv
```

### 2.2 Required Reference Genomes

| Species | Segment | Accession | Length (nt) | Geographic Region |
|---------|---------|-----------|-------------|-------------------|
| **Hantaan virus** | L | NC_005222 | 6,533 | Asia (Korea, China, Japan) |
| | M | NC_005219 | 3,651 | |
| | S | NC_005218 | 1,696 | |
| **Seoul virus** | L | NC_005238 | 6,544 | Worldwide (urban rodents) |
| | M | NC_005236 | 3,651 | |
| | S | NC_005237 | 1,715 | |
| **Dobrava-Belgrade virus** | L | NC_005233 | 6,550 | Europe |
| | M | NC_005234 | 3,643 | |
| | S | NC_005235 | 1,700 | |
| **Puumala virus** | L | NC_005224 | 6,530 | Europe |
| | M | NC_005223 | 3,623 | |
| | S | NC_005222 | 1,764 | |

### 2.3 Amplicon Primer Design (for high sensitivity)

**Based on:** Kim et al., Viruses 2021 (Hantaan virus multiplex PCR)

```
Total amplicons: 36 (L: 20, M: 11, S: 5)
Amplicon size: 400 bp
Overlap: 100 bp
Primer design tool: Primal Scheme (ARTIC Network)
```

**Example Primer Pool A (18 amplicons):**
```tsv
Name	Sequence	Length	Tm	Segment	Start	End
HTNV_L_1_LEFT	AGCTGACACATCAGTGTGCC	20	60.2	L	1	400
HTNV_L_1_RIGHT	GGTCAACAGTACCTGCCTAG	20	60.1	L	380	780
HTNV_L_3_LEFT	TGGACCTGATACCAAGGTCC	20	60.3	L	701	1100
...
```

### 2.4 Detection Parameters

- **3-segment concordance required:** L AND M AND S all detected
- **Minimum reads per segment:** 50 reads
- **Mapping tool:** Minimap2 or Kraken2
- **Amplicon mapping:** Align to expected amplicon coordinates
- **Segmented genome assembly:** Use medaka or Nanopolish for consensus

---

## 3. Eastern Equine Encephalitis Virus (EEEV) / Alphavirus Database

### 3.1 Directory Structure
```
/mnt/efs/databases/pmda/2024.1/alphavirus/
├── alphavirus_all.fasta        # All PMDA-listed alphaviruses
├── alphavirus_all.mmi
├── eeev_references.fasta       # EEEV-specific for phylogeny
├── eeev_phylogeny/
│   ├── eeev_alignment.fasta
│   └── eeev_tree.nwk
├── alphavirus_kraken2/
└── metadata.tsv
```

### 3.2 Required Reference Genomes

**PMDA-listed Alphaviruses (厚労省異種移植指針):**

| Virus | PMDA Line | Accession | Length (nt) | Poly(A)+ |
|-------|-----------|-----------|-------------|----------|
| **Eastern equine encephalitis virus (EEEV)** | 32/55 | NC_003899 | 11,841 | Yes |
| **Western equine encephalitis virus (WEEV)** | 32/55 | NC_003908 | 11,722 | Yes |
| **Venezuelan equine encephalitis virus (VEEV)** | 33/56 | NC_001449 | 11,441 | Yes |
| **Getah virus** | 18/41 | NC_003696 | 11,680 | Yes |
| Chikungunya virus | Reference | NC_004162 | 11,805 | Yes |
| Ross River virus | Reference | NC_001544 | 11,659 | Yes |

### 3.3 EEEV Lineage References

For phylogenetic lineage assignment (North American vs South American):

| Lineage | Representative Strain | Accession | Geographic Origin |
|---------|----------------------|-----------|-------------------|
| North American | EEEV FL93-939 | NC_003899 | USA (Florida) |
| South American | EEEV BeAn 436284 | KP765787 | Brazil |
| South American | EEEV Argentina | FJ402886 | Argentina |

### 3.4 Detection Parameters

- **Mapping tool:** Minimap2
- **Minimum reads:** ≥100 reads
- **Coverage:** ≥10× mean depth
- **Species assignment:** Best-hit reference
- **Phylogenetic confirmation:** Align nsP1 or E1 gene → FastTree
- **Cross-reactivity:** Distinguish between EEEV, WEEV, VEEV, Getah

---

## 4. Porcine Spumavirus Database

### 4.1 Directory Structure
```
/mnt/efs/databases/pmda/2024.1/spumavirus/
├── sfv_pol_gene.fasta          # Simian foamy virus pol genes
├── ffv_pol_gene.fasta          # Feline foamy virus pol genes
├── bfv_pol_gene.fasta          # Bovine foamy virus pol genes
├── spumavirus_all_pol.fasta    # All spumavirus pol genes
├── spumavirus_pol.mmi
├── perv_pol_gene.fasta         # For discrimination
├── spumavirus_blast/           # BLAST for sequence confirmation
└── metadata.tsv
```

### 4.2 Required Reference Sequences

**IMPORTANT:** No porcine spumavirus reference genome exists in NCBI. Use cross-genus references.

| Species | Gene | Accession | Length (bp) | Identity to pig (predicted) |
|---------|------|-----------|-------------|----------------------------|
| Simian foamy virus (SFV) | pol | NC_001364 (nt 2000-5000) | ~3,000 | 30-50% (degenerate) |
| Feline foamy virus (FFV) | pol | NC_001871 (nt 2100-5100) | ~3,000 | 30-50% |
| Bovine foamy virus (BFV) | pol | NC_001831 (nt 2050-5050) | ~3,000 | 30-50% |
| **PERV (for discrimination)** | pol | AF038600 (nt 4500-6500) | ~2,000 | **Porcine endogenous** |

### 4.3 Nested PCR Primer Sequences

**Degenerate primers targeting foamy virus pol gene (reverse transcriptase domain):**

```
# Outer primers (1st PCR)
FV-pol-F1: 5'-GGNCARATHGGNATGTTYGG-3' (degeneracy: 96-fold)
FV-pol-R1: 5'-CCRTCNCCRAANCCRTC-3' (degeneracy: 64-fold)
Expected product: ~800 bp

# Inner primers (2nd PCR)
FV-pol-F2: 5'-ATHGGNCARGGNTTYACNAC-3' (degeneracy: 96-fold)
FV-pol-R2: 5'-GTRTCNGTYTTRTCNCC-3' (degeneracy: 64-fold)
Expected product: ~400 bp

# Degeneracy code
N = A/T/G/C (4-fold)
R = A/G (2-fold, purine)
Y = C/T (2-fold, pyrimidine)
H = A/T/C (3-fold, not G)
```

### 4.4 Detection and Confirmation Strategy

1. **Nested PCR** from PBMC genomic DNA
   - Target: pol gene (most conserved region)
   - LOD: 1-10 copies/10⁵ PBMCs

2. **Agarose gel confirmation**
   - Expected band: ~400 bp
   - Negative control: Spumavirus-negative PBMC DNA
   - Positive control: SFV-positive primate DNA (if available)

3. **Sanger sequencing (MANDATORY)**
   - Sequence both strands
   - BLAST against NCBI nt database
   - Confirm identity: >70% to Spumavirus, <50% to PERV

4. **Phylogenetic analysis**
   - Align to SFV/FFV/BFV/PERV pol genes
   - Confirm clustering with Spumaretrovirinae (not Gammaretrovirus/PERV)

---

## 5. Quality Control and Validation

### 5.1 Database Validation Checklist

- [ ] All reference sequences downloaded and verified (MD5 checksums)
- [ ] No duplicate sequences (deduplicate by accession)
- [ ] Taxonomy IDs correctly assigned (NCBI Taxonomy)
- [ ] Minimap2 indices built and tested
- [ ] Kraken2 databases built and tested
- [ ] BLAST databases formatted and tested
- [ ] Metadata files complete (species, accession, length, notes)

### 5.2 Detection Validation (In Silico)

**Positive Controls (synthetic reads):**
- Generate simulated MinION reads (Badread or NanoSim)
- Spike into negative pig plasma metagenome
- Concentrations: 10, 50, 100, 500, 1,000 copies/mL
- Validate detection at ≥50 copies/mL

**Negative Controls:**
- Pathogen-free pig plasma metagenome
- No false positive detections

**Cross-Reactivity Testing:**
- Polyomavirus vs porcine parvovirus, PERV
- Hantavirus vs other bunyaviruses
- EEEV vs other alphaviruses (Getah, VEEV, WEEV)
- Spumavirus vs PERV (critical discrimination)

### 5.3 Database Update Schedule

- **Quarterly review:** Check NCBI for new sequences
- **Annual rebuild:** Rebuild all databases with latest references
- **Triggered updates:** If novel variant detected in clinical samples

---

## 6. Storage and Backup

### 6.1 Primary Storage
- **Location:** `/mnt/efs/databases/pmda/2024.1/`
- **File system:** AWS EFS (elastic, shared across EC2 instances)
- **Access:** Read-only for analysis EC2 instances
- **Size estimate:**
  - Polyomavirus: 50 MB
  - Hantavirus: 100 MB
  - Alphavirus: 200 MB
  - Spumavirus: 50 MB
  - **Total:** ~400 MB (plus Kraken2 indices: ~5 GB)

### 6.2 Backup Strategy
- **S3 backup:** `s3://minion-data/databases/pmda/`
- **Backup frequency:** Weekly automated sync
- **Versioning:** Enabled (S3 versioning)
- **Disaster recovery:** Restore from S3 to new EFS mount

---

## 7. Database Construction Scripts

### 7.1 Master Build Script

```bash
#!/bin/bash
# build_pmda_4virus_databases.sh

set -e

BASE_DIR="/mnt/efs/databases/pmda/2024.1"
mkdir -p $BASE_DIR

# 1. Polyomavirus
echo "Building polyomavirus database..."
bash scripts/databases/build_polyomavirus_db.sh $BASE_DIR/polyomavirus

# 2. Hantavirus
echo "Building hantavirus database..."
bash scripts/databases/build_hantavirus_db.sh $BASE_DIR/hantavirus

# 3. Alphavirus (EEEV)
echo "Building alphavirus database..."
bash scripts/databases/build_alphavirus_db.sh $BASE_DIR/alphavirus

# 4. Spumavirus
echo "Building spumavirus database..."
bash scripts/databases/build_spumavirus_db.sh $BASE_DIR/spumavirus

# 5. Sync to S3
echo "Backing up to S3..."
aws s3 sync $BASE_DIR s3://minion-data/databases/pmda/2024.1/ --delete

echo "Database build complete!"
```

---

## 8. Integration with Pipeline

### 8.1 Phase 4 (Pathogen Detection) Integration

**Current Phase 4 databases:**
- `/mnt/efs/databases/kraken2_standard/`
- `/mnt/efs/databases/blast_nt/`
- `/mnt/efs/databases/perv_references/`

**New Phase 4 additions:**
```python
# scripts/phase4_pathogen/detect_pmda_4viruses.py

PMDA_4VIRUS_DATABASES = {
    "polyomavirus": "/mnt/efs/databases/pmda/2024.1/polyomavirus/polyoma_all.mmi",
    "hantavirus": "/mnt/efs/databases/pmda/2024.1/hantavirus/hantavirus_all.mmi",
    "alphavirus": "/mnt/efs/databases/pmda/2024.1/alphavirus/alphavirus_all.mmi",
    "spumavirus": "/mnt/efs/databases/pmda/2024.1/spumavirus/spumavirus_all_pol.fasta"  # BLAST only
}
```

---

## 9. References

1. NCBI Nucleotide Database: https://www.ncbi.nlm.nih.gov/nuccore/
2. NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy/
3. Minimap2: https://github.com/lh3/minimap2
4. Kraken2: https://github.com/DerrickWood/kraken2
5. BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi
6. Protocol 11: MinION_Protocol_11_PMDA_4ウイルス高感度検出プロトコル.md

---

**Document Control:**
- Author: Claude Code
- Version: 1.0
- Date: 2025-11-12
- Next Review: 2025-12-12 (monthly)
- Approval: Pending
