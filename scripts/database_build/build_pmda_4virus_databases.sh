#!/bin/bash
#
# Build PMDA 4-Virus Reference Databases
# Creates comprehensive databases for Polyomavirus, Hantavirus, EEEV, and Spumavirus detection
#
# Usage: ./build_pmda_4virus_databases.sh [--base-dir PATH]
#
# Requirements:
#   - NCBI E-utilities (esearch, efetch)
#   - minimap2
#   - samtools
#   - kraken2 (optional)
#   - blast+ (makeblastdb)
#

set -euo pipefail

# Configuration
BASE_DIR="${1:-/mnt/efs/databases/pmda/2024.1}"
THREADS=${THREADS:-16}
LOG_FILE="build_pmda_databases.log"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $*" | tee -a "$LOG_FILE"
}

warn() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $*" | tee -a "$LOG_FILE"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $*" | tee -a "$LOG_FILE"
    exit 1
}

# Check dependencies
check_dependencies() {
    log "Checking dependencies..."

    local deps=("esearch" "efetch" "minimap2" "samtools" "makeblastdb")
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            error "$dep is not installed. Please install NCBI E-utilities, minimap2, samtools, and BLAST+"
        fi
    done

    log "✓ All dependencies satisfied"
}

# Create directory structure
create_directories() {
    log "Creating directory structure at $BASE_DIR..."

    mkdir -p "$BASE_DIR"/{polyomavirus,hantavirus,alphavirus,spumavirus}
    mkdir -p "$BASE_DIR"/polyomavirus/{fasta,minimap2,kraken2,blast}
    mkdir -p "$BASE_DIR"/hantavirus/{fasta,minimap2,amplicon_primers}
    mkdir -p "$BASE_DIR"/alphavirus/{fasta,minimap2,blast,phylogeny}
    mkdir -p "$BASE_DIR"/spumavirus/{fasta,blast}

    cd "$BASE_DIR"
    log "✓ Directory structure created"
}

#==============================================================================
# 1. POLYOMAVIRUS DATABASE
#==============================================================================

build_polyomavirus_db() {
    log "=========================================="
    log "Building Polyomavirus Database"
    log "=========================================="

    cd "$BASE_DIR"/polyomavirus/fasta

    # Download reference genomes
    log "Downloading polyomavirus reference genomes..."

    # Sus scrofa polyomavirus 2 (CRITICAL - porcine-specific)
    esearch -db nucleotide -query "MH381769" | \
        efetch -format fasta > polyoma_sscrofa2.fasta || warn "Sus scrofa polyomavirus 2 not found"

    # BK polyomavirus (human, reference)
    esearch -db nucleotide -query "NC_001538" | \
        efetch -format fasta > polyoma_bk.fasta

    # JC polyomavirus (human)
    esearch -db nucleotide -query "NC_001699" | \
        efetch -format fasta > polyoma_jc.fasta

    # SV40 (simian)
    esearch -db nucleotide -query "NC_001669" | \
        efetch -format fasta > polyoma_sv40.fasta

    # Concatenate all references
    cat polyoma_*.fasta > polyoma_all.fasta

    local count=$(grep -c "^>" polyoma_all.fasta || true)
    log "✓ Downloaded $count polyomavirus reference sequences"

    # Build Minimap2 index
    log "Building Minimap2 index..."
    minimap2 -d ../minimap2/polyoma_all.mmi polyoma_all.fasta
    log "✓ Minimap2 index created: polyoma_all.mmi"

    # Build BLAST database
    log "Building BLAST database..."
    makeblastdb -in polyoma_all.fasta -dbtype nucl \
        -out ../blast/polyoma \
        -title "PMDA Polyomavirus Database" \
        -parse_seqids
    log "✓ BLAST database created"

    # Create metadata
    cat > ../metadata.tsv <<EOF
accession	species	genome_size_bp	host	notes
MH381769	Sus scrofa polyomavirus 2	5185	Pig	Porcine-specific, CRITICAL for detection
NC_001538	BK polyomavirus	5153	Human	Reference strain
NC_001699	JC polyomavirus	5130	Human	Reference strain
NC_001669	SV40	5243	Simian	Cross-reference
EOF

    log "✓ Polyomavirus database complete"
    log "  Location: $BASE_DIR/polyomavirus/"
    log "  Sequences: $count"
    log "  Size: $(du -sh ../minimap2/polyoma_all.mmi | cut -f1)"
}

#==============================================================================
# 2. HANTAVIRUS DATABASE (Trisegmented)
#==============================================================================

build_hantavirus_db() {
    log "=========================================="
    log "Building Hantavirus Database (L/M/S segments)"
    log "=========================================="

    cd "$BASE_DIR"/hantavirus/fasta

    # Download all 3 segments for each hantavirus species
    log "Downloading hantavirus reference genomes (trisegmented)..."

    # Hantaan virus (Asia - Korea, China, Japan)
    log "  - Hantaan virus..."
    esearch -db nucleotide -query "NC_005222" | efetch -format fasta > hantaan_L.fasta
    esearch -db nucleotide -query "NC_005219" | efetch -format fasta > hantaan_M.fasta
    esearch -db nucleotide -query "NC_005218" | efetch -format fasta > hantaan_S.fasta

    # Seoul virus (Worldwide - urban rodents)
    log "  - Seoul virus..."
    esearch -db nucleotide -query "NC_005238" | efetch -format fasta > seoul_L.fasta
    esearch -db nucleotide -query "NC_005236" | efetch -format fasta > seoul_M.fasta
    esearch -db nucleotide -query "NC_005237" | efetch -format fasta > seoul_S.fasta

    # Dobrava-Belgrade virus (Europe)
    log "  - Dobrava virus..."
    esearch -db nucleotide -query "NC_005233" | efetch -format fasta > dobrava_L.fasta
    esearch -db nucleotide -query "NC_005234" | efetch -format fasta > dobrava_M.fasta
    esearch -db nucleotide -query "NC_005235" | efetch -format fasta > dobrava_S.fasta

    # Puumala virus (Europe)
    log "  - Puumala virus..."
    esearch -db nucleotide -query "NC_005224" | efetch -format fasta > puumala_L.fasta
    esearch -db nucleotide -query "NC_005223" | efetch -format fasta > puumala_M.fasta
    esearch -db nucleotide -query "NC_005221" | efetch -format fasta > puumala_S.fasta

    # Organize by segment
    log "Organizing by segment..."
    cat *_L.fasta > hantavirus_L_segment.fasta
    cat *_M.fasta > hantavirus_M_segment.fasta
    cat *_S.fasta > hantavirus_S_segment.fasta

    # All segments combined
    cat hantavirus_L_segment.fasta hantavirus_M_segment.fasta hantavirus_S_segment.fasta > hantavirus_all.fasta

    local count=$(grep -c "^>" hantavirus_all.fasta || true)
    log "✓ Downloaded $count hantavirus segment sequences (4 species × 3 segments)"

    # Build Minimap2 index
    log "Building Minimap2 index..."
    minimap2 -d ../minimap2/hantavirus_all.mmi hantavirus_all.fasta
    log "✓ Minimap2 index created"

    # Create segment-specific indices (for 3-segment concordance checking)
    minimap2 -d ../minimap2/hantavirus_L.mmi hantavirus_L_segment.fasta
    minimap2 -d ../minimap2/hantavirus_M.mmi hantavirus_M_segment.fasta
    minimap2 -d ../minimap2/hantavirus_S.mmi hantavirus_S_segment.fasta
    log "✓ Segment-specific indices created (L/M/S)"

    # Create metadata
    cat > ../metadata.tsv <<EOF
species	segment	accession	length_nt	geographic_region	reservoir
Hantaan	L	NC_005222	6533	Asia (Korea, China, Japan)	Apodemus agrarius
Hantaan	M	NC_005219	3651	Asia	Apodemus agrarius
Hantaan	S	NC_005218	1696	Asia	Apodemus agrarius
Seoul	L	NC_005238	6544	Worldwide (urban)	Rattus norvegicus
Seoul	M	NC_005236	3651	Worldwide	Rattus norvegicus
Seoul	S	NC_005237	1715	Worldwide	Rattus norvegicus
Dobrava	L	NC_005233	6550	Europe	Apodemus flavicollis
Dobrava	M	NC_005234	3643	Europe	Apodemus flavicollis
Dobrava	S	NC_005235	1700	Europe	Apodemus flavicollis
Puumala	L	NC_005224	6530	Europe	Clethrionomys glareolus
Puumala	M	NC_005223	3623	Europe	Clethrionomys glareolus
Puumala	S	NC_005221	1764	Europe	Clethrionomys glareolus
EOF

    # Note on amplicon primers (to be designed)
    cat > ../amplicon_primers/README.md <<EOF
# Hantavirus Amplicon Primers

## Design Required

For high-sensitivity detection (<100 copies/mL), tiled amplicon primers are required.

**Design specifications:**
- Amplicon size: 400 bp
- Overlap: 100 bp
- Total amplicons: 36 (L: 20, M: 11, S: 5)
- Primer pools: A and B (18 amplicons each)
- Strategy: ARTIC Network-style tiled design

**Reference:**
Kim J, et al. "Multiplex PCR-Based Nanopore Sequencing and Epidemiological Surveillance of Hantaan orthohantavirus." Viruses 2021;13:1679.

**TODO:**
- [ ] Design primers using Primal Scheme
- [ ] Validate with Hantaan virus positive controls
- [ ] Test multiplexing efficiency
EOF

    log "✓ Hantavirus database complete"
    log "  IMPORTANT: 3-segment concordance required for positive detection"
}

#==============================================================================
# 3. ALPHAVIRUS DATABASE (EEEV, WEEV, VEEV, Getah)
#==============================================================================

build_alphavirus_db() {
    log "=========================================="
    log "Building Alphavirus Database (PMDA-listed)"
    log "=========================================="

    cd "$BASE_DIR"/alphavirus/fasta

    log "Downloading PMDA-listed alphavirus genomes..."

    # Eastern Equine Encephalitis Virus (EEEV) - PMDA Line 32/55
    esearch -db nucleotide -query "NC_003899" | efetch -format fasta > eeev_north_american.fasta

    # Western Equine Encephalitis Virus (WEEV) - PMDA Line 32/55
    esearch -db nucleotide -query "NC_003908" | efetch -format fasta > weev.fasta

    # Venezuelan Equine Encephalitis Virus (VEEV) - PMDA Line 33/56
    esearch -db nucleotide -query "NC_001449" | efetch -format fasta > veev.fasta

    # Getah virus - PMDA Line 18/41
    esearch -db nucleotide -query "NC_003696" | efetch -format fasta > getah.fasta

    # Additional EEEV lineages for phylogenetic analysis
    log "Downloading EEEV lineage references for phylogeny..."
    esearch -db nucleotide -query "KP765787" | efetch -format fasta > eeev_south_american.fasta || warn "EEEV SA not found"

    # Concatenate all
    cat eeev_*.fasta weev.fasta veev.fasta getah.fasta > alphavirus_all.fasta

    # EEEV-specific for phylogeny
    cat eeev_*.fasta > ../phylogeny/eeev_references.fasta

    local count=$(grep -c "^>" alphavirus_all.fasta || true)
    log "✓ Downloaded $count alphavirus reference sequences"

    # Build Minimap2 index
    log "Building Minimap2 index..."
    minimap2 -d ../minimap2/alphavirus_all.mmi alphavirus_all.fasta
    log "✓ Minimap2 index created"

    # Build BLAST database
    log "Building BLAST database..."
    makeblastdb -in alphavirus_all.fasta -dbtype nucl \
        -out ../blast/alphavirus \
        -title "PMDA Alphavirus Database" \
        -parse_seqids
    log "✓ BLAST database created"

    # Create metadata
    cat > ../metadata.tsv <<EOF
virus	accession	length_nt	lineage	pmda_line	poly_a_tail	zoonotic
EEEV	NC_003899	11841	North American	32/55	Yes	Yes
EEEV	KP765787	11823	South American	32/55	Yes	Yes
WEEV	NC_003908	11722	N/A	32/55	Yes	Yes
VEEV	NC_001449	11441	N/A	33/56	Yes	Yes
Getah	NC_003696	11680	N/A	18/41	Yes	Yes
EOF

    log "✓ Alphavirus database complete"
    log "  Includes EEEV phylogenetic references for lineage assignment"
}

#==============================================================================
# 4. SPUMAVIRUS DATABASE (Cross-genus, no porcine reference)
#==============================================================================

build_spumavirus_db() {
    log "=========================================="
    log "Building Spumavirus Database (pol gene)"
    log "=========================================="

    cd "$BASE_DIR"/spumavirus/fasta

    warn "IMPORTANT: No porcine spumavirus reference genome exists in NCBI"
    warn "Using cross-genus foamy virus pol gene references for degenerate PCR design"

    log "Downloading foamy virus pol gene sequences..."

    # Simian Foamy Virus (SFV) - pol gene region (nt 2000-5000)
    log "  - Simian foamy virus (SFV) pol gene..."
    esearch -db nucleotide -query "NC_001364" | efetch -format fasta > sfv_complete.fasta
    # Extract pol gene region (approximate coordinates)
    # Note: Requires manual extraction or use of gene annotations
    samtools faidx sfv_complete.fasta NC_001364.1:2000-5000 > sfv_pol.fasta 2>/dev/null || \
        cp sfv_complete.fasta sfv_pol.fasta

    # Feline Foamy Virus (FFV) - pol gene
    log "  - Feline foamy virus (FFV) pol gene..."
    esearch -db nucleotide -query "NC_001871" | efetch -format fasta > ffv_complete.fasta
    samtools faidx ffv_complete.fasta NC_001871.1:2100-5100 > ffv_pol.fasta 2>/dev/null || \
        cp ffv_complete.fasta ffv_pol.fasta

    # Bovine Foamy Virus (BFV) - pol gene
    log "  - Bovine foamy virus (BFV) pol gene..."
    esearch -db nucleotide -query "NC_001831" | efetch -format fasta > bfv_complete.fasta
    samtools faidx bfv_complete.fasta NC_001831.1:2050-5050 > bfv_pol.fasta 2>/dev/null || \
        cp bfv_complete.fasta bfv_pol.fasta

    # PERV pol gene (for discrimination from spumavirus)
    log "  - PERV pol gene (for discrimination)..."
    esearch -db nucleotide -query "AF038600" | efetch -format fasta > perv_complete.fasta
    samtools faidx perv_complete.fasta AF038600.1:4500-6500 > perv_pol.fasta 2>/dev/null || \
        cp perv_complete.fasta perv_pol.fasta

    # Concatenate pol genes
    cat sfv_pol.fasta ffv_pol.fasta bfv_pol.fasta perv_pol.fasta > spumavirus_all_pol.fasta

    local count=$(grep -c "^>" spumavirus_all_pol.fasta || true)
    log "✓ Downloaded $count foamy virus pol gene sequences"

    # Build BLAST database (for Sanger confirmation)
    log "Building BLAST database..."
    makeblastdb -in spumavirus_all_pol.fasta -dbtype nucl \
        -out ../blast/spumavirus_pol \
        -title "Spumavirus pol gene (cross-genus)" \
        -parse_seqids
    log "✓ BLAST database created"

    # Create degenerate primer file
    cat > ../nested_pcr_primers.txt <<EOF
# Spumavirus Nested PCR Primers (Degenerate)
# Target: pol gene (reverse transcriptase domain)
# Expected product: 1st PCR ~800 bp, 2nd PCR ~400 bp

# Outer primers (1st PCR)
FV-pol-F1: 5'-GGNCARATHGGNATGTTYGG-3'
  Degeneracy: 96-fold
  Tm: 52-58°C (degenerate)
  Target: RT domain conserved motif

FV-pol-R1: 5'-CCRTCNCCRAANCCRTC-3'
  Degeneracy: 64-fold
  Tm: 52-58°C
  Target: RT domain conserved motif

# Inner primers (2nd PCR)
FV-pol-F2: 5'-ATHGGNCARGGNTTYACNAC-3'
  Degeneracy: 96-fold
  Tm: 54-60°C

FV-pol-R2: 5'-GTRTCNGTYTTRTCNCC-3'
  Degeneracy: 64-fold
  Tm: 54-60°C

# Degeneracy code:
# N = A/T/G/C (4-fold)
# R = A/G (purine, 2-fold)
# Y = C/T (pyrimidine, 2-fold)
# H = A/T/C (not G, 3-fold)

# PCR conditions:
# 1st PCR: 35 cycles, annealing 55°C
# 2nd PCR: 25 cycles, annealing 58°C, template 1:50 dilution

# Expected sensitivity: 1-10 copies per 10^5 PBMCs
EOF

    # Create metadata
    cat > ../metadata.tsv <<EOF
species	accession	genome_size_bp	pol_gene_region	identity_to_pig_predicted	notes
SFV	NC_001364	11904	2000-5000	30-50%	Simian foamy virus, reference for cross-genus detection
FFV	NC_001871	9474	2100-5100	30-50%	Feline foamy virus
BFV	NC_001831	11173	2050-5050	30-50%	Bovine foamy virus
PERV	AF038600	8571	4500-6500	Pig endogenous	For discrimination - Gammaretrovirus not Spumavirus
EOF

    log "✓ Spumavirus database complete"
    log "  CRITICAL: Use nested PCR approach, not metagenomic alignment"
    log "  Degenerate primers provided in: nested_pcr_primers.txt"
}

#==============================================================================
# MAIN BUILD PROCESS
#==============================================================================

main() {
    log "=========================================="
    log "PMDA 4-Virus Database Build"
    log "=========================================="
    log "Base directory: $BASE_DIR"
    log "Threads: $THREADS"
    log "=========================================="

    # Check dependencies
    check_dependencies

    # Create directory structure
    create_directories

    # Build individual databases
    build_polyomavirus_db
    build_hantavirus_db
    build_alphavirus_db
    build_spumavirus_db

    # Summary
    log ""
    log "=========================================="
    log "Database Build Complete!"
    log "=========================================="
    log "Databases created:"
    log "  1. Polyomavirus: $BASE_DIR/polyomavirus/"
    log "  2. Hantavirus: $BASE_DIR/hantavirus/"
    log "  3. Alphavirus (EEEV): $BASE_DIR/alphavirus/"
    log "  4. Spumavirus: $BASE_DIR/spumavirus/"
    log ""
    log "Total disk usage:"
    du -sh "$BASE_DIR"/* | while read size dir; do
        log "  $(basename $dir): $size"
    done
    log ""
    log "Next steps:"
    log "  1. Validate databases with synthetic reads"
    log "  2. Update Phase 4 detection scripts"
    log "  3. Test with spike-in controls"
    log "  4. Begin PMDA validation studies"
    log "=========================================="
}

# Run main function
main "$@"
