#!/bin/bash
# MinION Metagenomics Pipeline - Install BLAST+
# Installs BLAST+ and Diamond for sequence similarity searches

set -euo pipefail

INSTALL_DIR="${INSTALL_DIR:-/opt/tools}"
DB_DIR="${DB_DIR:-/opt/databases/blast}"
BLAST_VERSION="${BLAST_VERSION:-2.15.0}"
DIAMOND_VERSION="${DIAMOND_VERSION:-2.1.8}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

[[ $EUID -ne 0 ]] && log_error "Run as root"

mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR"
mkdir -p "$DB_DIR"

# Install BLAST+
log_info "Installing BLAST+ v${BLAST_VERSION}..."
wget -q "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
tar -xzf "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
cp "ncbi-blast-${BLAST_VERSION}+/bin/"* /usr/local/bin/
rm -rf "ncbi-blast-${BLAST_VERSION}+"*

# Install Diamond (faster protein alignment)
log_info "Installing Diamond v${DIAMOND_VERSION}..."
wget -q "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz"
tar -xzf diamond-linux64.tar.gz
mv diamond /usr/local/bin/
rm diamond-linux64.tar.gz

# Install BLAST database update script
cat > /usr/local/bin/update-blast-db << 'SCRIPT'
#!/bin/bash
set -euo pipefail

DB_DIR="${1:-/opt/databases/blast}"
DB_TYPE="${2:-nt}"  # nt, nr, refseq_rna, refseq_protein, swissprot

mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "Downloading BLAST database: $DB_TYPE..."
update_blastdb.pl --decompress --passive "$DB_TYPE"

echo "Database $DB_TYPE downloaded to $DB_DIR"
SCRIPT

chmod +x /usr/local/bin/update-blast-db

# Create BLAST search script
cat > /usr/local/bin/run-blast-search << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT -d DATABASE -o OUTPUT [-p PROGRAM] [-e EVALUE] [-t THREADS]"
    echo "  Programs: blastn, blastp, blastx, tblastn, tblastx"
    exit 1
}

PROGRAM="blastn"
EVALUE="1e-10"
THREADS=8
MAX_TARGETS=10

while getopts "i:d:o:p:e:t:m:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        d) DATABASE="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        p) PROGRAM="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MAX_TARGETS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${DATABASE:-}" || -z "${OUTPUT:-}" ]] && usage

mkdir -p "$(dirname "$OUTPUT")"

echo "Running BLAST search..."
echo "Program: $PROGRAM"
echo "Database: $DATABASE"
echo "E-value: $EVALUE"

$PROGRAM \
    -query "$INPUT" \
    -db "$DATABASE" \
    -out "$OUTPUT" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -evalue "$EVALUE" \
    -num_threads "$THREADS" \
    -max_target_seqs "$MAX_TARGETS"

# Generate summary
echo "Generating summary..."
awk '{print $2}' "$OUTPUT" | cut -d'|' -f4 | sort | uniq -c | sort -rn > "${OUTPUT%.txt}_species_summary.txt"

echo "BLAST search completed: $OUTPUT"
SCRIPT

chmod +x /usr/local/bin/run-blast-search

# Create Diamond search script
cat > /usr/local/bin/run-diamond-search << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT -d DATABASE -o OUTPUT [-m MODE] [-e EVALUE] [-t THREADS]"
    echo "  Modes: blastx (default), blastp"
    exit 1
}

MODE="blastx"
EVALUE="1e-10"
THREADS=8

while getopts "i:d:o:m:e:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        d) DATABASE="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        m) MODE="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${DATABASE:-}" || -z "${OUTPUT:-}" ]] && usage

mkdir -p "$(dirname "$OUTPUT")"

echo "Running Diamond $MODE search..."

diamond $MODE \
    -q "$INPUT" \
    -d "$DATABASE" \
    -o "$OUTPUT" \
    --evalue "$EVALUE" \
    --threads "$THREADS" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --max-target-seqs 10 \
    --sensitive

echo "Diamond search completed: $OUTPUT"
SCRIPT

chmod +x /usr/local/bin/run-diamond-search

# Create Diamond database build script
cat > /usr/local/bin/build-diamond-db << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_FASTA -o OUTPUT_DB"
    exit 1
}

while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" ]] && usage

echo "Building Diamond database..."
diamond makedb --in "$INPUT" --db "$OUTPUT" --threads 8

echo "Database created: ${OUTPUT}.dmnd"
SCRIPT

chmod +x /usr/local/bin/build-diamond-db

log_info "======================================"
log_info "BLAST+ and Diamond installed!"
log_info "Tools:"
log_info "  - BLAST+ v${BLAST_VERSION}"
log_info "  - Diamond v${DIAMOND_VERSION}"
log_info "Scripts:"
log_info "  - update-blast-db: Download NCBI databases"
log_info "  - run-blast-search: Run BLAST analysis"
log_info "  - run-diamond-search: Run Diamond search"
log_info "  - build-diamond-db: Create Diamond database"
log_info "======================================"