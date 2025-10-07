#!/bin/bash
# MinION Metagenomics Pipeline - Install Alignment Tools
# Installs alignment tools for long-read sequencing data

set -euo pipefail

INSTALL_DIR="${INSTALL_DIR:-/opt/tools}"
THREADS="${THREADS:-8}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

[[ $EUID -ne 0 ]] && log_error "Run as root"

# Install dependencies
log_info "Installing dependencies..."
if command -v dnf &> /dev/null; then
    dnf install -y gcc gcc-c++ make cmake zlib-devel bzip2-devel xz-devel libcurl-devel ncurses-devel
else
    apt-get update && apt-get install -y build-essential cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses5-dev
fi

mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR"

# Install Minimap2
log_info "Installing Minimap2..."
git clone https://github.com/lh3/minimap2
cd minimap2 && make -j"$THREADS"
cp minimap2 /usr/local/bin/
cd ..

# Install Samtools
log_info "Installing Samtools..."
wget -q https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=/usr/local && make -j"$THREADS" && make install
cd .. && rm -rf samtools-1.18*

# Install BCFtools
log_info "Installing BCFtools..."
wget -q https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
tar -xjf bcftools-1.18.tar.bz2
cd bcftools-1.18
./configure --prefix=/usr/local && make -j"$THREADS" && make install
cd .. && rm -rf bcftools-1.18*

# Install BWA (for short reads if needed)
log_info "Installing BWA..."
git clone https://github.com/lh3/bwa.git
cd bwa && make -j"$THREADS"
cp bwa /usr/local/bin/
cd ..

# Install Bowtie2 (optional, for short reads)
log_info "Installing Bowtie2..."
wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.5.2/bowtie2-2.5.2-linux-x86_64.zip
unzip -q bowtie2-2.5.2-linux-x86_64.zip
cp bowtie2-2.5.2-linux-x86_64/bowtie2* /usr/local/bin/
rm -rf bowtie2-2.5.2*

# Create alignment script
cat > /usr/local/bin/align-to-reference << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT -r REFERENCE -o OUTPUT [-t THREADS] [-p PRESET]"
    echo "  -p PRESET: map-ont (default), splice, asm5, asm10"
    exit 1
}

THREADS=8
PRESET="map-ont"

while getopts "i:r:o:t:p:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        r) REF="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PRESET="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${REF:-}" || -z "${OUTPUT:-}" ]] && usage

mkdir -p "$(dirname "$OUTPUT")"

echo "Aligning with minimap2..."
minimap2 -ax "$PRESET" -t "$THREADS" "$REF" "$INPUT" | \
    samtools sort -@ "$THREADS" -o "$OUTPUT"

samtools index -@ "$THREADS" "$OUTPUT"

echo "Generating statistics..."
samtools flagstat "$OUTPUT" > "${OUTPUT%.bam}.flagstat"
samtools idxstats "$OUTPUT" > "${OUTPUT%.bam}.idxstats"

echo "Alignment completed: $OUTPUT"
SCRIPT

chmod +x /usr/local/bin/align-to-reference

log_info "Alignment tools installed successfully!"