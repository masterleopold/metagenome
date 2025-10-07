#!/bin/bash
# MinION Metagenomics Pipeline - Install Assembly Tools
# Installs assemblers for long-read and metagenomic data

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
    dnf install -y gcc gcc-c++ make cmake boost-devel zlib-devel python3 python3-pip
else
    apt-get update && apt-get install -y build-essential cmake libboost-all-dev zlib1g-dev python3 python3-pip
fi

mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR"

# Install Flye (long-read assembler)
log_info "Installing Flye..."
pip3 install flye

# Install Canu (long-read assembler)
log_info "Installing Canu..."
wget -q https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz
tar -xJf canu-2.2.Linux-amd64.tar.xz
ln -sf "$INSTALL_DIR/canu-2.2/bin/canu" /usr/local/bin/
rm canu-2.2.Linux-amd64.tar.xz

# Install MEGAHIT (metagenome assembler)
log_info "Installing MEGAHIT..."
wget -q https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar -xzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
ln -sf "$INSTALL_DIR/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit" /usr/local/bin/
rm MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

# Install metaSPAdes (hybrid assembler)
log_info "Installing SPAdes/metaSPAdes..."
wget -q https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar -xzf SPAdes-3.15.5-Linux.tar.gz
ln -sf "$INSTALL_DIR/SPAdes-3.15.5-Linux/bin/spades.py" /usr/local/bin/
ln -sf "$INSTALL_DIR/SPAdes-3.15.5-Linux/bin/metaspades.py" /usr/local/bin/
rm SPAdes-3.15.5-Linux.tar.gz

# Install Raven (fast long-read assembler)
log_info "Installing Raven..."
git clone --recursive https://github.com/lbcb-sci/raven.git
cd raven && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make -j"$THREADS"
cp bin/raven /usr/local/bin/
cd ../.. && rm -rf raven

# Create assembly pipeline script
cat > /usr/local/bin/run-assembly << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT_DIR -a ASSEMBLER [-t THREADS] [-g GENOME_SIZE]"
    echo "  Assemblers: flye, canu, megahit, metaspades, raven"
    echo "  -g GENOME_SIZE: Expected genome size (e.g., 5m, 100m)"
    exit 1
}

ASSEMBLER="flye"
THREADS=8
GENOME_SIZE="5m"

while getopts "i:o:a:t:g:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        a) ASSEMBLER="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        g) GENOME_SIZE="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Running $ASSEMBLER assembly..."
echo "Input: $INPUT"
echo "Output: $OUTPUT"

case "$ASSEMBLER" in
    flye)
        flye --nano-raw "$INPUT" \
            --out-dir "$OUTPUT" \
            --threads "$THREADS" \
            --meta \
            --genome-size "$GENOME_SIZE"
        ;;
    canu)
        canu -p assembly -d "$OUTPUT" \
            genomeSize="$GENOME_SIZE" \
            -nanopore "$INPUT" \
            maxThreads="$THREADS"
        ;;
    megahit)
        # Convert FASTQ to FASTA if needed
        if [[ "$INPUT" == *.fastq* ]]; then
            seqtk seq -A "$INPUT" > "$OUTPUT/input.fasta"
            INPUT_FASTA="$OUTPUT/input.fasta"
        else
            INPUT_FASTA="$INPUT"
        fi
        megahit -r "$INPUT_FASTA" \
            -o "$OUTPUT" \
            -t "$THREADS" \
            --presets meta-large
        ;;
    metaspades)
        metaspades.py --nanopore "$INPUT" \
            -o "$OUTPUT" \
            -t "$THREADS" \
            -m 32
        ;;
    raven)
        raven -t "$THREADS" "$INPUT" > "$OUTPUT/assembly.fasta"
        ;;
    *)
        echo "Unknown assembler: $ASSEMBLER"
        exit 1
        ;;
esac

# Generate assembly statistics
echo "Generating assembly statistics..."
if [ -f "$OUTPUT/assembly.fasta" ]; then
    ASSEMBLY="$OUTPUT/assembly.fasta"
elif [ -f "$OUTPUT/contigs.fasta" ]; then
    ASSEMBLY="$OUTPUT/contigs.fasta"
elif [ -f "$OUTPUT/final.contigs.fa" ]; then
    ASSEMBLY="$OUTPUT/final.contigs.fa"
else
    echo "Warning: Could not find assembly file"
    exit 0
fi

# Calculate basic stats
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' "$ASSEMBLY" | \
    awk 'BEGIN{n=0;sum=0;} /^>/{n++;next} {sum+=$1; len[n]=$1} END{
        asort(len, sorted);
        total=0;
        for(i=1;i<=n;i++) total+=sorted[i];
        n50_sum=0;
        for(i=n;i>=1;i--){
            n50_sum+=sorted[i];
            if(n50_sum>=total/2){
                n50=sorted[i];
                break;
            }
        }
        print "Number of contigs:", n;
        print "Total length:", total;
        print "Longest contig:", sorted[n];
        print "N50:", n50;
        print "Mean length:", total/n;
    }' > "$OUTPUT/assembly_stats.txt"

cat "$OUTPUT/assembly_stats.txt"
echo "Assembly completed!"
SCRIPT

chmod +x /usr/local/bin/run-assembly

log_info "======================================"
log_info "Assembly tools installed!"
log_info "Tools:"
log_info "  - Flye: Long-read assembler"
log_info "  - Canu: High-quality long-read assembler"
log_info "  - MEGAHIT: Fast metagenome assembler"
log_info "  - SPAdes/metaSPAdes: Hybrid assembler"
log_info "  - Raven: Fast long-read assembler"
log_info "Script:"
log_info "  - run-assembly: Run assembly pipeline"
log_info "======================================"