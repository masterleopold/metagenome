#!/bin/bash
# MinION Metagenomics Pipeline - Install QC Tools
# Installs quality control tools for Nanopore sequencing data

set -euo pipefail

# Configuration
INSTALL_DIR="${INSTALL_DIR:-/opt/tools}"
FASTQC_VERSION="${FASTQC_VERSION:-0.12.1}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_root() {
    if [[ $EUID -ne 0 ]]; then
        log_error "This script must be run as root"
        exit 1
    fi
}

install_dependencies() {
    log_info "Installing system dependencies..."

    if command -v dnf &> /dev/null; then
        dnf install -y \
            java-11-amazon-corretto \
            python3 python3-pip \
            wget unzip \
            perl perl-GD
    elif command -v apt-get &> /dev/null; then
        apt-get update
        apt-get install -y \
            default-jre \
            python3 python3-pip \
            wget unzip \
            perl libgd-perl
    fi
}

install_fastqc() {
    log_info "Installing FastQC v${FASTQC_VERSION}..."

    mkdir -p "$INSTALL_DIR"
    cd "$INSTALL_DIR"

    # Download FastQC
    wget -q "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip"
    unzip -q "fastqc_v${FASTQC_VERSION}.zip"
    rm "fastqc_v${FASTQC_VERSION}.zip"

    # Make executable
    chmod +x FastQC/fastqc

    # Create symlink
    ln -sf "$INSTALL_DIR/FastQC/fastqc" /usr/local/bin/fastqc

    log_info "FastQC installed successfully"
}

install_nanoplot() {
    log_info "Installing NanoPlot..."

    pip3 install --upgrade pip
    pip3 install NanoPlot

    # Verify installation
    if command -v NanoPlot &> /dev/null; then
        log_info "NanoPlot installed successfully"
        NanoPlot --version
    else
        log_error "NanoPlot installation failed"
        exit 1
    fi
}

install_nanostat() {
    log_info "Installing NanoStat..."

    pip3 install nanostat

    if command -v NanoStat &> /dev/null; then
        log_info "NanoStat installed successfully"
    else
        log_error "NanoStat installation failed"
        exit 1
    fi
}

install_nanofilt() {
    log_info "Installing NanoFilt..."

    pip3 install nanofilt

    if command -v NanoFilt &> /dev/null; then
        log_info "NanoFilt installed successfully"
    else
        log_error "NanoFilt installation failed"
        exit 1
    fi
}

install_pycoqc() {
    log_info "Installing PycoQC..."

    pip3 install pycoQC

    if command -v pycoQC &> /dev/null; then
        log_info "PycoQC installed successfully"
        pycoQC --version
    else
        log_error "PycoQC installation failed"
        exit 1
    fi
}

install_multiqc() {
    log_info "Installing MultiQC..."

    pip3 install multiqc

    if command -v multiqc &> /dev/null; then
        log_info "MultiQC installed successfully"
        multiqc --version
    else
        log_error "MultiQC installation failed"
        exit 1
    fi
}

install_longqc() {
    log_info "Installing LongQC (for long-read specific QC)..."

    # Install dependencies
    pip3 install matplotlib scipy pandas scikit-learn

    # Clone and install LongQC
    cd "$INSTALL_DIR"
    git clone https://github.com/yfukasawa/LongQC.git
    cd LongQC

    # Create symlink
    ln -sf "$INSTALL_DIR/LongQC/longQC.py" /usr/local/bin/longqc

    log_info "LongQC installed"
}

create_qc_scripts() {
    log_info "Creating QC helper scripts..."

    # Comprehensive QC script
    cat > /usr/local/bin/run-nanopore-qc << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_DIR [-s SUMMARY_FILE] [-t THREADS]"
    echo ""
    echo "Options:"
    echo "  -i INPUT_FASTQ   Input FASTQ file (required)"
    echo "  -o OUTPUT_DIR    Output directory (required)"
    echo "  -s SUMMARY_FILE  Sequencing summary file (optional, for PycoQC)"
    echo "  -t THREADS       Number of threads (default: 4)"
    echo "  -h               Show this help"
    exit 1
}

THREADS=4

while getopts "i:o:s:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        s) SUMMARY="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "${INPUT:-}" ] || [ -z "${OUTPUT:-}" ]; then
    echo "Error: Input file and output directory are required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT"

echo "========================================="
echo "Running Nanopore QC Analysis"
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo "========================================="

# FastQC
echo "Running FastQC..."
mkdir -p "$OUTPUT/fastqc"
fastqc -o "$OUTPUT/fastqc" -t "$THREADS" "$INPUT"

# NanoPlot
echo "Running NanoPlot..."
NanoPlot \
    --fastq "$INPUT" \
    --outdir "$OUTPUT/nanoplot" \
    --prefix sample_ \
    --plots hex dot \
    --N50 \
    --threads "$THREADS" \
    --title "Nanopore QC Report" \
    --color darkblue

# NanoStat
echo "Running NanoStat..."
NanoStat --fastq "$INPUT" --threads "$THREADS" > "$OUTPUT/nanostats.txt"

# PycoQC (if summary file provided)
if [ -n "${SUMMARY:-}" ] && [ -f "$SUMMARY" ]; then
    echo "Running PycoQC..."
    pycoQC \
        --summary_file "$SUMMARY" \
        --output_file "$OUTPUT/pycoqc_report.html" \
        --min_pass_qual 7 \
        --min_pass_len 200
else
    echo "Skipping PycoQC (no summary file provided)"
fi

# Generate combined report
echo "Generating summary report..."
cat > "$OUTPUT/qc_summary.txt" << EOF
QC Analysis Summary
===================
Date: $(date)
Input: $INPUT

NanoStat Summary:
-----------------
$(cat "$OUTPUT/nanostats.txt")

Additional reports:
- FastQC: $OUTPUT/fastqc/
- NanoPlot: $OUTPUT/nanoplot/
EOF

if [ -n "${SUMMARY:-}" ]; then
    echo "- PycoQC: $OUTPUT/pycoqc_report.html" >> "$OUTPUT/qc_summary.txt"
fi

echo ""
echo "QC analysis completed!"
echo "Summary: $OUTPUT/qc_summary.txt"
echo "========================================="
SCRIPT

    chmod +x /usr/local/bin/run-nanopore-qc

    # Quality filtering script
    cat > /usr/local/bin/filter-nanopore-reads << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_FASTQ [options]"
    echo ""
    echo "Options:"
    echo "  -i INPUT   Input FASTQ file (required)"
    echo "  -o OUTPUT  Output FASTQ file (required)"
    echo "  -q QSCORE  Minimum quality score (default: 9)"
    echo "  -l LENGTH  Minimum read length (default: 200)"
    echo "  -L LENGTH  Maximum read length (default: none)"
    echo "  -h HEADCROP Remove n nucleotides from start"
    echo "  -t TAILCROP Remove n nucleotides from end"
    echo "  -h         Show this help"
    exit 1
}

MIN_Q=9
MIN_LEN=200
MAX_LEN=""
HEADCROP=""
TAILCROP=""

while getopts "i:o:q:l:L:h:t:" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        q) MIN_Q="$OPTARG" ;;
        l) MIN_LEN="$OPTARG" ;;
        L) MAX_LEN="$OPTARG" ;;
        h) HEADCROP="$OPTARG" ;;
        t) TAILCROP="$OPTARG" ;;
        *) usage ;;
    esac
done

if [ -z "${INPUT:-}" ] || [ -z "${OUTPUT:-}" ]; then
    usage
fi

echo "Filtering reads..."
echo "  Min quality: $MIN_Q"
echo "  Min length: $MIN_LEN"

# Build NanoFilt command
CMD="NanoFilt"
CMD="$CMD -q $MIN_Q"
CMD="$CMD -l $MIN_LEN"

if [ -n "$MAX_LEN" ]; then
    CMD="$CMD --maxlength $MAX_LEN"
fi

if [ -n "$HEADCROP" ]; then
    CMD="$CMD --headcrop $HEADCROP"
fi

if [ -n "$TAILCROP" ]; then
    CMD="$CMD --tailcrop $TAILCROP"
fi

# Run filtering
if [[ "$INPUT" == *.gz ]]; then
    gunzip -c "$INPUT" | $CMD | gzip > "$OUTPUT"
else
    $CMD "$INPUT" > "$OUTPUT"
fi

echo "Filtering completed!"
echo "Output: $OUTPUT"
SCRIPT

    chmod +x /usr/local/bin/filter-nanopore-reads

    log_info "QC scripts created:"
    log_info "  - /usr/local/bin/run-nanopore-qc"
    log_info "  - /usr/local/bin/filter-nanopore-reads"
}

main() {
    log_info "Installing Nanopore QC tools..."

    check_root
    install_dependencies
    install_fastqc
    install_nanoplot
    install_nanostat
    install_nanofilt
    install_pycoqc
    install_multiqc
    install_longqc
    create_qc_scripts

    log_info "======================================"
    log_info "QC tools installation completed!"
    log_info "======================================"
    log_info "Installed tools:"
    log_info "  - FastQC: General sequence QC"
    log_info "  - NanoPlot: Nanopore-specific plotting"
    log_info "  - NanoStat: Nanopore statistics"
    log_info "  - NanoFilt: Read filtering"
    log_info "  - PycoQC: Interactive QC reports"
    log_info "  - MultiQC: Aggregate QC reports"
    log_info "  - LongQC: Long-read specific QC"
    log_info ""
    log_info "Helper scripts:"
    log_info "  - run-nanopore-qc: Complete QC pipeline"
    log_info "  - filter-nanopore-reads: Quality filtering"
    log_info "======================================"
}

main "$@"