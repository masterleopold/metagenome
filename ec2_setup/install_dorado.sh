#!/bin/bash
# MinION Metagenomics Pipeline - Install Dorado
# Standalone script to install Oxford Nanopore Dorado basecaller

set -euo pipefail

# Configuration
INSTALL_DIR="${INSTALL_DIR:-/opt/ont}"
MODEL_DIR="${MODEL_DIR:-/opt/ont/models}"
INSTALL_MODELS="${INSTALL_MODELS:-true}"
CUDA_VERSION="${CUDA_VERSION:-12.2}"

# Colors for output
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

# Check if running as root
check_root() {
    if [[ $EUID -ne 0 ]]; then
        log_error "This script must be run as root"
        exit 1
    fi
}

# Check CUDA availability
check_cuda() {
    log_info "Checking CUDA installation..."

    if command -v nvidia-smi &> /dev/null; then
        log_info "NVIDIA GPU detected:"
        nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv

        if command -v nvcc &> /dev/null; then
            log_info "CUDA compiler version:"
            nvcc --version
        else
            log_warn "CUDA compiler (nvcc) not found. Dorado will use CPU mode."
            return 1
        fi
    else
        log_warn "No NVIDIA GPU detected. Dorado will run in CPU mode."
        return 1
    fi

    return 0
}

# Install system dependencies
install_dependencies() {
    log_info "Installing system dependencies..."

    if command -v dnf &> /dev/null; then
        # Amazon Linux 2023 / RHEL-based
        dnf install -y \
            wget \
            curl \
            tar \
            gzip \
            libgomp \
            libstdc++ \
            zlib
    elif command -v apt-get &> /dev/null; then
        # Ubuntu/Debian-based
        apt-get update
        apt-get install -y \
            wget \
            curl \
            tar \
            gzip \
            libgomp1 \
            libstdc++6 \
            zlib1g
    else
        log_error "Unsupported package manager. Please install dependencies manually."
        exit 1
    fi
}

# Download and install Dorado
install_dorado() {
    log_info "Installing Dorado..."

    # Create installation directory
    mkdir -p "$INSTALL_DIR/dorado"
    cd "$INSTALL_DIR/dorado"

    # Get latest Dorado version
    log_info "Fetching latest Dorado release..."
    DORADO_VERSION=$(curl -s https://api.github.com/repos/nanoporetech/dorado/releases/latest | \
                     grep '"tag_name":' | \
                     sed -E 's/.*"([^"]+)".*/\1/')

    if [ -z "$DORADO_VERSION" ]; then
        log_error "Failed to get Dorado version"
        exit 1
    fi

    log_info "Installing Dorado version: $DORADO_VERSION"

    # Determine architecture and download appropriate binary
    ARCH=$(uname -m)
    OS=$(uname -s | tr '[:upper:]' '[:lower:]')

    if [ "$ARCH" == "x86_64" ] && [ "$OS" == "linux" ]; then
        DORADO_BINARY="dorado-${DORADO_VERSION}-linux-x64.tar.gz"
    elif [ "$ARCH" == "aarch64" ] && [ "$OS" == "linux" ]; then
        DORADO_BINARY="dorado-${DORADO_VERSION}-linux-arm64.tar.gz"
    else
        log_error "Unsupported architecture: $ARCH on $OS"
        exit 1
    fi

    # Download Dorado
    DOWNLOAD_URL="https://github.com/nanoporetech/dorado/releases/download/${DORADO_VERSION}/${DORADO_BINARY}"
    log_info "Downloading from: $DOWNLOAD_URL"

    if ! wget -q --show-progress "$DOWNLOAD_URL"; then
        log_error "Failed to download Dorado"
        exit 1
    fi

    # Extract Dorado
    log_info "Extracting Dorado..."
    tar -xzf "$DORADO_BINARY"
    rm "$DORADO_BINARY"

    # Find the extracted directory
    DORADO_DIR=$(find . -maxdepth 1 -type d -name "dorado-*" | head -1)

    if [ -z "$DORADO_DIR" ]; then
        log_error "Failed to find extracted Dorado directory"
        exit 1
    fi

    # Create symlink for easy access
    ln -sf "$INSTALL_DIR/dorado/$DORADO_DIR/bin/dorado" /usr/local/bin/dorado

    # Verify installation
    if dorado --version &> /dev/null; then
        log_info "Dorado installed successfully:"
        dorado --version
    else
        log_error "Dorado installation verification failed"
        exit 1
    fi
}

# Download Dorado models
download_models() {
    if [ "$INSTALL_MODELS" != "true" ]; then
        log_info "Skipping model download (INSTALL_MODELS=false)"
        return
    fi

    log_info "Downloading Dorado models..."

    # Create models directory
    mkdir -p "$MODEL_DIR"
    cd "$MODEL_DIR"

    # List of models to download
    # Update these based on your sequencing kit and requirements
    MODELS=(
        "dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
        "dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
        "dna_r10.4.1_e8.2_400bps_fast@v4.3.0"
        "dna_r10.4.1_e8.2_260bps_sup@v4.3.0"
        "dna_r10.4.1_e8.2_260bps_hac@v4.3.0"
    )

    # RNA models (uncomment if needed)
    # MODELS+=(
    #     "rna004_130bps_sup@v3.0.1"
    #     "rna004_130bps_hac@v3.0.1"
    # )

    for model in "${MODELS[@]}"; do
        log_info "Downloading model: $model"
        if dorado download --model "$model" --directory "$MODEL_DIR"; then
            log_info "Model $model downloaded successfully"
        else
            log_warn "Failed to download model $model"
        fi
    done

    log_info "Models downloaded to: $MODEL_DIR"
    ls -la "$MODEL_DIR/"
}

# Install pod5 tools
install_pod5_tools() {
    log_info "Installing pod5 tools..."

    if command -v pip3 &> /dev/null; then
        pip3 install --upgrade pip
        pip3 install pod5
        log_info "pod5 tools installed successfully"
    else
        log_warn "pip3 not found. Skipping pod5 installation."
        log_warn "To install manually: pip3 install pod5"
    fi
}

# Create helper scripts
create_helper_scripts() {
    log_info "Creating helper scripts..."

    # Basecalling script
    cat > /usr/local/bin/dorado-basecall << 'SCRIPT'
#!/bin/bash
# Dorado basecalling helper script

set -euo pipefail

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -i INPUT    Input directory with FAST5/POD5 files (required)"
    echo "  -o OUTPUT   Output directory for FASTQ files (required)"
    echo "  -m MODEL    Model name (default: dna_r10.4.1_e8.2_400bps_sup@v4.3.0)"
    echo "  -d DEVICE   CUDA device (default: cuda:all)"
    echo "  -b BATCH    Batch size (default: auto)"
    echo "  -x          Enable duplex calling"
    echo "  -r          Enable RNA mode"
    echo "  -t TRIM     Trim strategy (none/primers/adapters/all)"
    echo "  -q MIN_Q    Minimum qscore filter (default: 9)"
    echo "  -h          Show this help"
    exit 1
}

# Defaults
MODEL="dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
DEVICE="cuda:all"
BATCH="0"  # Auto
DUPLEX=false
RNA=false
TRIM=""
MIN_QSCORE=""

while getopts "i:o:m:d:b:t:q:xrh" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        d) DEVICE="$OPTARG" ;;
        b) BATCH="$OPTARG" ;;
        t) TRIM="$OPTARG" ;;
        q) MIN_QSCORE="$OPTARG" ;;
        x) DUPLEX=true ;;
        r) RNA=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [ -z "${INPUT:-}" ] || [ -z "${OUTPUT:-}" ]; then
    echo "Error: Input and output directories are required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT"

# Build dorado command
if [ "$DUPLEX" = true ]; then
    CMD="dorado duplex"
else
    CMD="dorado basecaller"
fi

# Add model path
if [[ "$MODEL" == *"@"* ]]; then
    # Model with version specified
    MODEL_PATH="/opt/ont/models/$MODEL"
else
    # Try to find model in models directory
    MODEL_PATH="/opt/ont/models/$MODEL"
fi

CMD="$CMD $MODEL_PATH"

# Add input
CMD="$CMD $INPUT"

# Add device
CMD="$CMD --device $DEVICE"

# Add batch size if specified
if [ "$BATCH" != "0" ]; then
    CMD="$CMD --batchsize $BATCH"
fi

# Add trim if specified
if [ -n "$TRIM" ]; then
    CMD="$CMD --trim $TRIM"
fi

# Add min qscore if specified
if [ -n "$MIN_QSCORE" ]; then
    CMD="$CMD --min-qscore $MIN_QSCORE"
fi

# Add RNA flag if needed
if [ "$RNA" = true ]; then
    CMD="$CMD --rna"
fi

# Run basecalling
echo "Running: $CMD"
echo "Output will be saved to: $OUTPUT"

if [ "$DUPLEX" = true ]; then
    OUTPUT_FILE="$OUTPUT/basecalls_duplex.fastq"
else
    OUTPUT_FILE="$OUTPUT/basecalls.fastq"
fi

# Execute basecalling with progress
$CMD --emit-fastq 2>&1 | tee "$OUTPUT/basecalling.log" > "$OUTPUT_FILE"

# Compress output
echo "Compressing output..."
pigz -p 4 "$OUTPUT_FILE"

# Generate summary
echo "Generating summary..."
dorado summary "${OUTPUT_FILE}.gz" > "$OUTPUT/sequencing_summary.txt"

echo "Basecalling completed!"
echo "Output: ${OUTPUT_FILE}.gz"
echo "Summary: $OUTPUT/sequencing_summary.txt"
echo "Log: $OUTPUT/basecalling.log"
SCRIPT

    chmod +x /usr/local/bin/dorado-basecall

    # Model listing script
    cat > /usr/local/bin/dorado-list-models << 'SCRIPT'
#!/bin/bash
# List available Dorado models

echo "Available Dorado models in /opt/ont/models:"
echo "=========================================="
ls -la /opt/ont/models/ 2>/dev/null | grep -E '^d' | awk '{print $NF}'

echo ""
echo "To download additional models:"
echo "  dorado download --model MODEL_NAME"
echo ""
echo "Available models online:"
echo "  dorado download --list"
SCRIPT

    chmod +x /usr/local/bin/dorado-list-models

    log_info "Helper scripts created:"
    log_info "  - /usr/local/bin/dorado-basecall"
    log_info "  - /usr/local/bin/dorado-list-models"
}

# Create systemd service for GPU monitoring (optional)
create_monitoring_service() {
    if check_cuda; then
        cat > /etc/systemd/system/gpu-monitor.service << 'SERVICE'
[Unit]
Description=GPU Monitoring for Dorado
After=multi-user.target

[Service]
Type=simple
ExecStart=/usr/bin/nvidia-smi dmon -s pucvmet -d 60
Restart=always
User=root

[Install]
WantedBy=multi-user.target
SERVICE

        systemctl daemon-reload
        systemctl enable gpu-monitor.service
        log_info "GPU monitoring service created"
    fi
}

# Main installation
main() {
    log_info "Starting Dorado installation..."

    check_root
    install_dependencies
    check_cuda
    install_dorado
    download_models
    install_pod5_tools
    create_helper_scripts
    create_monitoring_service

    log_info "======================================"
    log_info "Dorado installation completed!"
    log_info "======================================"
    log_info "Dorado binary: /usr/local/bin/dorado"
    log_info "Models directory: $MODEL_DIR"
    log_info "Helper scripts:"
    log_info "  - dorado-basecall: Simplified basecalling"
    log_info "  - dorado-list-models: List available models"
    log_info ""
    log_info "Quick test:"
    log_info "  dorado --version"
    log_info ""
    log_info "Example usage:"
    log_info "  dorado-basecall -i /data/fast5 -o /data/fastq"
    log_info "======================================"
}

# Run main function
main "$@"