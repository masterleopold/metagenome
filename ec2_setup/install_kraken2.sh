#!/bin/bash
# MinION Metagenomics Pipeline - Install Kraken2
# Installs Kraken2, Bracken, and KrakenTools for taxonomic classification

set -euo pipefail

INSTALL_DIR="${INSTALL_DIR:-/opt/tools}"
DB_DIR="${DB_DIR:-/opt/databases/kraken2}"
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
    dnf install -y gcc gcc-c++ make perl wget rsync
else
    apt-get update && apt-get install -y build-essential perl wget rsync
fi

mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR"

# Install Kraken2
log_info "Installing Kraken2..."
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh "$INSTALL_DIR/kraken2_install"
ln -sf "$INSTALL_DIR/kraken2_install/kraken2" /usr/local/bin/
ln -sf "$INSTALL_DIR/kraken2_install/kraken2-build" /usr/local/bin/
ln -sf "$INSTALL_DIR/kraken2_install/kraken2-inspect" /usr/local/bin/
cd ..

# Install Bracken
log_info "Installing Bracken..."
git clone https://github.com/jenniferlu717/Bracken.git
cd Bracken
./install_bracken.sh
ln -sf "$INSTALL_DIR/Bracken/bracken" /usr/local/bin/
ln -sf "$INSTALL_DIR/Bracken/bracken-build" /usr/local/bin/
cd ..

# Install KrakenTools
log_info "Installing KrakenTools..."
git clone https://github.com/jenniferlu717/KrakenTools.git
ln -sf "$INSTALL_DIR/KrakenTools/kreport2mpa.py" /usr/local/bin/
ln -sf "$INSTALL_DIR/KrakenTools/kreport2krona.py" /usr/local/bin/
ln -sf "$INSTALL_DIR/KrakenTools/combine_kreports.py" /usr/local/bin/

# Install Krona for visualization
log_info "Installing Krona..."
wget -q https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar
tar -xf KronaTools-2.8.1.tar
cd KronaTools-2.8.1
./install.pl --prefix=/usr/local
./updateTaxonomy.sh
cd .. && rm -rf KronaTools-2.8.1*

# Create database download script
cat > /usr/local/bin/download-kraken2-db << 'SCRIPT'
#!/bin/bash
set -euo pipefail

DB_DIR="${1:-/opt/databases/kraken2}"
DB_TYPE="${2:-standard}"  # standard, viral, minusb, standard-8

echo "Downloading Kraken2 $DB_TYPE database to $DB_DIR..."
mkdir -p "$DB_DIR"

if [[ "$DB_TYPE" == "standard" ]]; then
    kraken2-build --standard --db "$DB_DIR" --threads 8
elif [[ "$DB_TYPE" == "viral" ]]; then
    kraken2-build --download-taxonomy --db "$DB_DIR"
    kraken2-build --download-library viral --db "$DB_DIR"
    kraken2-build --build --db "$DB_DIR" --threads 8
elif [[ "$DB_TYPE" == "minusb" ]]; then
    # PlusPFP minus bacteria
    kraken2-build --download-taxonomy --db "$DB_DIR"
    kraken2-build --download-library archaea --db "$DB_DIR"
    kraken2-build --download-library viral --db "$DB_DIR"
    kraken2-build --download-library plasmid --db "$DB_DIR"
    kraken2-build --download-library human --db "$DB_DIR"
    kraken2-build --download-library fungi --db "$DB_DIR"
    kraken2-build --download-library plant --db "$DB_DIR"
    kraken2-build --download-library protozoa --db "$DB_DIR"
    kraken2-build --build --db "$DB_DIR" --threads 8
else
    echo "Unknown database type: $DB_TYPE"
    exit 1
fi

# Build Bracken database
echo "Building Bracken database..."
bracken-build -d "$DB_DIR" -t 8 -k 35 -l 150

echo "Database ready at: $DB_DIR"
SCRIPT

chmod +x /usr/local/bin/download-kraken2-db

# Create classification script
cat > /usr/local/bin/run-kraken2-classification << 'SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT -d DATABASE -o OUTPUT_DIR [-t THREADS] [-c CONFIDENCE]"
    exit 1
}

THREADS=8
CONFIDENCE=0.05

while getopts "i:d:o:t:c:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        d) DATABASE="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        c) CONFIDENCE="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${DATABASE:-}" || -z "${OUTPUT:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Running Kraken2 classification..."
kraken2 --db "$DATABASE" \
    --threads "$THREADS" \
    --confidence "$CONFIDENCE" \
    --report "$OUTPUT/kraken2_report.txt" \
    --output "$OUTPUT/kraken2_output.txt" \
    --report-zero-counts \
    --use-names \
    "$INPUT"

# Bracken abundance estimation
echo "Running Bracken..."
for LEVEL in S G F; do
    bracken -d "$DATABASE" \
        -i "$OUTPUT/kraken2_report.txt" \
        -o "$OUTPUT/bracken_${LEVEL}.txt" \
        -r 150 \
        -l "$LEVEL" \
        -t 10
done

# Generate Krona plot
echo "Generating Krona visualization..."
kreport2krona.py -r "$OUTPUT/kraken2_report.txt" -o "$OUTPUT/krona.txt"
ktImportText "$OUTPUT/krona.txt" -o "$OUTPUT/krona.html"

echo "Classification completed!"
echo "Results in: $OUTPUT"
SCRIPT

chmod +x /usr/local/bin/run-kraken2-classification

log_info "======================================"
log_info "Kraken2 installation completed!"
log_info "Tools: kraken2, bracken, KrakenTools, Krona"
log_info "Scripts:"
log_info "  - download-kraken2-db: Download databases"
log_info "  - run-kraken2-classification: Run analysis"
log_info "======================================"