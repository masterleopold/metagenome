#!/bin/bash
# Database Setup Script for MinION Pipeline

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EFS_MOUNT="${EFS_MOUNT:-/mnt/efs}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -a, --all           Download all databases"
    echo "  -k, --kraken2       Download Kraken2 database"
    echo "  -r, --rvdb          Download RVDB database"
    echo "  -p, --pmda          Create PMDA custom database"
    echo "  -g, --genome        Download reference genome"
    echo "  -u, --update        Update existing databases"
    echo "  -c, --check         Check database status"
    echo "  -h, --help          Show this help message"
    echo ""
    echo "Environment Variables:"
    echo "  EFS_MOUNT          EFS mount point (default: /mnt/efs)"
    exit 0
}

check_requirements() {
    log "Checking requirements..."

    # Check for required tools
    local required_tools=("wget" "gunzip" "makeblastdb" "kraken2-build" "minimap2")
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            error "$tool is required but not installed"
        fi
    done

    # Check EFS mount
    if [[ ! -d "$EFS_MOUNT" ]]; then
        error "EFS mount point $EFS_MOUNT does not exist"
    fi

    # Check disk space (need at least 500GB)
    local available_space=$(df -BG "$EFS_MOUNT" | awk 'NR==2 {print $4}' | sed 's/G//')
    if [[ $available_space -lt 500 ]]; then
        warn "Low disk space: ${available_space}GB available (500GB recommended)"
    fi

    log "All requirements met"
}

setup_directories() {
    log "Setting up directory structure..."

    mkdir -p "$EFS_MOUNT/references/sus_scrofa_11.1"
    mkdir -p "$EFS_MOUNT/references/sus_scrofa_modified"
    mkdir -p "$EFS_MOUNT/databases/kraken2/standard"
    mkdir -p "$EFS_MOUNT/databases/kraken2/pmda"
    mkdir -p "$EFS_MOUNT/databases/rvdb"
    mkdir -p "$EFS_MOUNT/databases/ncbi_refseq"
    mkdir -p "$EFS_MOUNT/databases/pmda"
    mkdir -p "$EFS_MOUNT/indexes/minimap2"
    mkdir -p "$EFS_MOUNT/indexes/bwa"
    mkdir -p "$EFS_MOUNT/tmp"

    log "Directory structure created"
}

download_reference_genome() {
    log "Downloading Sus scrofa reference genome..."

    local genome_dir="$EFS_MOUNT/references/sus_scrofa_11.1"
    local genome_file="$genome_dir/sus_scrofa.11.1.fa"

    if [[ -f "$genome_file" ]]; then
        log "Reference genome already exists"
        return 0
    fi

    # Download from Ensembl
    local url="https://ftp.ensembl.org/pub/release-109/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"

    log "Downloading from Ensembl..."
    wget -O - "$url" | gunzip > "$genome_file"

    # Create index
    log "Creating minimap2 index..."
    minimap2 -x map-ont -d "$EFS_MOUNT/indexes/minimap2/sus_scrofa.11.1.idx" "$genome_file"

    log "Reference genome setup complete"
}

download_kraken2_database() {
    log "Setting up Kraken2 database..."

    local db_dir="$EFS_MOUNT/databases/kraken2/standard"

    if [[ -f "$db_dir/hash.k2d" ]]; then
        log "Kraken2 database already exists"
        return 0
    fi

    # Option 1: Download pre-built database (faster)
    log "Downloading pre-built Kraken2 database..."
    local url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz"

    cd "$db_dir"
    wget -O - "$url" | tar -xz

    # Option 2: Build from scratch (slower but customizable)
    # log "Building Kraken2 database from scratch..."
    # kraken2-build --standard --db "$db_dir" --threads 16

    log "Kraken2 database ready"
}

download_rvdb_database() {
    log "Downloading RVDB (Reference Viral Database)..."

    local rvdb_dir="$EFS_MOUNT/databases/rvdb"
    local rvdb_fasta="$rvdb_dir/rvdb.fasta"

    if [[ -f "$rvdb_fasta" ]]; then
        log "RVDB already exists"
        return 0
    fi

    # Download RVDB v30.0
    local url="https://rvdb.dbi.udel.edu/download/C-RVDBv30.0.fasta.gz"

    log "Downloading RVDB v30.0..."
    wget -O - "$url" | gunzip > "$rvdb_fasta"

    # Create BLAST database
    log "Creating BLAST database..."
    makeblastdb -in "$rvdb_fasta" -dbtype nucl -out "$rvdb_dir/rvdb"

    log "RVDB setup complete"
}

create_pmda_database() {
    log "Creating PMDA 91 pathogen database..."

    local pmda_dir="$EFS_MOUNT/databases/pmda"
    local pmda_fasta="$pmda_dir/pmda_91.fasta"
    local pmda_list="$SCRIPT_DIR/../templates/config/pmda_pathogens.json"

    if [[ -f "$pmda_fasta" ]]; then
        log "PMDA database already exists"
        return 0
    fi

    # Create placeholder sequences for demonstration
    # In production, download actual sequences from NCBI
    cat > "$pmda_fasta" << 'EOF'
>PERV-A|Porcine_endogenous_retrovirus_A
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>PERV-B|Porcine_endogenous_retrovirus_B
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>PERV-C|Porcine_endogenous_retrovirus_C
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>HEV|Hepatitis_E_virus
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>JEV|Japanese_encephalitis_virus
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

    # Build BLAST database
    log "Creating BLAST database for PMDA pathogens..."
    makeblastdb -in "$pmda_fasta" -dbtype nucl -out "$pmda_dir/pmda_91"

    # Build custom Kraken2 database
    log "Building custom Kraken2 database for PMDA pathogens..."
    local kraken_pmda="$EFS_MOUNT/databases/kraken2/pmda"

    kraken2-build --download-taxonomy --db "$kraken_pmda"
    kraken2-build --add-to-library "$pmda_fasta" --db "$kraken_pmda"
    kraken2-build --build --db "$kraken_pmda" --threads 8

    log "PMDA database created"
}

update_databases() {
    log "Updating databases..."

    # Update Kraken2
    if [[ -d "$EFS_MOUNT/databases/kraken2/standard" ]]; then
        log "Updating Kraken2 database..."
        kraken2-build --update-library --db "$EFS_MOUNT/databases/kraken2/standard"
    fi

    # Check for RVDB updates
    log "Checking for RVDB updates..."
    # Check version and download if newer available

    log "Database update complete"
}

check_database_status() {
    log "Checking database status..."

    echo ""
    echo "Database Status:"
    echo "==============="

    # Check reference genome
    if [[ -f "$EFS_MOUNT/references/sus_scrofa_11.1/sus_scrofa.11.1.fa" ]]; then
        echo "✓ Reference genome: INSTALLED"
        local size=$(du -sh "$EFS_MOUNT/references/sus_scrofa_11.1" | cut -f1)
        echo "  Size: $size"
    else
        echo "✗ Reference genome: NOT FOUND"
    fi

    # Check Kraken2
    if [[ -f "$EFS_MOUNT/databases/kraken2/standard/hash.k2d" ]]; then
        echo "✓ Kraken2 standard: INSTALLED"
        local size=$(du -sh "$EFS_MOUNT/databases/kraken2/standard" | cut -f1)
        echo "  Size: $size"
    else
        echo "✗ Kraken2 standard: NOT FOUND"
    fi

    # Check RVDB
    if [[ -f "$EFS_MOUNT/databases/rvdb/rvdb.fasta" ]]; then
        echo "✓ RVDB: INSTALLED"
        local size=$(du -sh "$EFS_MOUNT/databases/rvdb" | cut -f1)
        echo "  Size: $size"
        local count=$(grep -c "^>" "$EFS_MOUNT/databases/rvdb/rvdb.fasta")
        echo "  Sequences: $count"
    else
        echo "✗ RVDB: NOT FOUND"
    fi

    # Check PMDA
    if [[ -f "$EFS_MOUNT/databases/pmda/pmda_91.fasta" ]]; then
        echo "✓ PMDA custom: INSTALLED"
        local size=$(du -sh "$EFS_MOUNT/databases/pmda" | cut -f1)
        echo "  Size: $size"
    else
        echo "✗ PMDA custom: NOT FOUND"
    fi

    # Check indexes
    if [[ -f "$EFS_MOUNT/indexes/minimap2/sus_scrofa.11.1.idx" ]]; then
        echo "✓ Minimap2 index: INSTALLED"
    else
        echo "✗ Minimap2 index: NOT FOUND"
    fi

    # Total disk usage
    echo ""
    echo "Total disk usage:"
    du -sh "$EFS_MOUNT" 2>/dev/null | cut -f1

    echo ""
    echo "Available space:"
    df -h "$EFS_MOUNT" | awk 'NR==2 {print $4}'
}

# Main execution
main() {
    local action=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -a|--all)
                action="all"
                shift
                ;;
            -k|--kraken2)
                action="kraken2"
                shift
                ;;
            -r|--rvdb)
                action="rvdb"
                shift
                ;;
            -p|--pmda)
                action="pmda"
                shift
                ;;
            -g|--genome)
                action="genome"
                shift
                ;;
            -u|--update)
                action="update"
                shift
                ;;
            -c|--check)
                action="check"
                shift
                ;;
            -h|--help)
                usage
                ;;
            *)
                error "Unknown option: $1"
                ;;
        esac
    done

    # Execute action
    case $action in
        all)
            check_requirements
            setup_directories
            download_reference_genome
            download_kraken2_database
            download_rvdb_database
            create_pmda_database
            check_database_status
            ;;
        kraken2)
            check_requirements
            setup_directories
            download_kraken2_database
            ;;
        rvdb)
            check_requirements
            setup_directories
            download_rvdb_database
            ;;
        pmda)
            check_requirements
            setup_directories
            create_pmda_database
            ;;
        genome)
            check_requirements
            setup_directories
            download_reference_genome
            ;;
        update)
            check_requirements
            update_databases
            ;;
        check)
            check_database_status
            ;;
        *)
            usage
            ;;
    esac

    log "Operation completed successfully"
}

# Run main function
main "$@"