#!/bin/bash
# MinION Metagenomics Pipeline - Install Python Environment
# Sets up Python 3.11 with bioinformatics packages

set -euo pipefail

PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
VENV_DIR="${VENV_DIR:-/opt/minion/venv}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

[[ $EUID -ne 0 ]] && log_error "Run as root"

# Install Python 3.11
log_info "Installing Python ${PYTHON_VERSION}..."
if command -v dnf &> /dev/null; then
    dnf install -y python3.11 python3.11-devel python3.11-pip
else
    apt-get update
    apt-get install -y software-properties-common
    add-apt-repository -y ppa:deadsnakes/ppa
    apt-get update
    apt-get install -y python3.11 python3.11-dev python3.11-venv python3-pip
fi

# Create virtual environment
log_info "Creating virtual environment at $VENV_DIR..."
mkdir -p "$(dirname "$VENV_DIR")"
python3.11 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# Upgrade pip
pip install --upgrade pip setuptools wheel

# Install core scientific packages
log_info "Installing core scientific packages..."
pip install \
    numpy==1.24.3 \
    pandas==2.0.3 \
    scipy==1.11.1 \
    scikit-learn==1.3.0 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    plotly==5.15.0 \
    jupyter==1.0.0 \
    ipython==8.14.0

# Install bioinformatics packages
log_info "Installing bioinformatics packages..."
pip install \
    biopython==1.81 \
    pysam==0.21.0 \
    pysamstats==1.1.2 \
    pybedtools==0.9.0 \
    pyvcf==0.6.8 \
    HTSeq==2.0.3

# Install Nanopore-specific packages
log_info "Installing Nanopore-specific packages..."
pip install \
    ont-fast5-api==4.1.2 \
    pod5==0.3.2 \
    nanoplot==1.42.0 \
    nanostat==1.6.0 \
    nanofilt==2.8.0 \
    pycoqc==2.5.2 \
    nanoget==1.19.1 \
    nanomath==1.3.0

# Install metagenomics packages
log_info "Installing metagenomics packages..."
pip install \
    kraken-biom==1.2.0 \
    metaphlan==4.1.0 \
    humann==3.8 \
    biom-format==2.1.15

# Install workflow and reporting packages
log_info "Installing workflow packages..."
pip install \
    snakemake==7.32.4 \
    nextflow==0.1.0 \
    multiqc==1.15 \
    jinja2==3.1.2 \
    reportlab==4.0.4 \
    xlsxwriter==3.1.2

# Install AWS and cloud packages
log_info "Installing cloud packages..."
pip install \
    boto3==1.28.25 \
    awscli==1.29.25 \
    s3fs==2023.6.0

# Install additional utilities
log_info "Installing utility packages..."
pip install \
    tqdm==4.65.0 \
    click==8.1.6 \
    pyyaml==6.0.1 \
    toml==0.10.2 \
    python-dotenv==1.0.0 \
    psutil==5.9.5 \
    gpustat==1.1.1

# Create activation script
cat > /etc/profile.d/minion-venv.sh << SCRIPT
#!/bin/bash
# Activate MinION virtual environment
alias activate-minion="source $VENV_DIR/bin/activate"
echo "To activate MinION Python environment: activate-minion"
SCRIPT

chmod +x /etc/profile.d/minion-venv.sh

# Create Python analysis script template
cat > "$VENV_DIR/analysis_template.py" << 'TEMPLATE'
#!/usr/bin/env python3
"""
MinION Metagenomics Analysis Template
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
import pysam

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="MinION Metagenomics Analysis"
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        type=Path,
        help='Input file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        type=Path,
        help='Output directory'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help='Number of threads'
    )
    return parser.parse_args()


def process_fastq(fastq_path):
    """Process FASTQ file."""
    logger.info(f"Processing {fastq_path}")

    records = []
    for record in SeqIO.parse(fastq_path, "fastq"):
        records.append({
            'id': record.id,
            'length': len(record.seq),
            'mean_quality': np.mean(record.letter_annotations["phred_quality"])
        })

    df = pd.DataFrame(records)
    logger.info(f"Processed {len(df)} reads")

    return df


def main():
    """Main analysis pipeline."""
    args = parse_arguments()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Process input
    if args.input.suffix in ['.fastq', '.fq']:
        df = process_fastq(args.input)

        # Generate statistics
        stats = {
            'total_reads': len(df),
            'total_bases': df['length'].sum(),
            'mean_length': df['length'].mean(),
            'median_length': df['length'].median(),
            'mean_quality': df['mean_quality'].mean()
        }

        # Save results
        output_file = args.output / 'analysis_results.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"Results saved to {output_file}")

        # Print summary
        for key, value in stats.items():
            print(f"{key}: {value:.2f}")

    logger.info("Analysis completed")


if __name__ == '__main__':
    main()
TEMPLATE

# Create Jupyter notebook configuration
log_info "Configuring Jupyter..."
jupyter notebook --generate-config -y || true
cat >> ~/.jupyter/jupyter_notebook_config.py << 'CONFIG'
c.NotebookApp.ip = '0.0.0.0'
c.NotebookApp.open_browser = False
c.NotebookApp.port = 8888
c.NotebookApp.allow_remote_access = True
CONFIG

# Test imports
log_info "Testing package imports..."
python3.11 -c "
import numpy, pandas, biopython, pysam, nanoplot
print('Core packages imported successfully')
"

log_info "======================================"
log_info "Python environment installed!"
log_info "======================================"
log_info "Python version: ${PYTHON_VERSION}"
log_info "Virtual environment: $VENV_DIR"
log_info ""
log_info "To activate environment:"
log_info "  source $VENV_DIR/bin/activate"
log_info "  or use: activate-minion"
log_info ""
log_info "Installed packages:"
log_info "  - Scientific: numpy, pandas, scipy, sklearn"
log_info "  - Bioinformatics: biopython, pysam"
log_info "  - Nanopore: nanoplot, pycoqc, pod5"
log_info "  - Metagenomics: kraken-biom, metaphlan"
log_info "  - Workflow: snakemake, multiqc"
log_info "  - Cloud: boto3, awscli"
log_info ""
log_info "Template script: $VENV_DIR/analysis_template.py"
log_info "======================================"