#!/bin/bash
# MinION Metagenomics Pipeline - Build Analysis AMI
# Creates custom AMI with all bioinformatics tools for analysis phases
# Base: Amazon Linux 2023

set -euo pipefail

# ===== Configuration =====
REGION="${AWS_REGION:-ap-northeast-1}"
INSTANCE_TYPE="${INSTANCE_TYPE:-m5.2xlarge}"
BASE_AMI_ID="${BASE_AMI_ID:-}"  # Will be auto-detected if not provided
AMI_NAME="minion-analysis-ami-$(date +%Y%m%d%H%M%S)"
KEY_NAME="${KEY_NAME:-}"  # Optional: SSH key pair name

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ===== Functions =====
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_requirements() {
    log_info "Checking requirements..."

    # Check AWS CLI
    if ! command -v aws &> /dev/null; then
        log_error "AWS CLI is not installed"
        exit 1
    fi

    # Check AWS credentials
    if ! aws sts get-caller-identity &> /dev/null; then
        log_error "AWS credentials not configured"
        exit 1
    fi

    # Check jq
    if ! command -v jq &> /dev/null; then
        log_warn "jq not installed, installing..."
        if [[ "$OSTYPE" == "darwin"* ]]; then
            brew install jq
        else
            sudo yum install -y jq || sudo apt-get install -y jq
        fi
    fi
}

get_base_ami() {
    if [ -n "$BASE_AMI_ID" ]; then
        echo "$BASE_AMI_ID"
        return
    fi

    log_info "Finding latest Amazon Linux 2023 AMI..."

    BASE_AMI_ID=$(aws ec2 describe-images \
        --region "$REGION" \
        --owners amazon \
        --filters \
            "Name=name,Values=al2023-ami-*-x86_64" \
            "Name=state,Values=available" \
        --query 'sort_by(Images, &CreationDate)[-1].ImageId' \
        --output text)

    if [ "$BASE_AMI_ID" == "None" ] || [ -z "$BASE_AMI_ID" ]; then
        log_error "Could not find Amazon Linux 2023 AMI"
        exit 1
    fi

    log_info "Using AMI: $BASE_AMI_ID"
    echo "$BASE_AMI_ID"
}

create_user_data_script() {
    cat > /tmp/analysis_ami_setup.sh << 'SETUP_SCRIPT'
#!/bin/bash
set -euo pipefail

# Log all output
exec > >(tee -a /var/log/ami-build.log)
exec 2>&1

echo "========================================="
echo "Starting Analysis AMI Setup"
echo "Date: $(date)"
echo "========================================="

# ===== System Updates and Base Packages =====
echo "Updating system packages..."
dnf update -y
dnf install -y \
    git \
    wget \
    curl \
    gcc \
    gcc-c++ \
    make \
    cmake \
    autoconf \
    automake \
    libtool \
    zlib-devel \
    bzip2-devel \
    xz-devel \
    libcurl-devel \
    openssl-devel \
    ncurses-devel \
    python3-devel \
    java-11-amazon-corretto \
    htop \
    screen \
    tmux \
    tree \
    pigz \
    parallel \
    awscli \
    amazon-efs-utils

# ===== Create directory structure =====
mkdir -p /opt/tools
mkdir -p /opt/databases
mkdir -p /data/fastq
mkdir -p /data/analysis
mkdir -p /data/reports
mkdir -p /data/logs
mkdir -p /mnt/efs

# Set environment variables
cat >> /etc/environment << 'ENV'
PATH="/opt/tools/bin:/usr/local/bin:/usr/bin:/bin"
KRAKEN2_DB_PATH="/mnt/efs/kraken2"
BLAST_DB_PATH="/mnt/efs/blast"
ENV

source /etc/environment

# ===== Install Python 3.11 and packages =====
echo "Setting up Python environment..."
dnf install -y python3.11 python3.11-pip python3.11-devel

# Create virtual environment
python3.11 -m venv /opt/tools/venv
source /opt/tools/venv/bin/activate

# Install Python packages
pip install --upgrade pip
pip install \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    biopython \
    pysam \
    pysamstats \
    matplotlib \
    seaborn \
    plotly \
    jupyter \
    boto3 \
    awscli \
    pycoqc \
    nanoplot \
    nanostat \
    nanofilt \
    multiqc

# ===== Install Conda (for complex dependencies) =====
echo "Installing Miniconda..."
cd /tmp
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/tools/miniconda3
rm Miniconda3-latest-Linux-x86_64.sh

/opt/tools/miniconda3/bin/conda init bash
source ~/.bashrc

# Add bioconda channels
/opt/tools/miniconda3/bin/conda config --add channels defaults
/opt/tools/miniconda3/bin/conda config --add channels bioconda
/opt/tools/miniconda3/bin/conda config --add channels conda-forge

# ===== Install QC Tools =====
echo "Installing QC tools..."

# FastQC
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod +x FastQC/fastqc
ln -s /opt/tools/FastQC/fastqc /usr/local/bin/
rm fastqc_v0.12.1.zip

# NanoPlot (already installed via pip)
# PycoQC (already installed via pip)

# ===== Install Alignment Tools =====
echo "Installing alignment tools..."

# Minimap2
cd /opt/tools
git clone https://github.com/lh3/minimap2
cd minimap2
make
cp minimap2 /usr/local/bin/
cd ..

# Samtools
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=/opt/tools
make && make install
ln -s /opt/tools/bin/samtools /usr/local/bin/
cd ..
rm samtools-1.18.tar.bz2

# BCFtools
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
tar -xjf bcftools-1.18.tar.bz2
cd bcftools-1.18
./configure --prefix=/opt/tools
make && make install
ln -s /opt/tools/bin/bcftools /usr/local/bin/
cd ..
rm bcftools-1.18.tar.bz2

# BWA (for short read alignment if needed)
cd /opt/tools
git clone https://github.com/lh3/bwa.git
cd bwa
make
cp bwa /usr/local/bin/
cd ..

# ===== Install Kraken2 and Bracken =====
echo "Installing Kraken2..."
cd /opt/tools
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh /opt/tools/kraken2_install
ln -s /opt/tools/kraken2_install/kraken2 /usr/local/bin/
ln -s /opt/tools/kraken2_install/kraken2-build /usr/local/bin/
ln -s /opt/tools/kraken2_install/kraken2-inspect /usr/local/bin/
cd ..

# Bracken
cd /opt/tools
git clone https://github.com/jenniferlu717/Bracken.git
cd Bracken
./install_bracken.sh
ln -s /opt/tools/Bracken/bracken /usr/local/bin/
ln -s /opt/tools/Bracken/bracken-build /usr/local/bin/
cd ..

# KrakenTools
cd /opt/tools
git clone https://github.com/jenniferlu717/KrakenTools.git

# ===== Install BLAST+ =====
echo "Installing BLAST+..."
cd /opt/tools
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.15.0+-x64-linux.tar.gz
cp ncbi-blast-2.15.0+/bin/* /usr/local/bin/
rm -rf ncbi-blast-2.15.0+-x64-linux.tar.gz
rm -rf ncbi-blast-2.15.0+

# ===== Install Diamond (faster protein alignment) =====
echo "Installing Diamond..."
cd /opt/tools
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
mv diamond /usr/local/bin/
rm diamond-linux64.tar.gz

# ===== Install Assembly Tools =====
echo "Installing assembly tools..."

# Flye (long-read assembler)
/opt/tools/miniconda3/bin/conda install -y -c bioconda flye

# Canu (long-read assembler)
cd /opt/tools
wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz
tar -xJf canu-2.2.Linux-amd64.tar.xz
ln -s /opt/tools/canu-2.2/bin/canu /usr/local/bin/
rm canu-2.2.Linux-amd64.tar.xz

# MEGAHIT (for metagenome assembly)
cd /opt/tools
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar -xzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
ln -s /opt/tools/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit /usr/local/bin/
rm MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

# SPAdes (hybrid assembler)
cd /opt/tools
wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar -xzf SPAdes-3.15.5-Linux.tar.gz
ln -s /opt/tools/SPAdes-3.15.5-Linux/bin/spades.py /usr/local/bin/
ln -s /opt/tools/SPAdes-3.15.5-Linux/bin/metaspades.py /usr/local/bin/
rm SPAdes-3.15.5-Linux.tar.gz

# ===== Install Additional Metagenomics Tools =====
echo "Installing additional metagenomics tools..."

# MetaPhlAn
/opt/tools/miniconda3/bin/conda install -y -c bioconda metaphlan

# Centrifuge
cd /opt/tools
wget https://github.com/DaehwanKimLab/centrifuge/archive/v1.0.4.tar.gz
tar -xzf v1.0.4.tar.gz
cd centrifuge-1.0.4
make
make install prefix=/opt/tools
ln -s /opt/tools/bin/centrifuge /usr/local/bin/
cd ..
rm v1.0.4.tar.gz

# ===== Install Visualization Tools =====
echo "Installing visualization tools..."

# IGV (for local visualization if needed)
cd /opt/tools
wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.2_WithJava.zip
unzip IGV_Linux_2.16.2_WithJava.zip
rm IGV_Linux_2.16.2_WithJava.zip

# ===== Create Pipeline Scripts =====
echo "Creating pipeline utility scripts..."

# QC script
cat > /usr/local/bin/run-qc.sh << 'QC_SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_DIR"
    exit 1
}

while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "${INPUT:-}" ] || [ -z "${OUTPUT:-}" ]; then
    usage
fi

mkdir -p "$OUTPUT"

echo "Running QC analysis..."

# NanoPlot
NanoPlot --fastq "$INPUT" \
    --outdir "$OUTPUT/nanoplot" \
    --prefix sample_ \
    --plots hex dot pauvre kde \
    --N50 \
    --threads 8

# NanoStat
NanoStat --fastq "$INPUT" > "$OUTPUT/nanostats.txt"

# PycoQC (if sequencing summary is available)
# pycoqc -f "$SUMMARY" -o "$OUTPUT/pycoqc_report.html"

echo "QC analysis completed!"
QC_SCRIPT

chmod +x /usr/local/bin/run-qc.sh

# Host removal script
cat > /usr/local/bin/run-host-removal.sh << 'HOST_SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -r REFERENCE -o OUTPUT_DIR"
    exit 1
}

while getopts "i:r:o:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "${INPUT:-}" ] || [ -z "${REFERENCE:-}" ] || [ -z "${OUTPUT:-}" ]; then
    usage
fi

mkdir -p "$OUTPUT"

echo "Removing host sequences..."

# Align to host genome
minimap2 -ax map-ont "$REFERENCE" "$INPUT" \
    -t 8 \
    --secondary=no \
    | samtools view -bS - > "$OUTPUT/aligned.bam"

# Extract unmapped reads (non-host)
samtools view -b -f 4 "$OUTPUT/aligned.bam" > "$OUTPUT/unmapped.bam"
samtools fastq "$OUTPUT/unmapped.bam" | pigz > "$OUTPUT/non_host.fastq.gz"

# Calculate statistics
TOTAL=$(samtools view -c "$OUTPUT/aligned.bam")
UNMAPPED=$(samtools view -c -f 4 "$OUTPUT/aligned.bam")
MAPPED=$((TOTAL - UNMAPPED))

echo "Total reads: $TOTAL"
echo "Host reads: $MAPPED"
echo "Non-host reads: $UNMAPPED"
echo "Host removal rate: $(echo "scale=2; $MAPPED * 100 / $TOTAL" | bc)%"

echo "Host removal completed!"
HOST_SCRIPT

chmod +x /usr/local/bin/run-host-removal.sh

# Kraken2 script
cat > /usr/local/bin/run-kraken2.sh << 'KRAKEN_SCRIPT'
#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -d DATABASE -o OUTPUT_DIR"
    exit 1
}

while getopts "i:d:o:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        d) DATABASE="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "${INPUT:-}" ] || [ -z "${DATABASE:-}" ] || [ -z "${OUTPUT:-}" ]; then
    usage
fi

mkdir -p "$OUTPUT"

echo "Running Kraken2 classification..."

kraken2 --db "$DATABASE" \
    --threads 8 \
    --report "$OUTPUT/kraken2_report.txt" \
    --output "$OUTPUT/kraken2_output.txt" \
    --confidence 0.05 \
    --report-zero-counts \
    --use-names \
    "$INPUT"

# Run Bracken for abundance estimation
bracken -d "$DATABASE" \
    -i "$OUTPUT/kraken2_report.txt" \
    -o "$OUTPUT/bracken_species.txt" \
    -r 150 \
    -l S \
    -t 10

echo "Kraken2 analysis completed!"
KRAKEN_SCRIPT

chmod +x /usr/local/bin/run-kraken2.sh

# ===== CloudWatch Agent Configuration =====
echo "Configuring CloudWatch agent..."
cat > /opt/aws/amazon-cloudwatch-agent/etc/cloudwatch-config.json << 'CW_CONFIG'
{
  "metrics": {
    "namespace": "MinION/Analysis",
    "metrics_collected": {
      "cpu": {
        "measurement": [
          {
            "name": "cpu_usage_idle",
            "rename": "CPU_IDLE",
            "unit": "Percent"
          }
        ],
        "metrics_collection_interval": 60
      },
      "mem": {
        "measurement": [
          {
            "name": "mem_used_percent",
            "rename": "Memory_Utilization",
            "unit": "Percent"
          }
        ],
        "metrics_collection_interval": 60
      },
      "disk": {
        "measurement": [
          {
            "name": "used_percent",
            "rename": "Disk_Utilization",
            "unit": "Percent"
          }
        ],
        "resources": [
          "/data"
        ],
        "metrics_collection_interval": 60
      }
    }
  },
  "logs": {
    "logs_collected": {
      "files": {
        "collect_list": [
          {
            "file_path": "/data/logs/*.log",
            "log_group_name": "/aws/ec2/minion/analysis",
            "log_stream_name": "{instance_id}/analysis"
          }
        ]
      }
    }
  }
}
CW_CONFIG

# ===== Set permissions =====
chown -R ec2-user:ec2-user /data
chown -R ec2-user:ec2-user /opt/tools

# ===== Create MOTD =====
cat > /etc/motd << 'MOTD'
========================================
    MinION Analysis AMI
========================================
Tools installed:
  - QC: FastQC, NanoPlot, PycoQC
  - Alignment: Minimap2, Samtools, BWA
  - Classification: Kraken2, Bracken, Centrifuge
  - Assembly: Flye, Canu, MEGAHIT, SPAdes
  - Search: BLAST+, Diamond
  - Python: 3.11 with bioinformatics packages

Data directory: /data/

Quick commands:
  - run-qc.sh          : Run QC analysis
  - run-host-removal.sh: Remove host sequences
  - run-kraken2.sh     : Run taxonomic classification

========================================
MOTD

echo "AMI setup completed successfully!"
echo "System will reboot in 10 seconds..."
sleep 10
# reboot
SETUP_SCRIPT

    chmod +x /tmp/analysis_ami_setup.sh
}

launch_instance() {
    local ami_id=$1

    log_info "Launching EC2 instance..."

    # Use existing security group or create new one
    SG_ID=$(aws ec2 describe-security-groups \
        --region "$REGION" \
        --group-names "ami-builder-sg" \
        --query 'SecurityGroups[0].GroupId' \
        --output text 2>/dev/null || echo "")

    if [ -z "$SG_ID" ] || [ "$SG_ID" == "None" ]; then
        log_info "Creating security group..."
        SG_ID=$(aws ec2 create-security-group \
            --region "$REGION" \
            --group-name "ami-builder-sg" \
            --description "Temporary SG for AMI building" \
            --query 'GroupId' \
            --output text)

        if [ -n "$KEY_NAME" ]; then
            aws ec2 authorize-security-group-ingress \
                --region "$REGION" \
                --group-id "$SG_ID" \
                --protocol tcp \
                --port 22 \
                --cidr 0.0.0.0/0
        fi
    fi

    # Launch parameters
    LAUNCH_PARAMS=(
        --region "$REGION"
        --image-id "$ami_id"
        --instance-type "$INSTANCE_TYPE"
        --security-group-ids "$SG_ID"
        --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$AMI_NAME-builder}]"
        --user-data file:///tmp/analysis_ami_setup.sh
        --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":100,\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]"
    )

    if [ -n "$KEY_NAME" ]; then
        LAUNCH_PARAMS+=(--key-name "$KEY_NAME")
    fi

    INSTANCE_ID=$(aws ec2 run-instances "${LAUNCH_PARAMS[@]}" \
        --query 'Instances[0].InstanceId' \
        --output text)

    log_info "Instance launched: $INSTANCE_ID"
    echo "$INSTANCE_ID"
}

wait_for_instance() {
    local instance_id=$1

    log_info "Waiting for instance to be ready..."

    aws ec2 wait instance-running \
        --region "$REGION" \
        --instance-ids "$instance_id"

    log_info "Instance is running"

    log_info "Waiting for status checks to pass (this may take 10-15 minutes)..."
    aws ec2 wait instance-status-ok \
        --region "$REGION" \
        --instance-ids "$instance_id" \
        --max-attempts 40 || true

    log_info "Waiting for user data script to complete..."
    sleep 600  # 10 minutes for comprehensive software installation

    log_info "Instance should be ready"
}

create_ami() {
    local instance_id=$1

    log_info "Creating AMI from instance $instance_id..."

    AMI_ID=$(aws ec2 create-image \
        --region "$REGION" \
        --instance-id "$instance_id" \
        --name "$AMI_NAME" \
        --description "MinION Analysis AMI with bioinformatics tools" \
        --tag-specifications "ResourceType=image,Tags=[{Key=Name,Value=$AMI_NAME},{Key=Purpose,Value=MinION-Analysis}]" \
        --query 'ImageId' \
        --output text)

    log_info "AMI creation initiated: $AMI_ID"

    log_info "Waiting for AMI to be available..."
    aws ec2 wait image-available \
        --region "$REGION" \
        --image-ids "$AMI_ID"

    log_info "AMI created successfully: $AMI_ID"
    echo "$AMI_ID"
}

cleanup() {
    local instance_id=$1

    log_info "Cleaning up..."

    if [ -n "$instance_id" ]; then
        log_info "Terminating instance $instance_id..."
        aws ec2 terminate-instances \
            --region "$REGION" \
            --instance-ids "$instance_id"

        aws ec2 wait instance-terminated \
            --region "$REGION" \
            --instance-ids "$instance_id"
    fi

    log_info "Cleanup completed"
}

# ===== Main Execution =====
main() {
    log_info "Starting Analysis AMI build process"

    check_requirements

    BASE_AMI=$(get_base_ami)

    create_user_data_script

    INSTANCE_ID=$(launch_instance "$BASE_AMI")

    trap "cleanup $INSTANCE_ID" EXIT

    wait_for_instance "$INSTANCE_ID"

    FINAL_AMI_ID=$(create_ami "$INSTANCE_ID")

    cat << EOF

========================================
AMI Build Completed Successfully!
========================================
AMI ID: $FINAL_AMI_ID
AMI Name: $AMI_NAME
Region: $REGION

To use this AMI in Terraform:
  analysis_ami_id = "$FINAL_AMI_ID"

To launch an instance manually:
  aws ec2 run-instances \\
    --image-id $FINAL_AMI_ID \\
    --instance-type m5.xlarge \\
    --key-name YOUR_KEY_NAME

========================================
EOF
}

main "$@"