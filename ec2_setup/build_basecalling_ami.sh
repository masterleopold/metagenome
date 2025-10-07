#!/bin/bash
# MinION Metagenomics Pipeline - Build Basecalling AMI
# Creates custom AMI with GPU drivers, CUDA, and Dorado for basecalling
# Base: AWS Deep Learning AMI or Amazon Linux 2023 with GPU

set -euo pipefail

# ===== Configuration =====
REGION="${AWS_REGION:-ap-northeast-1}"
INSTANCE_TYPE="${INSTANCE_TYPE:-g4dn.xlarge}"
BASE_AMI_ID="${BASE_AMI_ID:-}"  # Will be auto-detected if not provided
AMI_NAME="minion-basecalling-ami-$(date +%Y%m%d%H%M%S)"
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

    # Check jq for JSON parsing
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

    log_info "Finding latest Deep Learning AMI..."

    # Find the latest AWS Deep Learning AMI with GPU support
    BASE_AMI_ID=$(aws ec2 describe-images \
        --region "$REGION" \
        --owners amazon \
        --filters \
            "Name=name,Values=Deep Learning AMI GPU PyTorch * (Amazon Linux 2023)*" \
            "Name=state,Values=available" \
        --query 'sort_by(Images, &CreationDate)[-1].ImageId' \
        --output text)

    if [ "$BASE_AMI_ID" == "None" ] || [ -z "$BASE_AMI_ID" ]; then
        log_error "Could not find Deep Learning AMI"
        exit 1
    fi

    log_info "Using AMI: $BASE_AMI_ID"
    echo "$BASE_AMI_ID"
}

create_user_data_script() {
    cat > /tmp/basecalling_ami_setup.sh << 'SETUP_SCRIPT'
#!/bin/bash
set -euo pipefail

# Log all output
exec > >(tee -a /var/log/ami-build.log)
exec 2>&1

echo "========================================="
echo "Starting Basecalling AMI Setup"
echo "Date: $(date)"
echo "========================================="

# ===== System Updates =====
echo "Updating system packages..."
dnf update -y
dnf install -y \
    git \
    wget \
    curl \
    htop \
    nvtop \
    screen \
    tmux \
    tree \
    pigz \
    parallel \
    awscli \
    amazon-efs-utils

# ===== NVIDIA Driver and CUDA Check =====
echo "Checking NVIDIA drivers..."
if ! nvidia-smi &> /dev/null; then
    echo "NVIDIA drivers not found, installing..."
    dnf install -y nvidia-driver-latest-dkms cuda-toolkit-12-2
fi

nvidia-smi
nvcc --version || echo "CUDA compiler not found"

# ===== Create directory structure =====
mkdir -p /opt/ont
mkdir -p /opt/ont/dorado
mkdir -p /opt/ont/models
mkdir -p /data/fast5
mkdir -p /data/fastq
mkdir -p /data/logs
mkdir -p /mnt/efs

# ===== Install Dorado =====
echo "Installing Dorado..."
cd /opt/ont/dorado

# Download latest Dorado release
DORADO_VERSION=$(curl -s https://api.github.com/repos/nanoporetech/dorado/releases/latest | grep tag_name | cut -d'"' -f4)
echo "Installing Dorado version: $DORADO_VERSION"

# Download Linux x86_64 binary
wget -q "https://github.com/nanoporetech/dorado/releases/download/${DORADO_VERSION}/dorado-${DORADO_VERSION}-linux-x64.tar.gz"
tar -xzf "dorado-${DORADO_VERSION}-linux-x64.tar.gz"
rm "dorado-${DORADO_VERSION}-linux-x64.tar.gz"

# Create symlink for easy access
ln -sf "/opt/ont/dorado/dorado-${DORADO_VERSION}-linux-x64/bin/dorado" /usr/local/bin/dorado

# Verify installation
dorado --version

# ===== Download Dorado Models =====
echo "Downloading Dorado models..."
cd /opt/ont/models

# Download commonly used models for MinION
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.3.0
dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.3.0
dorado download --model dna_r10.4.1_e8.2_400bps_fast@v4.3.0

# For RNA if needed
# dorado download --model rna004_130bps_sup@v3.0.1

echo "Models downloaded to /opt/ont/models"
ls -la /opt/ont/models/

# ===== Install Pod5 Tools =====
echo "Installing pod5 tools..."
pip3 install pod5

# ===== Install MinKNOW prerequisites (optional) =====
# Only if you want to run MinKNOW on the instance
# echo "Installing MinKNOW prerequisites..."
# wget -q https://cdn.oxfordnanoportal.com/apt/ont-repo.pub -O /tmp/ont-repo.pub
# apt-key add /tmp/ont-repo.pub
# echo "deb http://cdn.oxfordnanoportal.com/apt focal-stable non-free" > /etc/apt/sources.list.d/ont.list
# apt-get update
# apt-get install -y ont-minknow-core-minion-nc

# ===== Install Python analysis tools =====
echo "Installing Python analysis tools..."
pip3 install --upgrade pip
pip3 install \
    numpy \
    pandas \
    biopython \
    pysam \
    h5py \
    matplotlib \
    seaborn \
    plotly \
    boto3 \
    awscli

# ===== Install monitoring tools =====
echo "Installing monitoring tools..."
pip3 install \
    gpustat \
    py3nvml

# ===== Create monitoring script =====
cat > /usr/local/bin/gpu-monitor.sh << 'MONITOR'
#!/bin/bash
while true; do
    clear
    echo "========== GPU Status =========="
    nvidia-smi
    echo ""
    echo "========== Dorado Processes =========="
    ps aux | grep dorado | grep -v grep
    echo ""
    echo "========== Disk Usage =========="
    df -h /data
    echo ""
    echo "Press Ctrl+C to exit"
    sleep 5
done
MONITOR

chmod +x /usr/local/bin/gpu-monitor.sh

# ===== Create basecalling wrapper script =====
cat > /usr/local/bin/run-basecalling.sh << 'WRAPPER'
#!/bin/bash
set -euo pipefail

# Usage function
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-m MODEL] [-d DEVICE] [-b BATCH_SIZE]"
    echo "  -i INPUT_DIR     Directory containing FAST5/POD5 files"
    echo "  -o OUTPUT_DIR    Directory for output FASTQ files"
    echo "  -m MODEL         Model to use (default: dna_r10.4.1_e8.2_400bps_sup@v4.3.0)"
    echo "  -d DEVICE        CUDA device (default: cuda:0)"
    echo "  -b BATCH_SIZE    Batch size (default: 256)"
    echo "  -x               Enable duplex basecalling"
    exit 1
}

# Default values
MODEL="dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
DEVICE="cuda:0"
BATCH_SIZE=256
DUPLEX=false

# Parse arguments
while getopts "i:o:m:d:b:xh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        d) DEVICE="$OPTARG" ;;
        b) BATCH_SIZE="$OPTARG" ;;
        x) DUPLEX=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [ -z "${INPUT_DIR:-}" ] || [ -z "${OUTPUT_DIR:-}" ]; then
    echo "Error: Input and output directories are required"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="$OUTPUT_DIR/basecalling_$(date +%Y%m%d_%H%M%S).log"

# Start basecalling
echo "Starting basecalling..." | tee "$LOG_FILE"
echo "Input: $INPUT_DIR" | tee -a "$LOG_FILE"
echo "Output: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Model: $MODEL" | tee -a "$LOG_FILE"
echo "Device: $DEVICE" | tee -a "$LOG_FILE"
echo "Batch size: $BATCH_SIZE" | tee -a "$LOG_FILE"
echo "Duplex: $DUPLEX" | tee -a "$LOG_FILE"

# Run Dorado
if [ "$DUPLEX" = true ]; then
    dorado duplex \
        /opt/ont/models/"$MODEL" \
        "$INPUT_DIR" \
        --device "$DEVICE" \
        --batchsize "$BATCH_SIZE" \
        --emit-fastq \
        2>&1 | tee -a "$LOG_FILE" | \
        pigz -p 4 > "$OUTPUT_DIR/basecalls_duplex.fastq.gz"
else
    dorado basecaller \
        /opt/ont/models/"$MODEL" \
        "$INPUT_DIR" \
        --device "$DEVICE" \
        --batchsize "$BATCH_SIZE" \
        --emit-fastq \
        2>&1 | tee -a "$LOG_FILE" | \
        pigz -p 4 > "$OUTPUT_DIR/basecalls.fastq.gz"
fi

echo "Basecalling completed!" | tee -a "$LOG_FILE"

# Generate summary
echo "Generating summary..." | tee -a "$LOG_FILE"
if [ "$DUPLEX" = true ]; then
    dorado summary "$OUTPUT_DIR/basecalls_duplex.fastq.gz" > "$OUTPUT_DIR/sequencing_summary.txt"
else
    dorado summary "$OUTPUT_DIR/basecalls.fastq.gz" > "$OUTPUT_DIR/sequencing_summary.txt"
fi

echo "Summary saved to $OUTPUT_DIR/sequencing_summary.txt" | tee -a "$LOG_FILE"
WRAPPER

chmod +x /usr/local/bin/run-basecalling.sh

# ===== CloudWatch Agent Configuration =====
echo "Configuring CloudWatch agent..."
cat > /opt/aws/amazon-cloudwatch-agent/etc/cloudwatch-config.json << 'CW_CONFIG'
{
  "metrics": {
    "namespace": "MinION/Basecalling",
    "metrics_collected": {
      "nvidia_smi": {
        "measurement": [
          {
            "name": "utilization.gpu",
            "rename": "GPU_Utilization",
            "unit": "Percent"
          },
          {
            "name": "utilization.memory",
            "rename": "GPU_Memory_Utilization",
            "unit": "Percent"
          },
          {
            "name": "temperature.gpu",
            "rename": "GPU_Temperature",
            "unit": "None"
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
            "log_group_name": "/aws/ec2/minion/basecalling",
            "log_stream_name": "{instance_id}/basecalling"
          }
        ]
      }
    }
  }
}
CW_CONFIG

# Start CloudWatch agent
# /opt/aws/amazon-cloudwatch-agent/bin/amazon-cloudwatch-agent-ctl \
#     -a fetch-config \
#     -m ec2 \
#     -s \
#     -c file:/opt/aws/amazon-cloudwatch-agent/etc/cloudwatch-config.json

# ===== Create startup script =====
cat > /etc/systemd/system/basecalling-setup.service << 'SYSTEMD'
[Unit]
Description=MinION Basecalling Setup
After=network.target

[Service]
Type=oneshot
ExecStart=/bin/bash -c 'nvidia-smi && echo "GPU Ready"'
RemainAfterExit=true

[Install]
WantedBy=multi-user.target
SYSTEMD

systemctl enable basecalling-setup.service

# ===== Set permissions =====
chown -R ec2-user:ec2-user /data
chown -R ec2-user:ec2-user /opt/ont

# ===== Final setup message =====
cat > /etc/motd << 'MOTD'
========================================
    MinION Basecalling AMI
========================================
Dorado is installed at: /usr/local/bin/dorado
Models are located at: /opt/ont/models/
Data directory: /data/

Quick commands:
  - run-basecalling.sh   : Run basecalling pipeline
  - gpu-monitor.sh       : Monitor GPU usage
  - dorado --help        : Dorado help
  - nvidia-smi           : Check GPU status

For duplex basecalling:
  run-basecalling.sh -i /data/fast5 -o /data/fastq -x

========================================
MOTD

echo "AMI setup completed successfully!"
echo "System will reboot in 10 seconds..."
sleep 10
# reboot
SETUP_SCRIPT

    chmod +x /tmp/basecalling_ami_setup.sh
}

launch_instance() {
    local ami_id=$1

    log_info "Launching EC2 instance..."

    # Create security group if needed
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

        # Allow SSH if key is provided
        if [ -n "$KEY_NAME" ]; then
            aws ec2 authorize-security-group-ingress \
                --region "$REGION" \
                --group-id "$SG_ID" \
                --protocol tcp \
                --port 22 \
                --cidr 0.0.0.0/0
        fi
    fi

    # Prepare launch parameters
    LAUNCH_PARAMS=(
        --region "$REGION"
        --image-id "$ami_id"
        --instance-type "$INSTANCE_TYPE"
        --security-group-ids "$SG_ID"
        --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$AMI_NAME-builder}]"
        --user-data file:///tmp/basecalling_ami_setup.sh
        --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":200,\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]"
    )

    # Add key pair if provided
    if [ -n "$KEY_NAME" ]; then
        LAUNCH_PARAMS+=(--key-name "$KEY_NAME")
    fi

    # Launch instance
    INSTANCE_ID=$(aws ec2 run-instances "${LAUNCH_PARAMS[@]}" \
        --query 'Instances[0].InstanceId' \
        --output text)

    log_info "Instance launched: $INSTANCE_ID"
    echo "$INSTANCE_ID"
}

wait_for_instance() {
    local instance_id=$1

    log_info "Waiting for instance to be ready..."

    # Wait for instance to be running
    aws ec2 wait instance-running \
        --region "$REGION" \
        --instance-ids "$instance_id"

    log_info "Instance is running"

    # Wait for status checks to pass (this takes longer but ensures the instance is fully ready)
    log_info "Waiting for status checks to pass (this may take 10-15 minutes)..."
    aws ec2 wait instance-status-ok \
        --region "$REGION" \
        --instance-ids "$instance_id" \
        --max-attempts 40 || true  # Don't fail if this times out

    # Additional wait for user data to complete
    log_info "Waiting for user data script to complete..."
    sleep 300  # 5 minutes for script to run

    log_info "Instance should be ready"
}

create_ami() {
    local instance_id=$1

    log_info "Creating AMI from instance $instance_id..."

    AMI_ID=$(aws ec2 create-image \
        --region "$REGION" \
        --instance-id "$instance_id" \
        --name "$AMI_NAME" \
        --description "MinION Basecalling AMI with Dorado and GPU support" \
        --tag-specifications "ResourceType=image,Tags=[{Key=Name,Value=$AMI_NAME},{Key=Purpose,Value=MinION-Basecalling}]" \
        --query 'ImageId' \
        --output text)

    log_info "AMI creation initiated: $AMI_ID"

    # Wait for AMI to be available
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

    # Terminate instance
    if [ -n "$instance_id" ]; then
        log_info "Terminating instance $instance_id..."
        aws ec2 terminate-instances \
            --region "$REGION" \
            --instance-ids "$instance_id"

        aws ec2 wait instance-terminated \
            --region "$REGION" \
            --instance-ids "$instance_id"
    fi

    # Note: Not deleting security group as it might be reused
    log_info "Cleanup completed"
}

# ===== Main Execution =====
main() {
    log_info "Starting Basecalling AMI build process"

    # Check requirements
    check_requirements

    # Get base AMI
    BASE_AMI=$(get_base_ami)

    # Create user data script
    create_user_data_script

    # Launch instance
    INSTANCE_ID=$(launch_instance "$BASE_AMI")

    # Set trap to cleanup on exit
    trap "cleanup $INSTANCE_ID" EXIT

    # Wait for instance to be ready
    wait_for_instance "$INSTANCE_ID"

    # Create AMI
    FINAL_AMI_ID=$(create_ami "$INSTANCE_ID")

    # Output results
    cat << EOF

========================================
AMI Build Completed Successfully!
========================================
AMI ID: $FINAL_AMI_ID
AMI Name: $AMI_NAME
Region: $REGION

To use this AMI in Terraform:
  basecalling_ami_id = "$FINAL_AMI_ID"

To launch an instance manually:
  aws ec2 run-instances \\
    --image-id $FINAL_AMI_ID \\
    --instance-type g4dn.xlarge \\
    --key-name YOUR_KEY_NAME

========================================
EOF
}

# Run main function
main "$@"