#!/bin/bash
# MinION Metagenomics Pipeline - Phase 1: Basecalling
# Performs duplex basecalling with Dorado for Q30 accuracy

set -euo pipefail

# ===== Configuration =====
SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-/data/logs}"
WORK_DIR="${WORK_DIR:-/data/work}"

# Dorado settings
DORADO_BIN="${DORADO_BIN:-dorado}"
MODEL_PATH="${MODEL_PATH:-/opt/ont/models}"
DEFAULT_MODEL="dna_r10.4.1_e8.2_400bps_sup@v4.3.0"

# AWS settings
S3_RAW_BUCKET="${S3_RAW_BUCKET:-}"
S3_ANALYSIS_BUCKET="${S3_ANALYSIS_BUCKET:-}"

# Database settings
RUN_ID="${RUN_ID:-}"
DB_HOST="${DB_HOST:-}"
DB_NAME="${DB_NAME:-minion_metadata}"
DB_USER="${DB_USER:-minion_user}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# ===== Functions =====
log_info() {
    echo -e "${GREEN}[INFO $(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

log_warn() {
    echo -e "${YELLOW}[WARN $(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[ERROR $(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
    update_status "failed" "$1"
    exit 1
}

usage() {
    cat << EOF
Usage: $SCRIPT_NAME [OPTIONS]

Phase 1: Basecalling with Dorado (Duplex mode for Q30 accuracy)

Required Options:
    -i INPUT_DIR      Directory containing FAST5/POD5 files
    -o OUTPUT_DIR     Directory for output FASTQ files
    -r RUN_ID         Unique run identifier

Optional Options:
    -m MODEL          Dorado model (default: $DEFAULT_MODEL)
    -d DEVICE         CUDA device (default: cuda:all, use 'cpu' for CPU mode)
    -b BATCH_SIZE     Batch size (default: auto)
    -q MIN_QSCORE     Minimum quality score filter (default: 9)
    -t TRIM           Trim strategy (none/primers/adapters/all)
    -s S3_UPLOAD      Upload results to S3 (yes/no, default: yes)
    -x                Skip duplex calling (faster but lower accuracy)
    -h                Show this help message

Examples:
    # Basic duplex basecalling
    $SCRIPT_NAME -i /data/fast5 -o /data/fastq -r RUN001

    # CPU mode with quality filtering
    $SCRIPT_NAME -i /data/fast5 -o /data/fastq -r RUN001 -d cpu -q 10

    # Skip duplex for faster processing
    $SCRIPT_NAME -i /data/fast5 -o /data/fastq -r RUN001 -x

EOF
    exit 1
}

check_requirements() {
    log_info "Checking requirements..."

    # Check Dorado
    if ! command -v "$DORADO_BIN" &> /dev/null; then
        log_error "Dorado not found. Please install Dorado first."
    fi

    # Check CUDA if not CPU mode
    if [[ "${DEVICE:-cuda:all}" != "cpu" ]]; then
        if ! nvidia-smi &> /dev/null; then
            log_warn "NVIDIA GPU not detected. Switching to CPU mode."
            DEVICE="cpu"
        fi
    fi

    # Check model
    if [[ "$MODEL" == *"@"* ]]; then
        MODEL_DIR="$MODEL_PATH/$MODEL"
    else
        MODEL_DIR="$MODEL_PATH/$MODEL"
    fi

    if [[ ! -d "$MODEL_DIR" ]]; then
        log_warn "Model not found locally. Will attempt to download."
        $DORADO_BIN download --model "$MODEL" --directory "$MODEL_PATH" || \
            log_error "Failed to download model: $MODEL"
    fi

    # Check AWS CLI if S3 upload enabled
    if [[ "${S3_UPLOAD}" == "yes" ]]; then
        if ! command -v aws &> /dev/null; then
            log_error "AWS CLI not found but S3 upload requested"
        fi
    fi
}

update_status() {
    local status=$1
    local message="${2:-}"

    if [[ -n "$DB_HOST" ]] && [[ -n "$RUN_ID" ]]; then
        export PGPASSWORD="${DB_PASSWORD:-}"
        psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c \
            "UPDATE workflow_executions
             SET status = '$status',
                 updated_at = CURRENT_TIMESTAMP
                 ${message:+, error_message = '$message'}
             WHERE run_id = '$RUN_ID';" 2>/dev/null || true
    fi

    # Send SNS notification for failures
    if [[ "$status" == "failed" ]] && [[ -n "${SNS_TOPIC_ARN:-}" ]]; then
        aws sns publish \
            --topic-arn "$SNS_TOPIC_ARN" \
            --subject "Basecalling Failed - Run $RUN_ID" \
            --message "Basecalling failed for run $RUN_ID: $message" 2>/dev/null || true
    fi
}

perform_basecalling() {
    local input_dir=$1
    local output_dir=$2

    log_info "Starting basecalling..."
    log_info "Input: $input_dir"
    log_info "Output: $output_dir"
    log_info "Model: $MODEL"
    log_info "Device: $DEVICE"
    log_info "Duplex: $([ "$SKIP_DUPLEX" == "true" ] && echo "No" || echo "Yes")"

    # Build Dorado command
    if [[ "$SKIP_DUPLEX" == "true" ]]; then
        CMD="$DORADO_BIN basecaller"
    else
        CMD="$DORADO_BIN duplex"
    fi

    CMD="$CMD $MODEL_DIR $input_dir"
    CMD="$CMD --device $DEVICE"

    # Add optional parameters
    [[ -n "${BATCH_SIZE:-}" ]] && [[ "$BATCH_SIZE" != "auto" ]] && \
        CMD="$CMD --batchsize $BATCH_SIZE"

    [[ -n "${MIN_QSCORE:-}" ]] && \
        CMD="$CMD --min-qscore $MIN_QSCORE"

    [[ -n "${TRIM:-}" ]] && \
        CMD="$CMD --trim $TRIM"

    # Run basecalling
    local start_time=$(date +%s)
    update_status "basecalling"

    if [[ "$SKIP_DUPLEX" == "true" ]]; then
        OUTPUT_FILE="$output_dir/basecalls_simplex.fastq"
    else
        OUTPUT_FILE="$output_dir/basecalls_duplex.fastq"
    fi

    log_info "Executing: $CMD"

    # Run with progress monitoring
    $CMD --emit-fastq 2>&1 | tee "$output_dir/dorado.log" | \
        grep -E "Basecalled|samples/s|Duplex" | \
        while read line; do
            echo "  $line"
        done > "$OUTPUT_FILE"

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))

    log_info "Basecalling completed in $((duration/60)) minutes"

    # Compress output
    log_info "Compressing output..."
    pigz -p 4 "$OUTPUT_FILE"

    # Generate summary
    log_info "Generating summary..."
    $DORADO_BIN summary "${OUTPUT_FILE}.gz" > "$output_dir/sequencing_summary.txt"

    # Extract metrics
    local total_reads=$(grep -c "^@" "${OUTPUT_FILE}.gz" | head -1 || echo "0")
    local total_bases=$(awk '/^[ACGT]/{sum+=length($0)}END{print sum}' "${OUTPUT_FILE}.gz" || echo "0")

    log_info "Total reads: $total_reads"
    log_info "Total bases: $total_bases"

    # Update database with metrics
    if [[ -n "$DB_HOST" ]] && [[ -n "$RUN_ID" ]]; then
        export PGPASSWORD="${DB_PASSWORD:-}"
        psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c \
            "UPDATE workflow_executions
             SET total_reads = $total_reads,
                 total_bases_gb = $total_bases / 1000000000.0,
                 sequencing_duration_hours = $duration / 3600.0
             WHERE run_id = '$RUN_ID';" 2>/dev/null || true
    fi
}

upload_to_s3() {
    local output_dir=$1

    if [[ "$S3_UPLOAD" != "yes" ]] || [[ -z "$S3_ANALYSIS_BUCKET" ]]; then
        return
    fi

    log_info "Uploading results to S3..."

    # Upload FASTQ
    aws s3 cp "${OUTPUT_FILE}.gz" \
        "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/basecalling/" \
        --metadata "run_id=$RUN_ID,phase=basecalling"

    # Upload summary
    aws s3 cp "$output_dir/sequencing_summary.txt" \
        "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/basecalling/"

    # Upload logs
    aws s3 cp "$LOG_FILE" \
        "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/logs/"

    log_info "S3 upload completed"
}

cleanup() {
    log_info "Cleaning up temporary files..."

    # Remove uncompressed FASTQ if compression succeeded
    [[ -f "${OUTPUT_FILE}.gz" ]] && [[ -f "$OUTPUT_FILE" ]] && rm "$OUTPUT_FILE"

    # Clean work directory
    find "$WORK_DIR" -type f -name "*.tmp" -delete 2>/dev/null || true
}

# ===== Main Execution =====
main() {
    # Parse arguments
    MODEL="$DEFAULT_MODEL"
    DEVICE="cuda:all"
    BATCH_SIZE="auto"
    MIN_QSCORE="9"
    S3_UPLOAD="yes"
    SKIP_DUPLEX="false"
    TRIM=""

    while getopts "i:o:r:m:d:b:q:t:s:xh" opt; do
        case $opt in
            i) INPUT_DIR="$OPTARG" ;;
            o) OUTPUT_DIR="$OPTARG" ;;
            r) RUN_ID="$OPTARG" ;;
            m) MODEL="$OPTARG" ;;
            d) DEVICE="$OPTARG" ;;
            b) BATCH_SIZE="$OPTARG" ;;
            q) MIN_QSCORE="$OPTARG" ;;
            t) TRIM="$OPTARG" ;;
            s) S3_UPLOAD="$OPTARG" ;;
            x) SKIP_DUPLEX="true" ;;
            h) usage ;;
            *) usage ;;
        esac
    done

    # Validate required arguments
    [[ -z "${INPUT_DIR:-}" ]] && log_error "Input directory required (-i)"
    [[ -z "${OUTPUT_DIR:-}" ]] && log_error "Output directory required (-o)"
    [[ -z "${RUN_ID:-}" ]] && log_error "Run ID required (-r)"

    [[ ! -d "$INPUT_DIR" ]] && log_error "Input directory not found: $INPUT_DIR"

    # Setup directories
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$WORK_DIR"
    LOG_FILE="$LOG_DIR/basecalling_${RUN_ID}_$(date +%Y%m%d_%H%M%S).log"

    log_info "========================================="
    log_info "MinION Basecalling Pipeline - Phase 1"
    log_info "========================================="
    log_info "Run ID: $RUN_ID"
    log_info "Start time: $(date)"

    # Execute pipeline
    check_requirements
    perform_basecalling "$INPUT_DIR" "$OUTPUT_DIR"
    upload_to_s3 "$OUTPUT_DIR"
    cleanup

    update_status "completed"

    log_info "========================================="
    log_info "Basecalling completed successfully!"
    log_info "Output: ${OUTPUT_FILE}.gz"
    log_info "End time: $(date)"
    log_info "========================================="

    # Trigger next phase
    if [[ -n "${LAMBDA_ORCHESTRATOR_ARN:-}" ]]; then
        log_info "Triggering Phase 2 (QC)..."
        aws lambda invoke \
            --function-name "$LAMBDA_ORCHESTRATOR_ARN" \
            --payload "{\"run_id\":\"$RUN_ID\",\"phase\":2,\"input\":\"${OUTPUT_FILE}.gz\"}" \
            /dev/null 2>&1 || log_warn "Failed to trigger next phase"
    fi
}

# Trap for cleanup on exit
trap cleanup EXIT

# Run main function
main "$@"