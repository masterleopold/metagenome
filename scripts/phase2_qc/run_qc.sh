#!/bin/bash
# MinION Metagenomics Pipeline - Phase 2: Quality Control
# Performs comprehensive QC analysis on basecalled reads

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-/data/logs}"

# QC thresholds (PMDA requirements)
MIN_READS=100000
MIN_MEAN_Q=9
MIN_N50=200
MAX_FAILED_PCT=50

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_DIR -r RUN_ID [-t THREADS]"
    exit 1
}

THREADS=4
while getopts "i:o:r:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT" "$LOG_DIR"
LOG_FILE="$LOG_DIR/qc_${RUN_ID}_$(date +%Y%m%d_%H%M%S).log"

echo "Starting QC analysis..." | tee "$LOG_FILE"

# Run NanoPlot
NanoPlot --fastq "$INPUT" \
    --outdir "$OUTPUT/nanoplot" \
    --prefix "${RUN_ID}_" \
    --plots hex dot \
    --N50 \
    --threads "$THREADS" \
    --title "Run $RUN_ID QC Report" 2>&1 | tee -a "$LOG_FILE"

# Run NanoStat
NanoStat --fastq "$INPUT" --threads "$THREADS" > "$OUTPUT/nanostats.txt"

# Run FastQC
mkdir -p "$OUTPUT/fastqc"
fastqc -o "$OUTPUT/fastqc" -t "$THREADS" "$INPUT" 2>&1 | tee -a "$LOG_FILE"

# Extract key metrics
TOTAL_READS=$(grep "Number of reads:" "$OUTPUT/nanostats.txt" | awk '{print $NF}')
MEAN_Q=$(grep "Mean read quality:" "$OUTPUT/nanostats.txt" | awk '{print $NF}')
N50=$(grep "Read N50:" "$OUTPUT/nanostats.txt" | awk '{print $NF}')

# Check QC pass/fail
QC_PASS="true"
[[ $TOTAL_READS -lt $MIN_READS ]] && QC_PASS="false"
[[ $(echo "$MEAN_Q < $MIN_MEAN_Q" | bc) -eq 1 ]] && QC_PASS="false"
[[ $N50 -lt $MIN_N50 ]] && QC_PASS="false"

# Generate QC report
cat > "$OUTPUT/qc_summary.json" << EOF
{
  "run_id": "$RUN_ID",
  "qc_pass": $QC_PASS,
  "metrics": {
    "total_reads": $TOTAL_READS,
    "mean_qscore": $MEAN_Q,
    "n50": $N50
  },
  "thresholds": {
    "min_reads": $MIN_READS,
    "min_mean_q": $MIN_MEAN_Q,
    "min_n50": $MIN_N50
  },
  "timestamp": "$(date -Iseconds)"
}
EOF

# Update database
if [[ -n "${DB_HOST:-}" ]]; then
    python3 "${SCRIPT_DIR}/qc_check.py" \
        --summary "$OUTPUT/qc_summary.json" \
        --run-id "$RUN_ID" \
        --db-host "$DB_HOST"
fi

echo "QC analysis completed. Pass: $QC_PASS" | tee -a "$LOG_FILE"

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/qc/"
fi

# Exit with error if QC failed
[[ "$QC_PASS" == "false" ]] && exit 1
exit 0