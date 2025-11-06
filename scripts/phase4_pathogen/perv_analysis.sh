#!/bin/bash
# MinION Metagenomics Pipeline - PERV Detection and Analysis
# Specific detection of Porcine Endogenous Retroviruses (critical for xenotransplantation)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PERV_DB="${PERV_DB:-/mnt/efs/perv/perv_reference.fa}"
MIN_COVERAGE_BREADTH=80
MIN_COVERAGE_DEPTH=10

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_DIR -r RUN_ID [-t THREADS]"
    exit 1
}

THREADS=8
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

mkdir -p "$OUTPUT"

echo "Starting PERV-specific analysis..."

# Align to PERV reference sequences
minimap2 -ax map-ont "$PERV_DB" "$INPUT" -t "$THREADS" \
    --secondary=no 2>/dev/null | \
    samtools sort -@ "$THREADS" -o "$OUTPUT/perv_aligned.bam"

samtools index "$OUTPUT/perv_aligned.bam"

# Calculate coverage for each PERV type
for PERV_TYPE in PERV-A PERV-B PERV-C; do
    echo "Analyzing $PERV_TYPE..."

    # Extract reads mapping to specific PERV type
    samtools view -b "$OUTPUT/perv_aligned.bam" "$PERV_TYPE" > "$OUTPUT/${PERV_TYPE}.bam"

    # Calculate coverage
    samtools depth "$OUTPUT/${PERV_TYPE}.bam" > "$OUTPUT/${PERV_TYPE}.depth"

    # Calculate statistics
    READS=$(samtools view -c "$OUTPUT/${PERV_TYPE}.bam")

    if [ $READS -gt 0 ]; then
        AVG_DEPTH=$(awk '{sum+=$3} END {print sum/NR}' "$OUTPUT/${PERV_TYPE}.depth")
        BREADTH=$(awk '$3>0 {count++} END {print count/NR*100}' "$OUTPUT/${PERV_TYPE}.depth")
    else
        AVG_DEPTH=0
        BREADTH=0
    fi

    echo "$PERV_TYPE: $READS reads, ${AVG_DEPTH} avg depth, ${BREADTH}% breadth"

    # Check for full-length reads
    samtools view "$OUTPUT/${PERV_TYPE}.bam" | \
        awk '{if(length($10) > 7000) print}' > "$OUTPUT/${PERV_TYPE}_full_length.txt"
done

# Run PERV typing and recombination detection
python3 "${SCRIPT_DIR}/perv_typing.py" \
    --bam "$OUTPUT/perv_aligned.bam" \
    --output "$OUTPUT/perv_types.json" \
    --run-id "$RUN_ID"

python3 "${SCRIPT_DIR}/detect_recombinants.py" \
    --bam "$OUTPUT/perv_aligned.bam" \
    --output "$OUTPUT/perv_recombinants.json" \
    --run-id "$RUN_ID"

# Generate consensus sequences if coverage is sufficient
for PERV_TYPE in PERV-A PERV-B PERV-C; do
    READS=$(samtools view -c "$OUTPUT/${PERV_TYPE}.bam")
    if [ $READS -gt 100 ]; then
        samtools mpileup -uf "$PERV_DB" "$OUTPUT/${PERV_TYPE}.bam" | \
            bcftools call -c | \
            vcfutils.pl vcf2fq > "$OUTPUT/${PERV_TYPE}_consensus.fq"
    fi
done

# Create summary report
cat > "$OUTPUT/perv_summary.json" << EOF
{
  "run_id": "$RUN_ID",
  "perv_detected": $([ -s "$OUTPUT/PERV-A.bam" ] || [ -s "$OUTPUT/PERV-B.bam" ] || [ -s "$OUTPUT/PERV-C.bam" ] && echo "true" || echo "false"),
  "timestamp": "$(date -Iseconds)"
}
EOF

echo "PERV analysis completed!"

# Send critical alert if PERV detected
if [ -s "$OUTPUT/PERV-A.bam" ] || [ -s "$OUTPUT/PERV-B.bam" ] || [ -s "$OUTPUT/PERV-C.bam" ]; then
    echo "WARNING: PERV sequences detected! Review required for xenotransplantation safety."

    if [[ -n "${SNS_CRITICAL_TOPIC:-}" ]]; then
        aws sns publish \
            --topic-arn "$SNS_CRITICAL_TOPIC" \
            --subject "CRITICAL: PERV Detected - Run $RUN_ID" \
            --message "PERV sequences detected in run $RUN_ID. Immediate review required."
    fi
fi

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/perv/"
fi