#!/bin/bash
# MinION Metagenomics Pipeline - Phase 3: Host Removal
# Removes pig host sequences to retain only microbial reads

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HOST_GENOME="${HOST_GENOME:-/mnt/efs/host_genome/sus_scrofa_11.1.fa}"
MIN_IDENTITY=90  # Minimum identity for host match

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_DIR -r RUN_ID [-g GENOME] [-t THREADS]"
    exit 1
}

THREADS=8
while getopts "i:o:r:g:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        g) HOST_GENOME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Starting host removal..."
echo "Input: $INPUT"
echo "Host genome: $HOST_GENOME"

# Align to host genome using minimap2
minimap2 -ax map-ont "$HOST_GENOME" "$INPUT" -t "$THREADS" \
    --secondary=no 2>/dev/null | \
    samtools view -bS - > "$OUTPUT/aligned_to_host.bam"

# Extract unmapped reads (non-host)
samtools view -b -f 4 "$OUTPUT/aligned_to_host.bam" > "$OUTPUT/unmapped.bam"

# Convert unmapped BAM to FASTQ
samtools fastq "$OUTPUT/unmapped.bam" | pigz -p 4 > "$OUTPUT/non_host_reads.fastq.gz"

# Calculate statistics
TOTAL=$(samtools view -c "$OUTPUT/aligned_to_host.bam")
UNMAPPED=$(samtools view -c -f 4 "$OUTPUT/aligned_to_host.bam")
MAPPED=$((TOTAL - UNMAPPED))
REMOVAL_RATE=$(echo "scale=2; $MAPPED * 100 / $TOTAL" | bc)

# Save statistics
cat > "$OUTPUT/host_removal_stats.json" << EOF
{
  "run_id": "$RUN_ID",
  "total_reads": $TOTAL,
  "host_reads": $MAPPED,
  "non_host_reads": $UNMAPPED,
  "host_removal_rate": $REMOVAL_RATE,
  "timestamp": "$(date -Iseconds)"
}
EOF

echo "Host removal completed!"
echo "Total reads: $TOTAL"
echo "Host reads removed: $MAPPED ($REMOVAL_RATE%)"
echo "Non-host reads: $UNMAPPED"

# Clean up intermediate files
rm -f "$OUTPUT/aligned_to_host.bam" "$OUTPUT/unmapped.bam"

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 cp "$OUTPUT/non_host_reads.fastq.gz" \
        "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/host_removal/"
    aws s3 cp "$OUTPUT/host_removal_stats.json" \
        "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/host_removal/"
fi