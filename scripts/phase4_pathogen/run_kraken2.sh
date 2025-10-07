#!/bin/bash
# MinION Metagenomics Pipeline - Phase 4: Kraken2 Pathogen Detection
# Performs taxonomic classification using Kraken2

set -euo pipefail

KRAKEN2_DB="${KRAKEN2_DB:-/mnt/efs/kraken2/standard}"
CONFIDENCE=0.05
REPORT_ZERO=true

usage() {
    echo "Usage: $0 -i INPUT_FASTQ -o OUTPUT_DIR -r RUN_ID [-d DATABASE] [-c CONFIDENCE] [-t THREADS]"
    exit 1
}

THREADS=8
while getopts "i:o:r:d:c:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        d) KRAKEN2_DB="$OPTARG" ;;
        c) CONFIDENCE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Running Kraken2 classification..."
echo "Database: $KRAKEN2_DB"
echo "Confidence threshold: $CONFIDENCE"

# Run Kraken2
kraken2 --db "$KRAKEN2_DB" \
    --threads "$THREADS" \
    --confidence "$CONFIDENCE" \
    --report "$OUTPUT/kraken2_report.txt" \
    --output "$OUTPUT/kraken2_output.txt" \
    --report-zero-counts \
    --use-names \
    "$INPUT" 2>&1 | tee "$OUTPUT/kraken2.log"

# Run Bracken for abundance estimation
for LEVEL in S G F; do
    bracken -d "$KRAKEN2_DB" \
        -i "$OUTPUT/kraken2_report.txt" \
        -o "$OUTPUT/bracken_${LEVEL}.txt" \
        -r 150 \
        -l "$LEVEL" \
        -t 10 2>/dev/null || true
done

# Extract PMDA 91 pathogens
python3 "${SCRIPT_DIR}/extract_pmda_pathogens.py" \
    --report "$OUTPUT/kraken2_report.txt" \
    --output "$OUTPUT/pmda_pathogens.json" \
    --run-id "$RUN_ID"

# Generate Krona visualization
kreport2krona.py -r "$OUTPUT/kraken2_report.txt" -o "$OUTPUT/krona.txt"
ktImportText "$OUTPUT/krona.txt" -o "$OUTPUT/krona.html"

echo "Kraken2 classification completed!"

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/kraken2/"
fi