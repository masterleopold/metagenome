#!/bin/bash
# Phase 5: Pathogen Quantification
# Performs absolute and relative quantification of detected pathogens

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPIKE_IN="${SPIKE_IN:-PhiX174}"

usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -r RUN_ID [-s SPIKE_IN]"
    exit 1
}

while getopts "i:o:r:s:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        s) SPIKE_IN="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT_DIR:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Starting quantification analysis..."

# Quantify from Kraken2 results (all 91 PMDA pathogens)
if [ -f "$INPUT_DIR/kraken2/kraken2_report.txt" ]; then
    # Config file location (AMI or repo)
    CONFIG_FILE="${PMDA_CONFIG:-/opt/minion/templates/config/pmda_pathogens.json}"

    python3 "$SCRIPT_DIR/kraken_quantify.py" \
        --report "$INPUT_DIR/kraken2/kraken2_report.txt" \
        --output "$OUTPUT/kraken_quantification.json" \
        --run-id "$RUN_ID" \
        --config "$CONFIG_FILE"
fi

# Quantify from BLAST results
if [ -f "$INPUT_DIR/blast/blast_results.txt" ]; then
    python3 "$SCRIPT_DIR/blast_quantify.py" \
        --blast "$INPUT_DIR/blast/blast_results.txt" \
        --output "$OUTPUT/blast_quantification.json" \
        --run-id "$RUN_ID"
fi

# Spike-in normalization
if [ -f "$INPUT_DIR/kraken2/kraken2_report.txt" ]; then
    python3 "$SCRIPT_DIR/spike_in_normalization.py" \
        --report "$INPUT_DIR/kraken2/kraken2_report.txt" \
        --spike-in "$SPIKE_IN" \
        --output "$OUTPUT/normalized_quantification.json" \
        --run-id "$RUN_ID"
fi

# Calculate absolute copy numbers
python3 "$SCRIPT_DIR/absolute_copy_number.py" \
    --kraken "$OUTPUT/kraken_quantification.json" \
    --normalized "$OUTPUT/normalized_quantification.json" \
    --output "$OUTPUT/absolute_quantification.json" \
    --plasma-volume 10 \
    --run-id "$RUN_ID"

echo "Quantification completed"

# Generate summary
cat > "$OUTPUT/quantification_summary.json" << EOF
{
  "run_id": "$RUN_ID",
  "quantification_methods": ["kraken2", "blast", "spike_in"],
  "spike_in_control": "$SPIKE_IN",
  "timestamp": "$(date -Iseconds)"
}
EOF

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/quantification/"
fi