#!/bin/bash
# MinION Metagenomics Pipeline - Phase 4: Kraken2 Search
# Wrapper script for Kraken2 pathogen detection
# Called by Lambda trigger_pathogen_detection.py

set -euo pipefail

# Default parameters
KRAKEN2_DB="${KRAKEN2_DB:-/mnt/efs/databases/kraken2/standard}"
CONFIDENCE=0.1
MIN_READS=10
REPORT_ZERO=true

usage() {
    cat << EOF
Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -d DATABASE -r RUN_ID [-c CONFIDENCE] [-t THREADS]

Arguments:
    -i  Input directory containing filtered FASTQ files
    -o  Output directory for Kraken2 results
    -d  Kraken2 database path
    -r  Run ID
    -c  Confidence threshold (default: 0.1)
    -t  Number of threads (default: 16)
    -h  Show this help message

Example:
    $0 -i filtered/ -o kraken2/ -d /mnt/efs/databases/kraken2/standard -r RUN-2024-001
EOF
    exit 1
}

# Parse arguments
THREADS=16
while getopts "i:o:d:r:c:t:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) KRAKEN2_DB="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        c) CONFIDENCE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "${INPUT_DIR:-}" || -z "${OUTPUT_DIR:-}" || -z "${RUN_ID:-}" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Get script directory for helper scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Kraken2 Pathogen Detection"
echo "=========================================="
echo "Run ID: $RUN_ID"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Database: $KRAKEN2_DB"
echo "Confidence threshold: $CONFIDENCE"
echo "Threads: $THREADS"
echo "=========================================="

# Verify database exists
if [[ ! -d "$KRAKEN2_DB" ]]; then
    echo "ERROR: Kraken2 database not found: $KRAKEN2_DB"
    exit 1
fi

# Find input FASTQ files
FASTQ_FILES=$(find "$INPUT_DIR" -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" | head -1)

if [[ -z "$FASTQ_FILES" ]]; then
    echo "ERROR: No FASTQ files found in $INPUT_DIR"
    exit 1
fi

# Concatenate all FASTQ files if multiple
if [[ $(find "$INPUT_DIR" -name "*.fastq*" -o -name "*.fq*" | wc -l) -gt 1 ]]; then
    echo "Concatenating multiple FASTQ files..."
    COMBINED_FASTQ="$OUTPUT_DIR/combined_input.fastq"

    # Handle gzipped and uncompressed files
    find "$INPUT_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" -exec zcat {} \; > "$COMBINED_FASTQ"
    find "$INPUT_DIR" -name "*.fastq" -o -name "*.fq" ! -name "combined_input.fastq" -exec cat {} \; >> "$COMBINED_FASTQ"

    INPUT_FASTQ="$COMBINED_FASTQ"
else
    INPUT_FASTQ="$FASTQ_FILES"
fi

# Count input reads
echo "Counting input reads..."
if [[ "$INPUT_FASTQ" == *.gz ]]; then
    TOTAL_READS=$(($(zcat "$INPUT_FASTQ" | wc -l) / 4))
else
    TOTAL_READS=$(($(wc -l < "$INPUT_FASTQ") / 4))
fi

echo "Total input reads: $TOTAL_READS"

if [[ $TOTAL_READS -lt 1000 ]]; then
    echo "WARNING: Very low read count ($TOTAL_READS). Results may be unreliable."
fi

# Run Kraken2 classification
echo "Running Kraken2 classification..."
START_TIME=$(date +%s)

kraken2 \
    --db "$KRAKEN2_DB" \
    --threads "$THREADS" \
    --confidence "$CONFIDENCE" \
    --report "$OUTPUT_DIR/kraken2_report.txt" \
    --output "$OUTPUT_DIR/kraken2_output.txt" \
    --report-zero-counts \
    --use-names \
    --memory-mapping \
    "$INPUT_FASTQ" 2>&1 | tee "$OUTPUT_DIR/kraken2.log"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "Kraken2 classification completed in ${RUNTIME}s"

# Verify output
if [[ ! -f "$OUTPUT_DIR/kraken2_report.txt" ]]; then
    echo "ERROR: Kraken2 report not generated"
    exit 1
fi

# Extract classified read count
CLASSIFIED_READS=$(grep -P "^\s+[0-9]" "$OUTPUT_DIR/kraken2_report.txt" | head -1 | awk '{print $2}')
CLASSIFICATION_RATE=$(awk "BEGIN {printf \"%.2f\", ($CLASSIFIED_READS / $TOTAL_READS) * 100}")

echo "Classified reads: $CLASSIFIED_READS ($CLASSIFICATION_RATE%)"

# Run Bracken for abundance estimation
echo "Running Bracken for abundance estimation..."

# Species level
if [[ -f "$KRAKEN2_DB/database150mers.kmer_distrib" ]]; then
    for LEVEL in S G F; do
        echo "  Bracken level $LEVEL..."
        bracken \
            -d "$KRAKEN2_DB" \
            -i "$OUTPUT_DIR/kraken2_report.txt" \
            -o "$OUTPUT_DIR/bracken_${LEVEL}.txt" \
            -w "$OUTPUT_DIR/bracken_${LEVEL}_report.txt" \
            -r 150 \
            -l "$LEVEL" \
            -t "$MIN_READS" 2>/dev/null || {
                echo "WARNING: Bracken level $LEVEL failed (this is non-critical)"
            }
    done
else
    echo "WARNING: Bracken database not found, skipping abundance estimation"
fi

# Extract PMDA 91 pathogens
echo "Extracting PMDA-designated pathogens..."

python3 "${SCRIPT_DIR}/extract_pmda_pathogens.py" \
    --report "$OUTPUT_DIR/kraken2_report.txt" \
    --output "$OUTPUT_DIR/pmda_pathogens.json" \
    --run-id "$RUN_ID" \
    --min-reads "$MIN_READS" \
    --verbose || {
        echo "ERROR: Failed to extract PMDA pathogens"
        exit 1
    }

# Generate Krona visualization (if available)
if command -v kreport2krona.py &> /dev/null && command -v ktImportText &> /dev/null; then
    echo "Generating Krona visualization..."
    kreport2krona.py -r "$OUTPUT_DIR/kraken2_report.txt" -o "$OUTPUT_DIR/krona.txt"
    ktImportText "$OUTPUT_DIR/krona.txt" -o "$OUTPUT_DIR/krona.html"
    echo "Krona HTML report: $OUTPUT_DIR/krona.html"
else
    echo "WARNING: Krona tools not found, skipping visualization"
fi

# Create summary JSON
cat > "$OUTPUT_DIR/kraken2_summary.json" << EOF
{
  "run_id": "$RUN_ID",
  "database": "$KRAKEN2_DB",
  "total_reads": $TOTAL_READS,
  "classified_reads": $CLASSIFIED_READS,
  "classification_rate": $CLASSIFICATION_RATE,
  "confidence_threshold": $CONFIDENCE,
  "runtime_seconds": $RUNTIME,
  "timestamp": "$(date -Iseconds)"
}
EOF

# Print summary
echo ""
echo "=========================================="
echo "Kraken2 Search Summary"
echo "=========================================="
echo "Total reads: $TOTAL_READS"
echo "Classified reads: $CLASSIFIED_READS ($CLASSIFICATION_RATE%)"
echo "Runtime: ${RUNTIME}s"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="

# Check for critical findings
if [[ -f "$OUTPUT_DIR/pmda_pathogens.json" ]]; then
    CRITICAL_COUNT=$(python3 -c "import json; data=json.load(open('$OUTPUT_DIR/pmda_pathogens.json')); print(data.get('critical_pathogens_detected', 0))")

    if [[ $CRITICAL_COUNT -gt 0 ]]; then
        echo ""
        echo "⚠️  WARNING: $CRITICAL_COUNT critical PMDA pathogen(s) detected!"
        echo "Review $OUTPUT_DIR/pmda_pathogens.json for details"
        echo ""
    fi
fi

# Clean up temporary files
if [[ -f "$OUTPUT_DIR/combined_input.fastq" ]]; then
    rm "$OUTPUT_DIR/combined_input.fastq"
fi

echo "Kraken2 search completed successfully!"
exit 0
