#!/bin/bash
#
# Run NanoPlot QC analysis on FASTQ files
# Generates comprehensive quality reports for nanopore sequencing data
#

set -euo pipefail

# Default values
INPUT_DIR=""
OUTPUT_DIR=""
RUN_ID=""
THREADS=4

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--run-id)
            RUN_ID="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -r RUN_ID [-t THREADS]"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$RUN_ID" ]]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -r RUN_ID [-t THREADS]"
    exit 1
fi

# Validate input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "NanoPlot QC Analysis"
echo "=========================================="
echo "Run ID: $RUN_ID"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "=========================================="

# Find all FASTQ files
FASTQ_FILES=$(find "$INPUT_DIR" -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" | tr '\n' ' ')

if [[ -z "$FASTQ_FILES" ]]; then
    echo "ERROR: No FASTQ files found in $INPUT_DIR"
    exit 1
fi

echo "Found FASTQ files:"
echo "$FASTQ_FILES" | tr ' ' '\n'
echo ""

# Run NanoPlot
START_TIME=$(date +%s)

echo "Running NanoPlot..."
NanoPlot \
    --fastq $FASTQ_FILES \
    --outdir "$OUTPUT_DIR" \
    --prefix "${RUN_ID}_" \
    --threads $THREADS \
    --plots kde hex dot \
    --N50 \
    --title "QC Report: $RUN_ID" \
    --verbose

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "NanoPlot completed successfully!"
echo "Runtime: ${RUNTIME}s"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="

# List generated files
echo ""
echo "Generated files:"
ls -lh "$OUTPUT_DIR"

# Check if NanoStats.txt was created
if [[ ! -f "$OUTPUT_DIR/NanoStats.txt" ]]; then
    echo "ERROR: NanoStats.txt not found in output"
    exit 1
fi

echo ""
echo "NanoStats summary:"
cat "$OUTPUT_DIR/NanoStats.txt"

exit 0
