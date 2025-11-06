#!/bin/bash
# MinION Metagenomics Pipeline - Phase 4: BLAST Search
# Searches filtered reads against viral database (RVDB)
# Called by Lambda trigger_pathogen_detection.py

set -euo pipefail

# Default parameters
BLAST_DB="${BLAST_DB:-/mnt/efs/databases/rvdb/rvdb.fasta}"
EVALUE=1e-5
MAX_TARGET_SEQS=10
MAX_HSPS=1
TASK="blastn"

usage() {
    cat << EOF
Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -d DATABASE -r RUN_ID [-e EVALUE] [-t THREADS]

Arguments:
    -i  Input directory containing filtered FASTQ files
    -o  Output directory for BLAST results
    -d  BLAST database path
    -r  Run ID
    -e  E-value threshold (default: 1e-5)
    -t  Number of threads (default: 16)
    -h  Show this help message

Example:
    $0 -i filtered/ -o blast/ -d /mnt/efs/databases/rvdb/rvdb.fasta -r RUN-2024-001
EOF
    exit 1
}

# Parse arguments
THREADS=16
while getopts "i:o:d:r:e:t:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) BLAST_DB="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
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

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "BLAST Viral Search"
echo "=========================================="
echo "Run ID: $RUN_ID"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Database: $BLAST_DB"
echo "E-value: $EVALUE"
echo "Threads: $THREADS"
echo "=========================================="

# Verify database exists
if [[ ! -f "$BLAST_DB" ]]; then
    echo "ERROR: BLAST database not found: $BLAST_DB"
    exit 1
fi

# Check if database is indexed
if [[ ! -f "${BLAST_DB}.nhr" && ! -f "${BLAST_DB}.nal" ]]; then
    echo "WARNING: BLAST database not indexed. Creating index..."
    makeblastdb -in "$BLAST_DB" -dbtype nucl -parse_seqids
fi

# Find input FASTQ files
echo "Finding input files..."
FASTQ_FILES=$(find "$INPUT_DIR" -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz")

if [[ -z "$FASTQ_FILES" ]]; then
    echo "ERROR: No FASTQ files found in $INPUT_DIR"
    exit 1
fi

# Convert FASTQ to FASTA (BLAST requires FASTA)
echo "Converting FASTQ to FASTA..."
FASTA_INPUT="$OUTPUT_DIR/input_sequences.fasta"

{
    while IFS= read -r file; do
        if [[ "$file" == *.gz ]]; then
            zcat "$file"
        else
            cat "$file"
        fi
    done <<< "$FASTQ_FILES"
} | awk 'NR%4==1 {print ">"substr($0,2)} NR%4==2 {print}' > "$FASTA_INPUT"

# Count sequences
SEQ_COUNT=$(grep -c "^>" "$FASTA_INPUT" || true)
echo "Total sequences to search: $SEQ_COUNT"

if [[ $SEQ_COUNT -eq 0 ]]; then
    echo "ERROR: No sequences found in input"
    exit 1
fi

# Subsample if too many sequences (BLAST can be slow)
MAX_SEQUENCES=100000
if [[ $SEQ_COUNT -gt $MAX_SEQUENCES ]]; then
    echo "WARNING: Too many sequences ($SEQ_COUNT). Subsampling to $MAX_SEQUENCES..."
    SUBSAMPLED="$OUTPUT_DIR/input_sequences_subsampled.fasta"

    # Calculate sampling rate
    SAMPLE_RATE=$(awk "BEGIN {printf \"%.6f\", $MAX_SEQUENCES / $SEQ_COUNT}")

    seqtk sample "$FASTA_INPUT" "$SAMPLE_RATE" > "$SUBSAMPLED"
    FASTA_INPUT="$SUBSAMPLED"
    SEQ_COUNT=$(grep -c "^>" "$FASTA_INPUT")
    echo "Subsampled to $SEQ_COUNT sequences"
fi

# Run BLAST search
echo "Running BLAST search..."
START_TIME=$(date +%s)

blastn \
    -query "$FASTA_INPUT" \
    -db "$BLAST_DB" \
    -out "$OUTPUT_DIR/blast_results.txt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -evalue "$EVALUE" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -max_hsps "$MAX_HSPS" \
    -num_threads "$THREADS" \
    -task "$TASK" \
    2>&1 | tee "$OUTPUT_DIR/blast.log"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "BLAST search completed in ${RUNTIME}s"

# Check if results were generated
if [[ ! -f "$OUTPUT_DIR/blast_results.txt" ]]; then
    echo "ERROR: BLAST results not generated"
    exit 1
fi

# Count hits
HIT_COUNT=$(wc -l < "$OUTPUT_DIR/blast_results.txt")
echo "Total BLAST hits: $HIT_COUNT"

# Parse and filter high-quality hits
echo "Filtering high-quality hits (identity >= 90%, e-value <= 1e-5)..."

awk -F'\t' '$3 >= 90 && $11 <= 1e-5' "$OUTPUT_DIR/blast_results.txt" \
    > "$OUTPUT_DIR/blast_filtered.txt"

FILTERED_HITS=$(wc -l < "$OUTPUT_DIR/blast_filtered.txt")
echo "High-quality hits: $FILTERED_HITS"

# Extract unique viral families
echo "Extracting unique viral families..."

awk -F'\t' '{print $13}' "$OUTPUT_DIR/blast_filtered.txt" | \
    sed 's/.*\[\(.*\)\].*/\1/' | \
    sort | uniq -c | sort -rn > "$OUTPUT_DIR/viral_families.txt"

# Generate top hits summary
echo "Generating top hits summary..."

head -50 "$OUTPUT_DIR/blast_filtered.txt" > "$OUTPUT_DIR/top_50_hits.txt"

# Create JSON summary
cat > "$OUTPUT_DIR/blast_summary.json" << EOF
{
  "run_id": "$RUN_ID",
  "database": "$BLAST_DB",
  "total_sequences_searched": $SEQ_COUNT,
  "total_hits": $HIT_COUNT,
  "high_quality_hits": $FILTERED_HITS,
  "evalue_threshold": $EVALUE,
  "identity_filter": 90,
  "runtime_seconds": $RUNTIME,
  "timestamp": "$(date -Iseconds)"
}
EOF

# Extract viral species for reporting
if [[ $FILTERED_HITS -gt 0 ]]; then
    echo "Top viral detections:"
    head -20 "$OUTPUT_DIR/viral_families.txt" | while read -r count family; do
        echo "  $count hits: $family"
    done

    # Create detailed hits JSON for aggregation
    python3 << 'PYEOF' > "$OUTPUT_DIR/blast_hits.json"
import json
import sys

hits = []
try:
    with open("$OUTPUT_DIR/blast_filtered.txt") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 13:
                hit = {
                    'query': parts[0],
                    'subject': parts[1],
                    'identity': float(parts[2]),
                    'length': int(parts[3]),
                    'evalue': float(parts[10]),
                    'bitscore': float(parts[11]),
                    'description': parts[12]
                }
                hits.append(hit)
except Exception as e:
    print(f"Error parsing BLAST results: {e}", file=sys.stderr)
    sys.exit(1)

# Get top hits by bitscore
hits_sorted = sorted(hits, key=lambda x: x['bitscore'], reverse=True)

output = {
    'run_id': '$RUN_ID',
    'total_hits': len(hits),
    'top_hits': hits_sorted[:100]
}

print(json.dumps(output, indent=2))
PYEOF

else
    echo "No significant viral hits detected"

    cat > "$OUTPUT_DIR/blast_hits.json" << EOF
{
  "run_id": "$RUN_ID",
  "total_hits": 0,
  "top_hits": []
}
EOF
fi

# Clean up temporary files
rm -f "$OUTPUT_DIR/input_sequences.fasta"
rm -f "$OUTPUT_DIR/input_sequences_subsampled.fasta"

echo ""
echo "=========================================="
echo "BLAST Search Summary"
echo "=========================================="
echo "Sequences searched: $SEQ_COUNT"
echo "Total hits: $HIT_COUNT"
echo "High-quality hits: $FILTERED_HITS"
echo "Runtime: ${RUNTIME}s"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="

# Check for significant viral findings
if [[ $FILTERED_HITS -gt 100 ]]; then
    echo ""
    echo "⚠️  WARNING: High number of viral hits detected ($FILTERED_HITS)"
    echo "Review $OUTPUT_DIR/blast_filtered.txt for details"
    echo ""
fi

echo "BLAST search completed successfully!"
exit 0
