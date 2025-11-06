#!/bin/bash
# Phase 4: BLAST-based pathogen detection

set -euo pipefail

BLAST_DB="${BLAST_DB:-/mnt/efs/blast/nt}"
EVALUE="1e-10"
MAX_TARGETS=10

usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT_DIR -r RUN_ID [-d DB] [-e EVALUE] [-t THREADS]"
    exit 1
}

THREADS=8
while getopts "i:o:r:d:e:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        d) BLAST_DB="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Running BLAST search..."

# Convert FASTQ to FASTA
if [[ "$INPUT" == *.fastq* ]]; then
    seqtk seq -A "$INPUT" > "$OUTPUT/input.fasta"
    INPUT_FASTA="$OUTPUT/input.fasta"
else
    INPUT_FASTA="$INPUT"
fi

# Run BLASTN
blastn -query "$INPUT_FASTA" \
    -db "$BLAST_DB" \
    -out "$OUTPUT/blast_results.txt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -evalue "$EVALUE" \
    -num_threads "$THREADS" \
    -max_target_seqs "$MAX_TARGETS"

# Parse results for PMDA pathogens
python3 "${SCRIPT_DIR}/parse_blast_results.py" \
    --blast "$OUTPUT/blast_results.txt" \
    --output "$OUTPUT/blast_pmda_pathogens.json" \
    --run-id "$RUN_ID"

echo "BLAST analysis completed"

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/blast/"
fi