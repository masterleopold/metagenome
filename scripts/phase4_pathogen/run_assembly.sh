#!/bin/bash
# Phase 4: De novo assembly for pathogen discovery

set -euo pipefail

ASSEMBLER="${ASSEMBLER:-flye}"
GENOME_SIZE="${GENOME_SIZE:-5m}"

usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT_DIR -r RUN_ID [-a ASSEMBLER] [-t THREADS]"
    echo "  Assemblers: flye, canu, megahit"
    exit 1
}

THREADS=8
while getopts "i:o:r:a:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        a) ASSEMBLER="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Running $ASSEMBLER assembly..."

case "$ASSEMBLER" in
    flye)
        flye --nano-raw "$INPUT" \
            --out-dir "$OUTPUT" \
            --threads "$THREADS" \
            --meta \
            --genome-size "$GENOME_SIZE"
        ASSEMBLY="$OUTPUT/assembly.fasta"
        ;;
    canu)
        canu -p assembly -d "$OUTPUT" \
            genomeSize="$GENOME_SIZE" \
            -nanopore "$INPUT" \
            maxThreads="$THREADS"
        ASSEMBLY="$OUTPUT/assembly.contigs.fasta"
        ;;
    megahit)
        seqtk seq -A "$INPUT" > "$OUTPUT/input.fasta"
        megahit -r "$OUTPUT/input.fasta" \
            -o "$OUTPUT" \
            -t "$THREADS" \
            --presets meta-large
        ASSEMBLY="$OUTPUT/final.contigs.fa"
        ;;
esac

# Annotate contigs
echo "Annotating contigs..."
blastn -query "$ASSEMBLY" \
    -db "/mnt/efs/blast/nt" \
    -out "$OUTPUT/contig_annotations.txt" \
    -outfmt 6 \
    -num_threads "$THREADS" \
    -max_target_seqs 1 \
    -evalue 1e-10

# Generate assembly stats
python3 "${SCRIPT_DIR}/assembly_stats.py" \
    --assembly "$ASSEMBLY" \
    --output "$OUTPUT/assembly_stats.json" \
    --run-id "$RUN_ID"

echo "Assembly completed"

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/assembly/"
fi