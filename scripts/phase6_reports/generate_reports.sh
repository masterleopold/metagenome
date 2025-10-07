#!/bin/bash
# Phase 6: Report Generation
# Generates comprehensive PMDA-compliant reports

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -r RUN_ID [-f FORMAT]"
    echo "  FORMAT: pdf, html, json, all (default: all)"
    exit 1
}

FORMAT="all"
while getopts "i:o:r:f:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        f) FORMAT="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT_DIR:-}" || -z "${OUTPUT:-}" || -z "${RUN_ID:-}" ]] && usage

mkdir -p "$OUTPUT"

echo "Generating reports..."

# Generate PMDA checklist
python3 "$SCRIPT_DIR/generate_pmda_checklist.py" \
    --input-dir "$INPUT_DIR" \
    --output "$OUTPUT/pmda_checklist.json" \
    --run-id "$RUN_ID"

# Generate PDF report
if [[ "$FORMAT" == "pdf" ]] || [[ "$FORMAT" == "all" ]]; then
    python3 "$SCRIPT_DIR/generate_pdf_report.py" \
        --input-dir "$INPUT_DIR" \
        --output "$OUTPUT/${RUN_ID}_report.pdf" \
        --run-id "$RUN_ID"
fi

# Generate HTML report
if [[ "$FORMAT" == "html" ]] || [[ "$FORMAT" == "all" ]]; then
    python3 "$SCRIPT_DIR/generate_html_report.py" \
        --input-dir "$INPUT_DIR" \
        --output "$OUTPUT/${RUN_ID}_report.html" \
        --run-id "$RUN_ID"
fi

# Generate JSON report
if [[ "$FORMAT" == "json" ]] || [[ "$FORMAT" == "all" ]]; then
    python3 "$SCRIPT_DIR/generate_pmda_report.py" \
        --input-dir "$INPUT_DIR" \
        --output-dir "$OUTPUT" \
        --format json \
        --run-id "$RUN_ID"
fi

# Create summary
cat > "$OUTPUT/report_summary.json" << EOF
{
  "run_id": "$RUN_ID",
  "report_date": "$(date -Iseconds)",
  "formats_generated": ["pdf", "html", "json"],
  "pmda_compliant": true,
  "files": {
    "pdf": "${RUN_ID}_report.pdf",
    "html": "${RUN_ID}_report.html",
    "json": "${RUN_ID}_pmda_report.json",
    "checklist": "pmda_checklist.json"
  }
}
EOF

echo "Reports generated in $OUTPUT"

# Upload to S3
if [[ -n "${S3_ANALYSIS_BUCKET:-}" ]]; then
    aws s3 sync "$OUTPUT" "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/reports/"

    # Generate presigned URLs for report access
    PDF_URL=$(aws s3 presign "s3://$S3_ANALYSIS_BUCKET/runs/$RUN_ID/reports/${RUN_ID}_report.pdf" --expires-in 604800)
    echo "Report URL (valid for 7 days): $PDF_URL"
fi

# Send completion notification
if [[ -n "${SNS_TOPIC_ARN:-}" ]]; then
    aws sns publish \
        --topic-arn "$SNS_TOPIC_ARN" \
        --subject "Analysis Complete - Run $RUN_ID" \
        --message "Pathogen screening analysis completed for run $RUN_ID. Reports available at: $PDF_URL"
fi