# Exit Code Standards

This document defines the standardized exit codes used across all MinION pipeline scripts.

## Exit Code Convention

| Code | Meaning | Usage |
|------|---------|-------|
| 0 | Success | Normal completion, no issues found |
| 1 | Error | Script failed due to technical error (file not found, invalid input, execution failure) |
| 2 | Warning | Script succeeded but found critical condition requiring attention (PERV detected, critical pathogen found) |

## Rationale

**Exit Code 0 (Success)**:
- Script completed successfully
- No critical findings
- Workflow can proceed to next phase

**Exit Code 1 (Error)**:
- Technical failure preventing script completion
- Invalid parameters or missing required files
- Database connection failures
- Computational errors

Examples:
- Input file not found
- Database not accessible
- Insufficient resources
- Parse errors

**Exit Code 2 (Warning/Critical Finding)**:
- Script completed successfully
- Analysis detected critical condition
- Requires immediate human attention
- Workflow may need to halt

Examples:
- PERV sequences detected (highest priority)
- PMDA critical pathogen detected
- QC thresholds not met but analysis completed
- High host contamination

## Implementation Examples

### Python Scripts

```python
import sys

# Success - no findings
sys.exit(0)

# Error - technical failure
if not input_file.exists():
    print(f"ERROR: Input file not found: {input_file}")
    sys.exit(1)

# Warning - critical finding
if perv_detected:
    print("[CRITICAL] PERV detected")
    sys.exit(2)
elif critical_pathogen_count > 0:
    print("[WARNING] Critical pathogen detected")
    sys.exit(2)
```

### Bash Scripts

```bash
#!/bin/bash
set -euo pipefail

# Error - file not found
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Error - command failed
if ! some_command; then
    echo "ERROR: Command failed"
    exit 1
fi

# Warning - critical finding
if [[ $CRITICAL_COUNT -gt 0 ]]; then
    echo "[WARNING] Critical pathogens detected"
    exit 2
fi

# Success
exit 0
```

## Priority Handling

When multiple conditions exist, use this priority order:

1. **Technical errors (exit 1)**: Always exit immediately on technical failures
2. **Critical findings (exit 2)**: Complete analysis, then exit with code 2
3. **Success (exit 0)**: Only if no errors or critical findings

Example:
```python
# Check for technical errors first
if not validate_inputs():
    sys.exit(1)

# Run analysis
results = run_analysis()

# Check for critical findings
if results.perv_detected:
    sys.exit(2)
elif results.critical_count > 0:
    sys.exit(2)

# Success
sys.exit(0)
```

## Workflow Integration

Pipeline orchestration (Lambda functions) should handle exit codes as follows:

```python
result = subprocess.run(script_cmd, capture_output=True)

if result.returncode == 0:
    # Success - continue to next phase
    trigger_next_phase()
elif result.returncode == 2:
    # Warning - send alert but may continue based on policy
    send_critical_alert()
    # May halt workflow for PERV detection
    if 'PERV' in result.stdout:
        halt_workflow()
else:  # returncode == 1 or other
    # Error - halt workflow and alert
    send_error_alert()
    halt_workflow()
```

## Scripts Using Exit Codes

### Phase 1 - Basecalling
- `basecall_duplex.sh`: 0=success, 1=error

### Phase 2 - QC
- `qc_check.py`: 0=pass, 1=fail
- `run_qc.sh`: 0=pass, 1=fail

### Phase 3 - Host Removal
- `remove_host.sh`: 0=success, 1=error

### Phase 4 - Pathogen Detection
- `kraken2_search.sh`: 0=success, 1=error
- `extract_pmda_pathogens.py`: 0=success, 1=error, 2=critical findings
- `blast_search.sh`: 0=success, 1=error
- `pmda_targeted_search.py`: 0=success, 1=error, 2=critical findings
- `aggregate_results.py`: 0=success, 1=critical findings, 2=PERV detected

### Phase 5 - Quantification
- `absolute_copy_number.py`: 0=success, 1=error
- `run_quantification.sh`: 0=success, 1=error

### Phase 6 - Reporting
- `generate_pmda_report.py`: 0=success, 1=error
- `generate_reports.sh`: 0=success, 1=error

## Testing Exit Codes

```bash
# Test success path
./script.sh && echo "Success (exit 0)"

# Test error handling
./script.sh || echo "Failed with exit code $?"

# Test all paths
./script.sh
EXIT_CODE=$?
case $EXIT_CODE in
    0) echo "Success" ;;
    1) echo "Error" ;;
    2) echo "Critical finding" ;;
    *) echo "Unexpected exit code: $EXIT_CODE" ;;
esac
```

## Notes

- **Never use exit codes > 2** for consistency
- **Always document exit behavior** in script headers
- **Set `set -euo pipefail`** in bash scripts to fail fast
- **Test both success and failure paths** in unit tests
- **Log exit reason** before exiting for debugging
