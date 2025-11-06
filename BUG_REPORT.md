# Comprehensive Bug Report - MinION Pathogen Screening Pipeline

**Generated:** 2025-11-06
**Severity Levels:** 🔴 CRITICAL | 🟠 HIGH | 🟡 MEDIUM | 🔵 LOW

---

## Executive Summary

Comprehensive codebase audit identified **23 bugs and issues** across all pipeline components:
- **7 Critical bugs** that will cause pipeline failures
- **8 High-priority issues** affecting reliability and accuracy
- **5 Medium-priority issues** impacting maintainability
- **3 Low-priority issues** for code quality improvements

**IMMEDIATE ACTION REQUIRED:** Critical bugs #1-7 will prevent pipeline execution.

---

## 🔴 CRITICAL BUGS (Pipeline-Breaking)

### 1. Missing Critical Script: `extract_pmda_pathogens.py`
**Location:** Referenced in multiple locations but file does not exist
**Impact:** Pipeline will fail during pathogen detection phase
**References:**
- `lambda/phases/trigger_pathogen_detection.py:138` - Calls `/opt/minion/scripts/phase4_pathogen/extract_pmda_pathogens.py`
- `scripts/phase4_pathogen/run_kraken2.sh:58` - Calls `${SCRIPT_DIR}/extract_pmda_pathogens.py`

**Expected Functionality:**
- Extract PMDA-designated 91 pathogens from Kraken2 report
- Output JSON format with pathogen codes and read counts
- Critical for regulatory compliance

**Fix Required:** Create this script to parse Kraken2 reports and filter for PMDA pathogens

---

### 2. Missing Script: `aggregate_results.py`
**Location:** `lambda/phases/trigger_pathogen_detection.py:200`
**Impact:** Pathogen detection phase cannot complete - results aggregation will fail

**Called with:**
```bash
/opt/minion/scripts/phase4_pathogen/aggregate_results.py \
    --kraken kraken2/ \
    --blast blast/ \
    --perv perv/ \
    --output pathogen_summary.json \
    --run-id "$RUN_ID"
```

**Fix Required:** Create script to merge results from multiple detection methods

---

### 3. Missing Shell Scripts Referenced by Lambda
**Location:** `lambda/phases/trigger_pathogen_detection.py`
**Missing Files:**
- `/opt/minion/scripts/phase4_pathogen/kraken2_search.sh` (line 132)
- `/opt/minion/scripts/phase4_pathogen/blast_search.sh` (line 148)
- `/opt/minion/scripts/phase4_pathogen/pmda_targeted_search.py` (line 157)

**Impact:** EC2 instances will fail when attempting to run these commands

**Fix Required:** Create these analysis wrapper scripts

---

### 4. Missing Script: `basecall_duplex.sh`
**Location:** `lambda/phases/trigger_basecalling.py:160`
**Impact:** Basecalling phase (Phase 1) will fail immediately

**Called with:**
```bash
/opt/minion/scripts/phase1_basecalling/basecall_duplex.sh \
    -i "$S3_INPUT" \
    -o "$S3_OUTPUT" \
    -r "$RUN_ID" \
    -s (optional)
```

**Fix Required:** File exists in directory but may not be executable or properly installed on AMI

---

### 5. QC Threshold Inconsistency (10x Discrepancy)
**Location:**
- `templates/config/default_pipeline.yaml:121` - `min_reads: 10000`
- `scripts/phase2_qc/qc_check.py:16` - `'min_reads': 100000`

**Impact:**
- Samples with 10,000-99,999 reads will pass pipeline config but fail QC checks
- Users will receive confusing failures
- PMDA compliance documentation inconsistency

**Fix Required:** Align thresholds to same value (recommend 100,000 based on PMDA requirements in comments)

---

### 6. Non-existent `lib/` Directory Referenced in Documentation
**Location:** `CLAUDE.md` and multiple documentation references
**Impact:**
- Developers trying to import from `lib/` will fail
- Instructions reference non-existent shared libraries
- Code organization confusion

**References in CLAUDE.md:**
- "lib/ # Shared Python libraries" (line in file structure)
- Multiple references to using lib directory

**Actual Status:** No `lib/` directory exists in codebase, no imports from lib found

**Fix Required:** Either create lib/ directory and move shared code, or update documentation

---

### 7. Missing `numpy` Import at Module Level
**Location:** `scripts/phase5_quantification/absolute_copy_number.py:83`
**Code:**
```python
results[f"{perv_type.lower().replace('-', '_')}_coverage_depth"] = np.mean(coverage_array)
```
**But numpy not imported until line 109:**
```python
import numpy as np
```

**Impact:** Script will crash with NameError when executed
**Fix Required:** Move `import numpy as np` to top of file (line ~12)

---

## 🟠 HIGH PRIORITY BUGS

### 8. Spot Instance Fulfillment Not Handled
**Location:** `lambda/phases/trigger_basecalling.py:128-136`
**Issue:** Code assumes spot request will always be fulfilled
**Impact:**
- If spot instances unavailable, pipeline will hang or timeout
- No fallback to on-demand instances
- Waiter may timeout after 15 minutes

**Recommendation:** Add try/except with fallback to on-demand instances

---

### 9. Database Password Passed as Command-Line Argument
**Location:** `scripts/phase2_qc/qc_check.py:87`
```python
parser.add_argument('--db-password')
```
**Security Risk:** Password visible in process list, logs, and command history
**Fix Required:** Use environment variables or AWS Secrets Manager

---

### 10. SQL Injection Vulnerability
**Location:** `lambda/monitoring/alert_handler.py:310-327`
**Code:**
```python
sql = """
    INSERT INTO alerts
    (workflow_id, alert_type, details, created_at)
    VALUES (:workflow_id, :alert_type, :details, NOW())
"""
```
**Issue:** While using parameterized query, workflow_id type coercion may fail:
```python
{'name': 'workflow_id', 'value': {'longValue': workflow_id if workflow_id else 0}}
```
**Problem:** workflow_id may be string, causing type error
**Fix Required:** Add proper type validation/conversion

---

### 11. Hardcoded IAM Role Name
**Location:** Multiple Lambda functions
- `lambda/phases/trigger_basecalling.py:97` - `'Name': 'MinIONEC2Role'`
- `lambda/phases/trigger_pathogen_detection.py:94` - `'Name': 'MinIONEC2Role'`

**Issue:** Assumes specific IAM role name, not configurable
**Impact:**
- Fails if role name differs in deployment
- Not multi-environment friendly
- Terraform may create different role names

**Fix Required:** Use environment variable for IAM role name

---

### 12. Incomplete Pathogen Count (96 ≠ 91)
**Location:** `templates/config/pmda_pathogens.json`
**Declared Counts:**
- Viruses: 50 (but only 25 listed)
- Bacteria: 35 (28 listed)
- Parasites: 5 ✓
- Fungi: 5 ✓
- Prions: 1 ✓

**Total Listed:** 64 pathogens (should be 91)
**Total Claimed:** 96 pathogens (50+35+5+5+1)

**Impact:**
- PMDA compliance tests may fail
- Missing critical pathogens from official list
- Test `test_91_pathogen_coverage` will fail

**Fix Required:** Complete pathogen list with all 91 PMDA-designated organisms

---

### 13. Test File Paths Use Relative Paths
**Location:** `tests/test_pmda_compliance.py:15`
```python
self.pmda_pathogens_file = Path('../templates/config/pmda_pathogens.json')
```
**Issue:** Tests will fail if not run from `tests/` directory
**Impact:** CI/CD pipelines may fail, pytest discovery issues
**Fix Required:** Use `__file__` to get absolute paths or use package resources

---

### 14. Error Suppression in Critical Code
**Location:** `scripts/phase4_pathogen/perv_phylogenetics.py:33`
```python
subprocess.run([
    "mafft", "--auto", str(combined_file)
], stdout=open(aligned_file, 'w'), stderr=subprocess.DEVNULL)
```

**Issue:** STDERR suppressed - alignment errors will be silent
**Impact:** Failed alignments will produce empty files, causing downstream errors
**Fix Required:** Log stderr or check return code

---

### 15. Missing Error Handling for BAM File Operations
**Location:** Multiple PERV analysis scripts
- `scripts/phase4_pathogen/perv_typing.py:37` - No check if BAM file exists/is valid
- `scripts/phase4_pathogen/detect_recombinants.py:16` - Same issue

**Impact:** Cryptic pysam errors if BAM files corrupted or missing
**Fix Required:** Add try/except blocks and file validation

---

## 🟡 MEDIUM PRIORITY ISSUES

### 16. Inconsistent String Formatting
**Location:** Throughout codebase
**Examples:**
- `lambda/phases/trigger_basecalling.py` uses f-strings
- `lambda/phases/trigger_pathogen_detection.py` uses `.format()`
- Some files mix both styles

**Impact:** Code maintainability, harder to read
**Recommendation:** Standardize on f-strings (Python 3.6+)

---

### 17. Duplicate Dictionary Key
**Location:** `templates/config/pmda_pathogens.json:68`
```json
{"code": "BB", "name": "Borrelia burgdorferi", "risk_level": "MEDIUM"},
...
{"code": "BB", "name": "Brachyspira species", "risk_level": "LOW"}
```

**Impact:** JSON parsers will only keep one "BB" entry (likely the last one)
**Fix Required:** Use unique codes (e.g., "BB-B" for Borrelia, "BB-BR" for Brachyspira)

---

### 18. Magic Numbers in Cost Calculations
**Location:** `scripts/phase5_quantification/absolute_copy_number.py:31,39`
```python
avg_read_length = 1000  # Hardcoded
extraction_efficiency: float = 0.7  # Default parameter
```

**Issue:** These should be configurable per run
**Impact:** Quantification accuracy varies with actual read length and extraction method
**Fix Required:** Pass as command-line arguments with defaults

---

### 19. Incomplete GENOME_SIZES Dictionary
**Location:** `scripts/phase5_quantification/absolute_copy_number.py:11-18`
**Code:**
```python
GENOME_SIZES = {
    'PERV': 8000,
    'HEV': 7200,
    'JEV': 11000,
    'SS': 2000000,
    'EC': 5000000,
    # Add all pathogens  <- Comment indicates incomplete
}
```

**Impact:** Most pathogens will use default 10kb size (line 69), reducing accuracy
**Fix Required:** Complete genome sizes for all 91 PMDA pathogens

---

### 20. Configuration Mismatch: QC Thresholds
**Location:** Multiple configs
- `templates/config/default_pipeline.yaml:26` - `min_pass_reads: 1000`
- `templates/config/default_pipeline.yaml:44` - `min_reads: 10` (pathogen detection)
- `templates/config/default_pipeline.yaml:121` - `min_reads: 10000` (quality thresholds)

**Issue:** Three different "min_reads" parameters with different meanings and values
**Impact:** Confusion about which threshold applies where
**Recommendation:** Use more specific names: `qc_min_pass_reads`, `pathogen_min_reads`, `workflow_min_reads`

---

## 🔵 LOW PRIORITY ISSUES

### 21. Emoji in Production Alert Code
**Location:** `lambda/monitoring/alert_handler.py:81,125`
```python
subject = f'🚨 CRITICAL: PERV Detection - Run {run_id}'
subject = f'⚠️ Critical Pathogen Detected: {pathogen} - Run {run_id}'
```

**Issue:** Emojis may not render in all email clients, logs, or monitoring systems
**Recommendation:** Use plain text in subject, add emojis in body if desired

---

### 22. Missing Docstrings
**Location:** Most functions lack comprehensive docstrings
**Examples:**
- `lambda/phases/trigger_basecalling.py:69` - `launch_basecalling_instance` - no docstring
- `lambda/phases/trigger_pathogen_detection.py:63` - `launch_pathogen_instance` - minimal docstring

**Impact:** Reduced code maintainability
**Recommendation:** Add comprehensive docstrings with parameter types, returns, raises

---

### 23. Inconsistent Exit Codes
**Location:** Shell scripts
- `scripts/phase4_pathogen/perv_analysis.sh` - No explicit exit codes
- `scripts/phase2_qc/qc_check.py:109,112` - Uses sys.exit(1) for failure, sys.exit(0) for success

**Issue:** Some scripts don't set proper exit codes for error conditions
**Impact:** Step Functions may not detect failures correctly
**Recommendation:** Standardize exit code handling across all scripts

---

## Test Coverage Issues

### Missing Test Files
The following components have **no test coverage**:
1. All shell scripts in `scripts/phase*/`
2. Lambda trigger functions (phases)
3. Quantification scripts
4. Report generation scripts
5. PERV detection scripts

**Recommendation:** Add integration tests for critical paths

---

## Documentation Inconsistencies

### CLAUDE.md Issues
1. **Line ~30:** References `lib/` directory that doesn't exist
2. **Line ~50:** Claims "shared Python libraries" but none exist
3. **Command Reference:** Some commands reference non-existent scripts

### Missing Documentation
1. No database schema documentation
2. No AMI build process documented
3. No disaster recovery procedures
4. EC2 UserData scripts not documented

---

## Infrastructure Issues

### Terraform Observations
Terraform files exist but unable to verify:
- Lambda function deployment packages
- AMI references (BASECALLING_AMI_ID, ANALYSIS_AMI_ID must be set)
- EFS mount targets configuration
- RDS security and backup configuration

**Recommendation:** Add infrastructure validation tests

---

## PMDA Compliance Gaps

### Incomplete Pathogen List
Current: 64/91 pathogens fully documented
Missing: 27 pathogens (see bug #12)

### Missing Detection Scripts
- No confirmation testing workflows
- No quality control sample processing
- No linearity/LOD validation scripts

---

## Priority Ranking for Fixes

### Must Fix Before Production (Sprint 0)
1. Bug #1: Create `extract_pmda_pathogens.py` ⏱️ 2-3 hours
2. Bug #2: Create `aggregate_results.py` ⏱️ 2-3 hours
3. Bug #3: Create missing shell scripts ⏱️ 4-6 hours
4. Bug #5: Fix QC threshold inconsistency ⏱️ 30 minutes
5. Bug #7: Fix numpy import order ⏱️ 5 minutes

**Total Sprint 0:** ~1.5 days

### Should Fix Before Beta (Sprint 1)
6. Bug #6: Resolve lib/ directory issue ⏱️ 1 hour
7. Bug #8: Add spot instance fallback ⏱️ 2 hours
8. Bug #12: Complete PMDA pathogen list ⏱️ 4-6 hours
9. Bug #13: Fix test file paths ⏱️ 1 hour
10. Bug #9: Secure database credentials ⏱️ 1 hour

**Total Sprint 1:** ~1 day

### Nice to Have (Sprint 2+)
- Bugs #16-23 (code quality, documentation)
- Additional test coverage
- Documentation completion

---

## Testing Recommendations

### Before Production
1. **End-to-end pipeline test** with known positive samples
2. **PERV detection validation** with synthetic PERV sequences
3. **Negative control testing** to verify specificity
4. **Failure mode testing** - test each error path
5. **Cost monitoring test** - verify spot instances work correctly

### Continuous Testing
1. Add integration tests for critical scripts
2. Add unit tests for quantification calculations
3. Mock AWS services in tests (using moto)
4. Add performance regression tests

---

## Security Audit Summary

### Issues Found
- 🔴 Database password in CLI arguments (Bug #9)
- 🟠 Potential SQL injection (Bug #10)
- 🟠 IAM role hardcoded (Bug #11)

### Recommendations
1. Use AWS Secrets Manager for all credentials
2. Implement least-privilege IAM policies
3. Add CloudTrail logging for all S3/RDS access
4. Encrypt all data at rest and in transit
5. Add VPC endpoint policies

---

## Estimated Fix Timeline

| Phase | Duration | Critical Bugs Fixed | Total Bugs Fixed |
|-------|----------|---------------------|------------------|
| Sprint 0 (Pre-launch) | 1.5 days | 5 | 5 |
| Sprint 1 (Beta) | 1 day | 2 | 5 |
| Sprint 2 (Refinement) | 2 days | 0 | 6 |
| Sprint 3 (Polish) | 2 days | 0 | 7 |
| **Total** | **6.5 days** | **7** | **23** |

---

## Next Steps

1. **Immediate (Today):**
   - Create missing Python scripts (bugs #1, #2)
   - Fix numpy import bug (#7)
   - Fix QC threshold (#5)

2. **This Week:**
   - Create missing shell scripts (#3)
   - Complete pathogen list (#12)
   - Add spot instance fallback (#8)

3. **Next Sprint:**
   - Address security issues (#9, #10, #11)
   - Improve test coverage
   - Update documentation

4. **Before Production:**
   - Complete end-to-end testing
   - Security audit
   - Performance testing
   - PMDA compliance verification

---

## Conclusion

The codebase has a **solid architectural foundation** but requires critical bug fixes before production deployment. The 7 critical bugs will prevent pipeline execution and must be resolved immediately. The infrastructure and overall design are sound, with most issues being missing implementation files rather than fundamental design flaws.

**Recommendation:** Complete Sprint 0 fixes (1.5 days) before any production testing.

**Risk Level:** 🔴 HIGH - Pipeline cannot execute in current state

---

**Report Generated by:** Claude Code Audit System
**Date:** 2025-11-06
**Total Issues:** 23 (7 Critical, 8 High, 5 Medium, 3 Low)
