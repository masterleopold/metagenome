# LOD Validation Protocol: PCV2, PCV3, TTV, PPV Detection with Step 2.5

**Protocol Type**: Analytical Validation Study
**Study Phase**: Method Validation
**Target Analytes**: PCV2, PCV3, TTV, PPV (Circular & Linear ssDNA Viruses)
**Study Duration**: 3-6 months (3 months core + 3 months extended)
**Version**: 1.0
**Date**: 2025-11-13
**Status**: Ready for Implementation

---

## Table of Contents

1. [Study Objectives](#study-objectives)
2. [Regulatory Context](#regulatory-context)
3. [Study Design Overview](#study-design-overview)
4. [Materials and Reagents](#materials-and-reagents)
5. [Sample Preparation](#sample-preparation)
6. [Experimental Design](#experimental-design)
7. [Data Collection](#data-collection)
8. [Statistical Analysis](#statistical-analysis)
9. [Success Criteria](#success-criteria)
10. [Timeline and Resources](#timeline-and-resources)
11. [Appendices](#appendices)

---

## Study Objectives

### Primary Objective
Determine the **Limit of Detection (LOD)** for 4 PMDA-designated pathogens using Protocol 12 v2.1 with Step 2.5 (circular DNA linearization + ssDNA→dsDNA conversion).

**Target LOD**: 100-500 copies/mL plasma for all 4 pathogens

### Secondary Objectives
1. Confirm **≥95% detection rate** at target LOD (PMDA PPA requirement)
2. Assess **day-to-day reproducibility** (inter-day CV <20%)
3. Evaluate **operator-to-operator reproducibility** (inter-operator CV <20%)
4. Determine **LOD95** (concentration with 95% detection probability)
5. Validate **quantification accuracy** (R² >0.90 vs. known input)

### Tertiary Objectives
1. Compare **with vs. without Step 2.5** (demonstrate necessity)
2. Assess **circular vs. linear ssDNA** performance differences
3. Document **failure modes** for troubleshooting guide

---

## Regulatory Context

### PMDA Requirements
- **異種移植指針 別添2**: All 91 designated pathogens must be detectable
- **PPA (Positive Percent Agreement)**: ≥95% at screening concentration
- **NPA (Negative Percent Agreement)**: ≥98% (false positive rate <2%)
- **LOD documentation**: Serial dilution studies with ≥20 replicates per concentration

### Target Pathogens (Special Considerations)

| Pathogen | PMDA Classification | Genome Structure | Size | Clinical Significance |
|----------|---------------------|------------------|------|----------------------|
| **PCV2** | Special Management #3 | Circular ssDNA | 1.7 kb | Post-weaning multisystemic wasting syndrome |
| **PCV3** | Special Management #3 | Circular ssDNA | 2.0 kb | Cardiac/multisystemic inflammation |
| **TTV** | Standard #40 | Circular ssDNA | 3.8 kb | Ubiquitous, low pathogenicity |
| **PPV** | Standard #1 | Linear ssDNA | 5.0 kb | Reproductive failure, fetal infection |

**Note**: PCV2 and PCV3 are **Special Management Pathogens** requiring enhanced surveillance.

---

## Study Design Overview

### Three-Phase Approach

**Phase 1: Core LOD Determination** (3 months)
- Serial dilution studies (50-2000 copies/mL, 8 concentrations)
- 20 replicates per concentration per pathogen
- 3 operators × 10 days = 30 independent runs
- Total samples: 4 pathogens × 8 concentrations × 20 replicates = **640 samples**

**Phase 2: Reproducibility & Robustness** (2 months)
- Inter-day reproducibility: 10 consecutive days, 3 replicates/day
- Inter-operator reproducibility: 3 operators, 10 replicates each
- Robustness testing: Variable conditions (±20% reagent variation)
- Total samples: **240 samples**

**Phase 3: Comparative & Extended Studies** (1 month, optional)
- With vs. without Step 2.5 comparison
- Mixed pathogen samples (multiple ssDNA viruses)
- Plasma matrix effects (different pig sources)
- Total samples: **120 samples**

**Grand Total**: 640 + 240 + 120 = **1000 samples** (full validation)

---

## Materials and Reagents

### Reference Standards (Critical)

**Option 1: Synthetic DNA (Recommended for Safety)**
- **Vendor**: Twist Bioscience or Integrated DNA Technologies (IDT)
- **Design**: Full-length genome sequences in appropriate topology
  - PCV2: Circular ssDNA construct, 1767 bp (GenBank AF027217)
  - PCV3: Circular ssDNA construct, 2000 bp (GenBank MF318988)
  - TTV: Circular ssDNA construct, 3800 bp (GenBank AB017610)
  - PPV: Linear ssDNA fragment, 5000 bp (GenBank NC_001718)
- **Quality control**: Sanger sequencing verification (full length)
- **Quantification**: Qubit dsDNA HS + ddPCR absolute quantification
- **Cost**: ~¥400,000 for 4 constructs (5 μg each)

**Option 2: Natural Isolates (Higher Biosafety Requirements)**
- **Source**: Pig plasma samples with confirmed PCV2/PCV3/TTV/PPV infection
- **Requirement**: BSL-2+ facility, approved by institutional biosafety committee
- **Quantification**: RT-qPCR with plasmid standards (triplicate)
- **Risk**: Natural variation, potential co-infections

**Recommendation**: Use **Option 1 (Synthetic DNA)** for controlled, reproducible validation.

---

### Reagents and Kits

**Core Protocol 12 v2.1 Reagents**:
- Zymo Quick-cfDNA/cfRNA Serum & Plasma Kit (50 preps) - ¥400,000
- NEBNext Microbiome DNA Enrichment Kit (48 rxn) - ¥384,000
- **Step 2.5 Reagents** (for 640 samples):
  - DNase I (RNase-free, NEB M0303, 1000 U) × 5 - ¥60,000
  - Klenow Fragment (exo-, NEB M0212, 500 U) × 5 - ¥75,000
  - Random Hexamers (NEB S1230, 300 μg) × 3 - ¥30,000
  - dNTP Set (NEB N0446, 8 μmol each) × 5 - ¥25,000
- SQK-LSK114 Ligation Sequencing Kit (10 kits) - ¥250,000
- Flow Cells R9.4.1 (54 cells) - ¥6,750,000

**QC Reagents**:
- Qubit dsDNA HS Assay Kit (500 assays) - ¥45,000
- Agilent High Sensitivity DNA Kit (25 chips) - ¥125,000
- AMPure XP beads (450 mL) - ¥180,000

**Plasma Matrix**:
- Pathogen-free pig plasma (SPF Yucatan, 2 L) - ¥100,000

**Total Reagent Cost**: ~¥8,824,000 (~¥8.8M for 1000-sample validation)

---

### Equipment Requirements

**Essential Equipment**:
- MinION Mk1D device (1 unit) - Available
- PCR machines with gradient function (2 units) - Available
- Magnetic separation racks (96-well format) - Available
- Qubit 4 Fluorometer - Available
- Bioanalyzer 2100 - Available
- Vortex, centrifuges, pipettes (P2-P1000) - Available
- -80°C freezer (sample storage) - Available
- -20°C freezer (reagent storage) - Available

**Optional Equipment**:
- ddPCR system (QX200, for absolute quantification) - ¥15,000,000
- Fragment Analyzer (QC automation) - ¥20,000,000

---

## Sample Preparation

### Reference Standard Preparation

**Step 1: Stock Solution Preparation** (Day -30)

1. **Resuspend synthetic DNA constructs**:
   - PCV2 circular ssDNA (5 μg lyophilized)
   - PCV3 circular ssDNA (5 μg lyophilized)
   - TTV circular ssDNA (5 μg lyophilized)
   - PPV linear ssDNA (5 μg lyophilized)

2. **Initial quantification** (Triplicate):
   - Qubit dsDNA HS assay
   - Bioanalyzer DNA HS chip
   - Calculate molecular weight: MW = (size in bp) × 650 Da/bp

3. **Calculate copy number**:
   ```
   copies/μL = (ng/μL × 6.022×10²³) / (MW in Da × 10⁹)

   Example for PCV2 (1767 bp, 10 ng/μL):
   MW = 1767 × 650 = 1,148,550 Da
   copies/μL = (10 × 6.022×10²³) / (1.14855×10⁶ × 10⁹)
              = 5.24×10⁹ copies/μL
   ```

4. **Prepare working stocks**:
   - Dilute to 1×10⁷ copies/μL in TE buffer (10 mM Tris, 1 mM EDTA, pH 8.0)
   - Aliquot into 50 μL single-use tubes
   - Store at -80°C (stable for 2 years)
   - Avoid freeze-thaw cycles (max 2 cycles)

**Step 2: Plasma Matrix Spiking** (Each experiment day)

1. **Thaw pathogen-free pig plasma** (on ice, 30 min)

2. **Thaw reference standard aliquot** (on ice, 15 min)

3. **Prepare serial dilutions in plasma**:
   - Concentration range: 50, 100, 200, 300, 500, 1000, 1500, 2000 copies/mL
   - Volume: 1 mL per concentration per pathogen
   - Mix thoroughly: vortex 5 sec, invert 10×

4. **Prepare negative control**: Pathogen-free plasma (no spike)

5. **Prepare positive control**: 1000 copies/mL (mid-range)

**Step 3: Extraction and Processing**

Follow Protocol 12 v2.1:
1. cfDNA extraction (Zymo Kit, 3h)
2. CpG host depletion (NEBNext, 2h)
3. **Step 2.5**: Circular DNA linearization + ssDNA→dsDNA conversion (2.5h)
4. DNA library preparation (LSK114, 4h)
5. MinION sequencing (24h)

---

## Experimental Design

### Phase 1: Core LOD Determination (3 Months)

#### Concentration Series Design

| Concentration (copies/mL) | Expected Detection Rate | Rationale | Replicates |
|---------------------------|-------------------------|-----------|------------|
| 2000 | 100% | Well above LOD | 20 |
| 1500 | 100% | Above LOD | 20 |
| 1000 | 100% | Target upper LOD | 20 |
| 500 | ≥95% | **Target LOD** | 20 |
| 300 | 80-95% | Near LOD | 20 |
| 200 | 60-80% | Below LOD | 20 |
| 100 | 40-60% | Well below LOD | 20 |
| 50 | <40% | Limit of blank | 20 |
| 0 (Negative control) | 0% | Specificity check | 20 |

**Total**: 9 concentrations × 20 replicates × 4 pathogens = **720 samples**

#### Experimental Schedule (12 Weeks)

**Week 1-3: PCV2 (Circular ssDNA, 1.7 kb)**
- 3 operators (Operator A, B, C)
- Each operator: 60 samples (9 concentrations × 6-7 replicates)
- 5 samples/day × 4 days/week = 20 samples/week
- Total: 180 samples over 3 weeks

**Week 4-6: PCV3 (Circular ssDNA, 2.0 kb)**
- Same schedule as PCV2
- Total: 180 samples over 3 weeks

**Week 7-9: TTV (Circular ssDNA, 3.8 kb)**
- Same schedule as PCV2
- Total: 180 samples over 3 weeks

**Week 10-12: PPV (Linear ssDNA, 5.0 kb)**
- Same schedule as PCV2
- Total: 180 samples over 3 weeks

**Phase 1 Total**: 720 samples, 12 weeks

---

### Phase 2: Reproducibility & Robustness (2 Months)

#### 2.1 Inter-Day Reproducibility (4 Weeks)

**Objective**: Confirm CV <20% across multiple days

**Design**:
- **Concentration**: 500 copies/mL (target LOD)
- **Schedule**: 10 consecutive working days
- **Replicates**: 3 per day per pathogen
- **Operator**: Single operator (Operator A)
- **Total**: 10 days × 3 replicates × 4 pathogens = **120 samples**

**Metrics**:
- Mean detection rate per day
- Inter-day CV% for detection rate
- Inter-day CV% for quantification (if detected)

---

#### 2.2 Inter-Operator Reproducibility (4 Weeks)

**Objective**: Confirm CV <20% across different operators

**Design**:
- **Concentration**: 500 copies/mL (target LOD)
- **Operators**: 3 independent operators (A, B, C)
- **Replicates**: 10 per operator per pathogen
- **Total**: 3 operators × 10 replicates × 4 pathogens = **120 samples**

**Metrics**:
- Mean detection rate per operator
- Inter-operator CV% for detection rate
- Inter-operator CV% for quantification

---

#### 2.3 Robustness Testing (2 Weeks, Optional)

**Objective**: Assess sensitivity to variable conditions

**Variables Tested** (±20% from nominal):
1. DNase I concentration: 0.004, 0.005, 0.006 U
2. Linearization time: 4, 5, 6 minutes
3. Klenow incubation time: 48, 60, 72 minutes
4. Random hexamer concentration: 4, 5, 6 μM

**Design**:
- **Concentration**: 500 copies/mL (PCV2 only, as representative)
- **Conditions**: 4 variables × 3 levels = 12 conditions
- **Replicates**: 5 per condition
- **Total**: 12 conditions × 5 replicates = **60 samples**

**Success Criteria**: LOD remains 100-500 copies/mL across all conditions

---

### Phase 3: Comparative & Extended Studies (1 Month, Optional)

#### 3.1 With vs. Without Step 2.5 Comparison

**Objective**: Demonstrate necessity of Step 2.5

**Design**:
- **Concentration**: 500 copies/mL
- **Conditions**:
  - WITH Step 2.5 (linearization + second-strand synthesis)
  - WITHOUT Step 2.5 (direct to library prep)
- **Pathogens**: All 4 (PCV2, PCV3, TTV, PPV)
- **Replicates**: 10 per condition per pathogen
- **Total**: 2 conditions × 10 replicates × 4 pathogens = **80 samples**

**Expected Result**:
- WITH Step 2.5: ≥95% detection
- WITHOUT Step 2.5: <5% detection (demonstrates >19× improvement)

---

#### 3.2 Mixed Pathogen Samples

**Objective**: Assess simultaneous detection of multiple ssDNA viruses

**Design**:
- **Spike combinations**:
  - PCV2 + PCV3 (500 copies/mL each)
  - PCV2 + TTV (500 copies/mL each)
  - PCV3 + PPV (500 copies/mL each)
  - All 4 pathogens (500 copies/mL each)
- **Replicates**: 5 per combination
- **Total**: 4 combinations × 5 replicates = **20 samples**

**Success Criteria**: Detection rate ≥95% for each pathogen in mixture

---

#### 3.3 Plasma Matrix Effects

**Objective**: Assess LOD variability across different plasma sources

**Design**:
- **Plasma sources**: 3 different SPF pig colonies
- **Concentration**: 500 copies/mL (PCV2 only)
- **Replicates**: 5 per source
- **Total**: 3 sources × 5 replicates = **15 samples**

**Success Criteria**: Detection rate ≥95% across all sources (no matrix interference)

---

## Data Collection

### Primary Data Points (Per Sample)

**Pre-Sequencing QC**:
1. cfDNA extraction yield (ng)
2. Post-CpG depletion yield (ng)
3. Post-Step 2.5.1 (linearization) yield (ng)
4. Post-Step 2.5.2 (second-strand synthesis) yield (ng)
5. Library concentration (ng/μL)
6. Bioanalyzer size distribution (bp range)

**Sequencing Metrics**:
7. Total reads generated
8. Q30 percentage (%)
9. Mean read length (bp)
10. Active pores during run

**Detection Metrics**:
11. Pathogen detected (Yes/No)
12. Read count mapping to pathogen genome
13. Genome coverage (%)
14. Mean sequencing depth (×)
15. Quantification (copies/mL, back-calculated)

### Data Recording Template

**Sample Processing Record**:
```
Sample ID: LODVAL-[Pathogen]-[Conc]-[Rep]-[Date]
Example: LODVAL-PCV2-0500-R01-20251113

Operator: _______________
Date: _______________
Spike Concentration (nominal): _______ copies/mL
Plasma Source: _______________

Extraction & Processing:
├─ cfDNA yield (ng): _______
├─ Post-CpG yield (ng): _______
├─ Post-linearization yield (ng): _______
├─ Post-synthesis yield (ng): _______
├─ Library concentration (ng/μL): _______
└─ Bioanalyzer peak size (bp): _______

Sequencing:
├─ Run ID: _______________
├─ Flow cell ID: _______________
├─ Total reads: _______
├─ Q30%: _______
├─ Mean read length (bp): _______
└─ Active pores: _______

Detection:
├─ Pathogen detected: ☐ Yes  ☐ No
├─ Read count: _______
├─ Genome coverage (%): _______
├─ Mean depth (×): _______
└─ Quantification (copies/mL): _______

QC Status: ☐ Pass  ☐ Fail
Operator Signature: _______________
Reviewer Signature: _______________
```

---

## Statistical Analysis

### LOD Calculation Methods

#### Method 1: Probit Analysis (Preferred)

**Principle**: Model detection probability as function of concentration using probit regression

**Software**: R (package: `MASS`), GraphPad Prism, or SPSS

**Steps**:
1. Code detection: 1 = detected, 0 = not detected
2. Fit probit model: Φ⁻¹(P) = β₀ + β₁ × log₁₀(concentration)
3. Calculate LOD95: Concentration with 95% detection probability
4. Calculate 95% confidence interval using likelihood ratio test

**Output**:
- LOD95 (copies/mL) ± 95% CI
- Probit plot (S-curve)

**Example R code**:
```R
library(MASS)

# Data: concentration (copies/mL), detected (1/0)
data <- data.frame(
  concentration = rep(c(50, 100, 200, 300, 500, 1000, 1500, 2000), each=20),
  detected = c(...)  # Binary detection results
)

# Fit probit model
model <- glm(detected ~ log10(concentration), family=binomial(link="probit"), data=data)

# Calculate LOD95
LOD95 <- 10^((qnorm(0.95) - coef(model)[1]) / coef(model)[2])
print(paste("LOD95 =", round(LOD95, 0), "copies/mL"))

# 95% CI using confint()
ci <- confint(model)
```

---

#### Method 2: Hit Rate Analysis (Alternative)

**Principle**: Empirical hit rate at each concentration

**Steps**:
1. Calculate hit rate: (# detected) / (# total) at each concentration
2. Identify lowest concentration with hit rate ≥95%
3. Calculate 95% CI using Clopper-Pearson exact method

**Example**:
```
Concentration: 500 copies/mL
Replicates: 20
Detected: 19
Hit rate: 19/20 = 95%
95% CI: [75.1%, 99.9%] (Clopper-Pearson)
```

**Advantage**: Simple, no assumptions about concentration-response relationship
**Disadvantage**: Requires exact 95% threshold, less precise than probit

---

### Reproducibility Analysis

#### Inter-Day CV% Calculation

**Formula**:
```
CV% = (SD / Mean) × 100%

Where:
SD = √[Σ(xi - x̄)² / (n-1)]
x̄ = Mean hit rate across days
n = Number of days (10)
```

**Example**:
```
Day 1: 100% (3/3)
Day 2: 100% (3/3)
Day 3: 67% (2/3)
Day 4: 100% (3/3)
...
Day 10: 100% (3/3)

Mean = 93.3%
SD = 10.5%
CV% = (10.5 / 93.3) × 100% = 11.3% ✓ (Pass: <20%)
```

---

#### Inter-Operator CV% Calculation

**Formula**: Same as inter-day CV%

**Example**:
```
Operator A: 90% (9/10)
Operator B: 100% (10/10)
Operator C: 90% (9/10)

Mean = 93.3%
SD = 5.8%
CV% = (5.8 / 93.3) × 100% = 6.2% ✓ (Pass: <20%)
```

---

### Quantification Accuracy Analysis

**Objective**: Assess R² >0.90 (linearity) between observed and expected concentrations

**Method**: Linear regression
```
Y = β₀ + β₁X

Where:
Y = Observed copies/mL (back-calculated from read counts)
X = Expected copies/mL (nominal spike concentration)
```

**Metrics**:
- R² (coefficient of determination): Target ≥0.90
- Slope: Target 0.8-1.2 (within ±20% of 1.0)
- Intercept: Target close to 0
- Residual plot: No systematic bias

---

## Success Criteria

### Primary Success Criteria (Must Pass All)

| Criterion | Target | Analysis Method | Pass/Fail |
|-----------|--------|-----------------|-----------|
| **LOD95 for PCV2** | 100-500 copies/mL | Probit analysis | ☐ Pass ☐ Fail |
| **LOD95 for PCV3** | 100-500 copies/mL | Probit analysis | ☐ Pass ☐ Fail |
| **LOD95 for TTV** | 100-500 copies/mL | Probit analysis | ☐ Pass ☐ Fail |
| **LOD95 for PPV** | 100-500 copies/mL | Probit analysis | ☐ Pass ☐ Fail |
| **Detection rate at 500 copies/mL** | ≥95% (all 4 pathogens) | Hit rate | ☐ Pass ☐ Fail |
| **Specificity (negative control)** | 0% false positive | Hit rate | ☐ Pass ☐ Fail |

---

### Secondary Success Criteria (Recommended)

| Criterion | Target | Analysis Method | Pass/Fail |
|-----------|--------|-----------------|-----------|
| **Inter-day CV%** | <20% | ANOVA | ☐ Pass ☐ Fail |
| **Inter-operator CV%** | <20% | ANOVA | ☐ Pass ☐ Fail |
| **Quantification R²** | ≥0.90 | Linear regression | ☐ Pass ☐ Fail |
| **Quantification slope** | 0.8-1.2 | Linear regression | ☐ Pass ☐ Fail |

---

### Tertiary Success Criteria (Optional)

| Criterion | Target | Analysis Method | Pass/Fail |
|-----------|--------|-----------------|-----------|
| **Robustness (variable conditions)** | LOD remains 100-500 | Probit analysis | ☐ Pass ☐ Fail |
| **With vs. without Step 2.5** | >19× improvement | Hit rate comparison | ☐ Pass ☐ Fail |
| **Mixed pathogen detection** | ≥95% for each | Hit rate | ☐ Pass ☐ Fail |
| **Matrix effect (3 sources)** | ≥95% across sources | Hit rate | ☐ Pass ☐ Fail |

---

### Decision Tree

```
┌─────────────────────────────────────────────────┐
│  Validation Decision Tree                       │
└─────────────────────────────────────────────────┘

All 4 pathogens: LOD95 = 100-500 copies/mL?
├─ YES → All 4 pathogens: Detection ≥95% at 500 copies/mL?
│        ├─ YES → Specificity 0% (no false positives)?
│        │        ├─ YES → ✅ PRIMARY CRITERIA MET
│        │        │        └─ Proceed to secondary criteria
│        │        └─ NO  → ❌ FAIL: Optimize specificity
│        └─ NO  → ❌ FAIL: Optimize sensitivity
└─ NO  → ❌ FAIL: Re-optimize Step 2.5 parameters

Secondary Criteria:
├─ Inter-day CV% <20%?
├─ Inter-operator CV% <20%?
├─ R² ≥0.90?
└─ Slope 0.8-1.2?

All YES → ✅ VALIDATION COMPLETE → Deploy to production
Any NO  → ⚠️ WARNING: Document deviations, deploy with caution
```

---

## Timeline and Resources

### Phase 1: Core LOD Determination (12 Weeks)

**Personnel**:
- 3 laboratory technicians (full-time, 12 weeks)
- 1 bioinformatics analyst (50% FTE, 12 weeks)
- 1 project manager (25% FTE, 12 weeks)

**Workload**:
- 720 samples ÷ 12 weeks = 60 samples/week
- 60 samples ÷ 3 operators = 20 samples/operator/week
- 20 samples ÷ 5 days = 4 samples/operator/day

**Cost**:
- Reagents: ¥8,800,000
- Personnel (3 technicians × 12 weeks × ¥50,000/week): ¥1,800,000
- Bioinformatics (1 analyst × 6 weeks × ¥80,000/week): ¥480,000
- Project management (1 PM × 3 weeks × ¥100,000/week): ¥300,000
- **Phase 1 Total**: ¥11,380,000

---

### Phase 2: Reproducibility (8 Weeks)

**Personnel**:
- 3 laboratory technicians (full-time, 8 weeks)
- 1 bioinformatics analyst (50% FTE, 8 weeks)
- 1 statistician (25% FTE, 8 weeks)

**Workload**:
- 240 samples ÷ 8 weeks = 30 samples/week
- Lighter workload allows for detailed QC and documentation

**Cost**:
- Reagents: Already included in Phase 1 budget
- Personnel (3 technicians × 8 weeks × ¥50,000/week): ¥1,200,000
- Bioinformatics (1 analyst × 4 weeks × ¥80,000/week): ¥320,000
- Statistics (1 statistician × 2 weeks × ¥80,000/week): ¥160,000
- **Phase 2 Total**: ¥1,680,000

---

### Phase 3: Comparative Studies (4 Weeks, Optional)

**Personnel**:
- 2 laboratory technicians (full-time, 4 weeks)
- 1 bioinformatics analyst (50% FTE, 4 weeks)

**Workload**:
- 120 samples ÷ 4 weeks = 30 samples/week

**Cost**:
- Reagents: Already included in Phase 1 budget
- Personnel (2 technicians × 4 weeks × ¥50,000/week): ¥400,000
- Bioinformatics (1 analyst × 2 weeks × ¥80,000/week): ¥160,000
- **Phase 3 Total**: ¥560,000

---

### Grand Total Budget

| Phase | Duration | Samples | Cost |
|-------|----------|---------|------|
| Phase 1 (Core LOD) | 12 weeks | 720 | ¥11,380,000 |
| Phase 2 (Reproducibility) | 8 weeks | 240 | ¥1,680,000 |
| Phase 3 (Comparative) | 4 weeks | 120 | ¥560,000 |
| **Total** | **24 weeks (6 months)** | **1080** | **¥13,620,000** |

**Note**: Reagent costs dominate (¥8.8M / ¥13.6M = 65%)

---

### Gantt Chart

```
Month 1:  [========PCV2========]
Month 2:  [========PCV3========]
Month 3:  [========TTV=========]
Month 4:  [========PPV=========]
Month 5:  [==Reproducibility===]
Month 6:  [==Robustness==][Comparative]
```

---

## Appendices

### Appendix A: Sample Size Justification

**Power Analysis for LOD95 Estimation**:
- Target precision: ±100 copies/mL (95% CI width <200 copies/mL)
- Expected hit rate at LOD95: 95%
- Probit analysis requires ≥15 replicates per concentration for robust estimation
- **Selected: 20 replicates** (provides 33% margin for sample failures)

**Regulatory Precedent**:
- FDA Guidance: ≥20 replicates at LOD (Bacteriological Analytical Manual)
- CLSI EP17-A2: ≥20 replicates for LOD confirmation studies

---

### Appendix B: Reference Sequences

**PCV2 (1767 bp, Circular ssDNA)**:
- GenBank: AF027217.1
- Order from: Twist Bioscience or IDT
- Topology: Circular construct (plasmid backbone with PCV2 insert)
- QC: Sanger sequencing primers spanning entire genome

**PCV3 (2000 bp, Circular ssDNA)**:
- GenBank: MF318988.1
- Order from: Twist Bioscience or IDT
- Topology: Circular construct

**TTV (3800 bp, Circular ssDNA)**:
- GenBank: AB017610.1 (prototype isolate)
- Order from: Twist Bioscience or IDT
- Topology: Circular construct

**PPV (5000 bp, Linear ssDNA)**:
- GenBank: NC_001718.1 (Kresse strain)
- Order from: IDT (linear gBlock)
- Topology: Linear fragment (no circularization needed)

---

### Appendix C: Data Analysis Scripts

**R Script for Probit Analysis**:
```R
# Load required libraries
library(MASS)
library(ggplot2)

# Load data (CSV format: Sample_ID, Concentration, Detected)
data <- read.csv("LOD_validation_PCV2.csv")

# Fit probit model
model <- glm(Detected ~ log10(Concentration),
             family = binomial(link = "probit"),
             data = data)

# Summary
summary(model)

# Calculate LOD95
LOD95 <- 10^((qnorm(0.95) - coef(model)[1]) / coef(model)[2])
cat("LOD95 =", round(LOD95, 0), "copies/mL\n")

# 95% Confidence interval
ci <- confint(model)
LOD95_lower <- 10^((qnorm(0.95) - ci[1,2]) / ci[2,2])
LOD95_upper <- 10^((qnorm(0.95) - ci[1,1]) / ci[2,1])
cat("95% CI: [", round(LOD95_lower, 0), ",", round(LOD95_upper, 0), "]\n")

# Plot probit curve
concentrations <- seq(log10(50), log10(2000), length.out=100)
predictions <- predict(model,
                       newdata = data.frame(Concentration = 10^concentrations),
                       type = "response")

plot_data <- data.frame(
  Concentration = 10^concentrations,
  Probability = predictions
)

ggplot() +
  geom_point(data = data, aes(x = Concentration, y = Detected), alpha=0.3) +
  geom_line(data = plot_data, aes(x = Concentration, y = Probability), color="blue", size=1.2) +
  geom_hline(yintercept = 0.95, linetype="dashed", color="red") +
  geom_vline(xintercept = LOD95, linetype="dashed", color="red") +
  scale_x_log10() +
  labs(title = "PCV2 LOD Probit Curve",
       x = "Concentration (copies/mL, log scale)",
       y = "Detection Probability") +
  theme_minimal()

ggsave("PCV2_LOD_probit_curve.png", width=8, height=6, dpi=300)
```

---

### Appendix D: Report Template

**LOD Validation Study Report**

**Title**: Limit of Detection Validation for PCV2, PCV3, TTV, PPV Using Protocol 12 v2.1

**Study Period**: [Start Date] to [End Date]

**1. Executive Summary**
- Study objectives
- Key findings (LOD95 for each pathogen)
- Success criteria met/not met
- Recommendations

**2. Introduction**
- Regulatory context (PMDA requirements)
- Protocol 12 v2.1 overview
- Study scope and limitations

**3. Materials and Methods**
- Reference standards (source, preparation, QC)
- Sample preparation (spiking, extraction, Step 2.5)
- Experimental design (Phase 1-3 details)
- Statistical analysis methods

**4. Results**
- **4.1 Phase 1: LOD Determination**
  - Probit curves for each pathogen
  - LOD95 ± 95% CI (table)
  - Hit rate tables
  - Representative sequencing data
- **4.2 Phase 2: Reproducibility**
  - Inter-day CV% (table, ANOVA)
  - Inter-operator CV% (table, ANOVA)
  - Quantification accuracy (scatter plot, R²)
- **4.3 Phase 3: Comparative Studies** (if performed)
  - With vs. without Step 2.5 comparison
  - Mixed pathogen results
  - Matrix effect assessment

**5. Discussion**
- Interpretation of results
- Comparison with target criteria
- Limitations and future work

**6. Conclusions**
- Summary statement on validation success
- Protocol readiness for production deployment

**7. Appendices**
- Raw data tables
- QC plots (Bioanalyzer traces)
- Operator competency records
- Deviation logs

---

### Appendix E: Risk Assessment

**Potential Risks and Mitigation**:

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| **Reference standard degradation** | Medium | High | Aliquot single-use tubes, limit freeze-thaw, QC at each use |
| **Operator error (pipetting)** | Medium | Medium | Competency training, duplicate measurements |
| **Contamination** | Low | High | Dedicated workspace, negative controls, regular cleaning |
| **Flow cell failure** | Medium | Medium | Purchase extra flow cells (10% buffer), QC before use |
| **Bioinformatics pipeline error** | Low | Medium | Validate pipeline with known standards, version control |
| **Sample mix-up** | Low | High | Barcode system, double-check labeling, independent verification |

---

## Document Control

**Version**: 1.0
**Date**: 2025-11-13
**Author**: MinION Pipeline Development Team
**Approved by**: [Pending]
**Review Date**: Upon completion of Phase 1 (3 months)

**Revision History**:
| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | 2025-11-13 | Initial protocol | Pipeline Team |

---

**End of Validation Protocol**
