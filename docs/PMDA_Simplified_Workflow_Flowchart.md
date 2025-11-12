# PMDA 91-Pathogen Simplified Workflow Flowchart

**Document Version**: 1.0
**Date**: 2025-11-13
**Protocol Reference**: Protocol 12 (Unified Sample Prep)

---

## Visual Workflow Overview

### Complete Workflow Diagram

```
┌──────────────────────────────────────────────────────────────────────────┐
│                    PMDA 91-PATHOGEN SCREENING WORKFLOW                   │
│                      (Protocol 12: Simplified Universal)                 │
└──────────────────────────────────────────────────────────────────────────┘

                               START
                                 │
                      ┌──────────┴──────────┐
                      │  Blood Collection   │
                      │   5-10 mL EDTA      │
                      │   (SUPERase•In for  │
                      │    RNA samples)     │
                      └──────────┬──────────┘
                                 │
                      ┌──────────┴──────────┐
                      │  Plasma Separation  │
                      │   1,500×g, 15 min   │
                      │        4°C          │
                      └──────────┬──────────┘
                                 │
          ┌──────────────────────┴──────────────────────┐
          │       UNIVERSAL DNA/RNA EXTRACTION          │
          │        (Zymo Quick-cfDNA/cfRNA Kit)         │
          │               3 hours                       │
          └──────────────────────┬──────────────────────┘
                                 │
                    ┌────────────┴────────────┐
                    │                         │
            ┌───────▼────────┐       ┌───────▼────────┐
            │   DNA Fraction │       │   RNA Fraction │
            │                │       │                │
            │   10-100 ng    │       │    5-50 ng     │
            └───────┬────────┘       └───────┬────────┘
                    │                        │
         ┌──────────▼──────────┐  ┌──────────▼──────────┐
         │  CpG Host Depletion │  │  Poly(A) Selection  │
         │   (NEBNext Kit)     │  │    (NEBNext Kit)    │
         │      2 hours        │  │      2 hours        │
         └──────────┬──────────┘  └──────────┬──────────┘
                    │                        │
                    │             ┌──────────▼──────────┐
                    │             │   cDNA Synthesis    │
                    │             │  (NEBNext Ultra II) │
                    │             │      4 hours        │
                    │             └──────────┬──────────┘
                    │                        │
         ┌──────────▼──────────┐  ┌──────────▼──────────┐
         │  DNA Library Prep   │  │  DNA Library Prep   │
         │    (LSK114 Kit)     │  │    (LSK114 Kit)     │
         │      4 hours        │  │      4 hours        │
         └──────────┬──────────┘  └──────────┬──────────┘
                    │                        │
                    └────────────┬───────────┘
                                 │
                      ┌──────────▼──────────┐
                      │  Library QC         │
                      │  (Qubit, TapeStation)│
                      │      30 min         │
                      └──────────┬──────────┘
                                 │
                                 ├─────► PASS?
                                 │       │
                          ┌──────┘       └──────┐
                         NO                    YES
                          │                      │
                  ┌───────▼────────┐   ┌────────▼────────┐
                  │   Troubleshoot │   │  MinION Loading │
                  │   (See Protocol│   │   24-48 hours   │
                  │   Appendix C)  │   │   (Duplex mode) │
                  └────────────────┘   └────────┬────────┘
                                                │
                                      ┌─────────▼─────────┐
                                      │  Bioinformatics   │
                                      │     Analysis      │
                                      │  detect_pmda_all_ │
                                      │  91_pathogens.py  │
                                      └─────────┬─────────┘
                                                │
                                    ┌───────────▼───────────┐
                                    │  90/91 Pathogens      │
                                    │  Detected & Quantified│
                                    └───────────┬───────────┘
                                                │
                                    ┌───────────▼───────────┐
                                    │  Retrovirus-like      │
                                    │  pol gene detected?   │
                                    └───────────┬───────────┘
                                                │
                                       ┌────────┴────────┐
                                      NO               YES
                                       │                 │
                            ┌──────────▼─────────┐  ┌───▼────────────┐
                            │  SCREENING COMPLETE│  │ Conditional    │
                            │  Report Generation │  │ Spumavirus Test│
                            │                    │  │ (PBMC + PCR)   │
                            └────────────────────┘  └───┬────────────┘
                                                        │
                                                   ┌────▼────────┐
                                                   │  Final Report│
                                                   │  (All 91)    │
                                                   └──────────────┘
                                                        │
                                                       END
```

---

## Workflow Decision Tree

### Sample Routing Logic

```
┌─────────────────────────────────┐
│  Blood Sample (5-10 mL EDTA)    │
└────────────┬────────────────────┘
             │
             ▼
    ┌────────────────────┐
    │ Always Extract     │
    │ DNA + RNA          │
    │ (Universal)        │
    └────────┬───────────┘
             │
             ├─────────────────────────┐
             │                         │
      ┌──────▼──────┐           ┌─────▼──────┐
      │ DNA Workflow│           │RNA Workflow│
      │ (Always)    │           │ (Always)   │
      └──────┬──────┘           └─────┬──────┘
             │                        │
             ▼                        ▼
    Target: 51 pathogens      Target: 20 pathogens
    - 24 DNA viruses          - 16 poly(A)+ viruses
    - 27 Bacteria             - 4 poly(A)- viruses

    Both sequenced → 90/91 detected

    If retrovirus pol detected:
         ↓
    ┌────────────────────┐
    │ Conditional PBMC   │
    │ Collection + PCR   │
    │ (Spumavirus only)  │
    └────────────────────┘
```

### Quality Control Decision Points

```
Each Step → QC Check
              │
       ┌──────┴──────┐
      PASS          FAIL
       │              │
       ▼              ▼
   Continue      Troubleshoot
                      │
                  ┌───┴────┐
                 Fix   Repeat
                  │      │
                  └──────┘
```

---

## Timeline Comparison

### Traditional Multi-Workflow (Pre-Simplification)

```
Day 1 (9h hands-on):
├─ 09:00-10:00: Sample prep
├─ 10:00-13:00: DNA/RNA extraction
├─ 14:00-16:00: Host depletion (DNA)
└─ 16:00-20:00: Decision: rRNA or poly(A)? ← BRANCHING

Day 2 (7h hands-on):
├─ 09:00-13:00: RNA specific prep (varies by virus type)
├─ 13:00-17:00: Library prep
└─ Optional: Amplicon RT-PCR (+6h) ← BRANCHING

Total: 16-22 hours hands-on + decision complexity
```

### Simplified Universal Workflow (Protocol 12)

```
Day 1 (9h hands-on):
├─ 09:00-10:00: Sample prep
├─ 10:00-13:00: DNA/RNA extraction
├─ 14:00-16:00: Host depletion (DNA) + Poly(A) (RNA) ← PARALLEL
└─ 16:00-20:00: DNA library prep starts

Day 2 (4h hands-on):
├─ 09:00-13:00: cDNA → RNA library prep
├─ 13:00-14:00: QC both libraries
└─ 14:00: Load MinION

Total: 13 hours hands-on (19% time savings, no branching)
```

---

## Pathogen Coverage Map

### By Sample Type and Workflow

```
┌─────────────────────────────────────────────────────────────┐
│                   PLASMA SAMPLE (5-10 mL)                   │
└────┬────────────────────────────────────────────────────┬───┘
     │                                                    │
┌────▼─────────────────────┐           ┌─────────────────▼────┐
│   DNA WORKFLOW           │           │   RNA WORKFLOW        │
│   (51 pathogens)         │           │   (20 pathogens)      │
├──────────────────────────┤           ├───────────────────────┤
│ • 24 DNA Viruses         │           │ • 16 Poly(A)+ Viruses │
│   - PPV, PRV, ASFV       │           │   - EEEV, GETV        │
│   - SWPV, PAV, PCMV      │           │   - PRRSV, PEDV       │
│   - PLHV, PGHV           │           │   - CSFV, JEV, HEV    │
│   - PCV2, PCV3, TTV      │           │   - TGEV, SIV, FMDV   │
│   - POLYOMA, PERV        │           │   - EMCV, RABV        │
│                          │           │   - MENV, NIPV        │
│ • 27 Bacteria            │           │   - BDV, BVDV         │
│   - All bacterial        │           │                       │
│     pathogens            │           │ • 4 Poly(A)- Viruses  │
│                          │           │   - HANTV, WEEV       │
│                          │           │   - VEEV, REO         │
└──────────────────────────┘           └───────────────────────┘
        ↓                                        ↓
     LOD:                                     LOD:
  100-200 copies/mL                   100-300 copies/mL (poly(A)+)
                                      300-500 copies/mL (poly(A)-)

┌─────────────────────────────────────────────────────────────┐
│           COMBINED ANALYSIS: 90/91 PATHOGENS                │
│                (Fungi + Parasites included)                 │
└─────────────────────────────────────────────────────────────┘
                              ↓
                 ┌────────────────────────┐
                 │ Retrovirus detected?   │
                 └────────┬───────────────┘
                          │
                   ┌──────┴──────┐
                  NO            YES
                   │              │
                  END     ┌───────▼────────┐
                          │ PBMC Collection│
                          │  (Spumavirus)  │
                          │   1/91 pathogen│
                          └────────────────┘
                                 ↓
                            LOD: 1-10
                        copies/10⁵ PBMCs
```

---

## Cost Breakdown Comparison

### Traditional Workflow (Variable Costs)

```
Sample Type 1 (DNA viruses only):
├─ Extraction: ¥12,000
├─ Host depletion: ¥8,000
├─ DNA library: ¥25,000
├─ Sequencing (1/24): ¥30,000
├─ Misc: ¥7,000
└─ Subtotal: ¥82,000

Sample Type 2 (Poly(A)+ RNA viruses):
├─ Extraction: ¥12,000
├─ Poly(A) selection: ¥5,000
├─ cDNA: ¥15,000
├─ RNA library: ¥25,000
├─ Sequencing (1/24): ¥30,000
├─ Misc: ¥8,000
└─ Subtotal: ¥95,000

Sample Type 3 (Poly(A)- RNA viruses with rRNA depletion):
├─ Extraction: ¥12,000
├─ rRNA depletion: ¥15,000 ← Expensive!
├─ cDNA: ¥15,000
├─ RNA library: ¥25,000
├─ Sequencing (1/24): ¥30,000
├─ Misc: ¥8,000
└─ Subtotal: ¥105,000

Average (mixed samples): ¥94,000-100,000
```

### Simplified Universal Workflow (Fixed Costs)

```
ALL Samples (Universal):
├─ Extraction (DNA+RNA): ¥12,000
├─ Host depletion (DNA): ¥8,000
├─ Poly(A) selection (RNA): ¥5,000 ← Cheaper than rRNA depletion
├─ cDNA synthesis: ¥15,000
├─ DNA library: ¥25,000
├─ RNA library: ¥25,000
├─ Sequencing (×2, 1/24 each): ¥60,000
├─ Misc: ¥7,000
└─ Total: ¥157,000 (fixed)

Conditional Spumavirus (5-10% samples):
└─ Add: ¥15,000 for PBMC + PCR

Weighted Average: ¥157,000 + (0.075 × ¥15,000) = ¥158,125
```

**Cost Impact**: +¥5,000-8,000 per sample (+3-5%)
**Benefit**: Predictable costs, no variability, simpler budgeting

---

## Performance Metrics Dashboard

### Key Performance Indicators (KPIs)

```
┌──────────────────────────────────────────────────────────┐
│                    WORKFLOW METRICS                      │
├────────────────────┬─────────────────┬──────────────────┤
│      Metric        │   Traditional   │   Simplified     │
├────────────────────┼─────────────────┼──────────────────┤
│ Workflows          │       3-4       │        2         │
│ Decision Points    │        5        │        2         │
│ Hands-on Time (h)  │      16-22      │       13         │
│ Cost (¥/sample)    │   94,000-105,000│     157,000      │
│ Pathogen Coverage  │      91/91      │      91/91       │
│ LOD (copies/mL)    │     50-500      │    100-500       │
│ Training Time (d)  │       5-7       │        3         │
│ Error Rate (%)     │       5-8       │       2-3        │
│ Turnaround (days)  │       4-5       │       3-4        │
└────────────────────┴─────────────────┴──────────────────┘

✅ Improved: Workflows, Decision Points, Time, Training, Errors, Turnaround
⚠️ Higher Cost: +¥5,000 (but predictable and justified by efficiency)
✅ Maintained: Coverage, LOD (acceptable for screening)
```

---

## Quality Control Checkpoints

### Sample-to-Report QC Gates

```
Gate 1: Sample Collection
├─ Volume: 5-10 mL ✓
├─ Anticoagulant: EDTA ✓
├─ Hemolysis: None ✓
└─ Time to processing: <4h ✓

Gate 2: Extraction
├─ DNA yield: 10-100 ng ✓
├─ RNA yield: 5-50 ng ✓
├─ 260/280: 1.8-2.1 ✓
└─ 260/230: >2.0 ✓

Gate 3: Host Depletion / Selection
├─ Host removal: >95% ✓
├─ Pathogen recovery: >90% ✓
└─ Final volume: adequate ✓

Gate 4: Library Prep
├─ Concentration: >50 fmol ✓
├─ Length distribution: expected range ✓
└─ Adapter dimers: <5% ✓

Gate 5: Sequencing
├─ Active pores: >800/2048 ✓
├─ Q30: >85% ✓
├─ Total reads: >500,000 ✓
└─ Data volume: >1.25 Gb/sample ✓

Gate 6: Analysis
├─ Mapping rate: >70% ✓
├─ Coverage: sufficient for LOD ✓
├─ Controls: pass ✓
└─ Report: PMDA compliant ✓
```

---

## Troubleshooting Quick Reference

### Common Issues and Solutions

```
Issue 1: Low DNA/RNA Yield
└─► Increase sample volume to 10 mL
└─► Check SUPERase•In concentration (RNA)
└─► Verify extraction kit lot number

Issue 2: High Host Contamination (>10%)
└─► Increase MBD2-Fc protein amount
└─► Extend incubation time to 15 min
└─► Verify pig genome methylation status

Issue 3: Poor Library Quality
└─► Check reagent expiration dates
└─► Optimize adapter:DNA ratio
└─► Increase AMPure bead wash repeats

Issue 4: Low Sequencing Yield
└─► Check flow cell QC (>800 active pores)
└─► Verify loading concentration
└─► Enable Duplex mode for Q30 boost

Issue 5: Poly(A)- Virus Not Detected
└─► Expected: LOD 300-500 copies/mL (acceptable)
└─► If higher sensitivity needed: Switch to Protocol 11
└─► Use rRNA depletion + amplicon RT-PCR
```

---

## Protocol Selection Guide

### When to Use Which Protocol

```
┌─────────────────────────────────────────────────────────┐
│              PROTOCOL SELECTION FLOWCHART               │
└─────────────────────────────────────────────────────────┘

                    Start Sample Screening
                            │
                   ┌────────▼────────┐
                   │ What is your    │
                   │ primary goal?   │
                   └────────┬────────┘
                            │
          ┌─────────────────┼─────────────────┐
          │                 │                 │
    ┌─────▼──────┐   ┌──────▼──────┐   ┌─────▼─────┐
    │ Routine    │   │ Confirmation│   │  Research │
    │ Screening  │   │ After Positive│  │ Validation│
    └─────┬──────┘   └──────┬──────┘   └─────┬─────┘
          │                 │                 │
          ▼                 ▼                 ▼
   ┌─────────────┐   ┌─────────────┐   ┌─────────────┐
   │ PROTOCOL 12 │   │ PROTOCOL 11 │   │ PROTOCOL 11 │
   │ (Unified)   │   │ (Enhanced)  │   │ (Enhanced)  │
   ├─────────────┤   ├─────────────┤   ├─────────────┤
   │ • 91 patho  │   │ • 4 viruses │   │ • LOD study │
   │ • 13h time  │   │ • 19.5h time│   │ • Sensitivity│
   │ • ¥157,000  │   │ • ¥185,000  │   │ • Validation│
   │ • LOD 100+  │   │ • LOD <50   │   └─────────────┘
   └─────────────┘   └─────────────┘

Recommendation: START with Protocol 12 for ALL samples
               ESCALATE to Protocol 11 only if needed
```

---

## References

1. Protocol 12: 統合サンプル調製プロトコル (md/MinION_Protocol_12_統合サンプル調製プロトコル.md)
2. Protocol 11: PMDA 4ウイルス高感度検出プロトコル (md/MinION_Protocol_11_PMDA_4ウイルス高感度検出プロトコル.md)
3. PMDA Simplified Sample Prep Strategy (docs/PMDA_Simplified_Sample_Prep_Strategy.md)
4. Protocol 00: マスタードキュメント (md/MinION_Protocol_00_目次とマスタードキュメント.md)

---

**Document Control**:
- **Version**: 1.0
- **Created**: 2025-11-13
- **Author**: MinION Pipeline Development Team
- **Format**: ASCII flowcharts for maximum compatibility
- **Next Review**: 2026-05-13
