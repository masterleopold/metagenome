# Step 2.5 Staff Training Guide: Circular & ssDNA Virus Detection

**Document Type**: Staff Training Guide
**Target Audience**: Laboratory Technicians, Research Scientists
**Training Duration**: 1 day (8 hours: 3h theory + 5h hands-on)
**Version**: 1.0
**Date**: 2025-11-13
**Related Protocol**: MinION Protocol 12 v2.1, Step 2.5

---

## Training Objectives

By the end of this training, participants will be able to:

1. **Explain** why circular and ssDNA viruses require special handling (theoretical understanding)
2. **Identify** the 4 PMDA pathogens requiring Step 2.5 (PCV2, PCV3, TTV, PPV)
3. **Perform** DNase I ultra-low concentration linearization independently
4. **Execute** Klenow Fragment second-strand synthesis with >95% success rate
5. **Troubleshoot** common issues (no linearization, incomplete synthesis, over-digestion)
6. **Interpret** QC results (Bioanalyzer traces, Qubit quantification)
7. **Document** Step 2.5 in compliance with ALCOA+ principles

---

## Prerequisites

### Required Knowledge
- ✅ Basic molecular biology (DNA structure, PCR, gel electrophoresis)
- ✅ Experience with AMPure XP bead purification
- ✅ Familiarity with Protocol 12 Steps 1-2 (extraction, CpG depletion)
- ✅ Qubit and Bioanalyzer operation

### Required Certifications
- ✅ BSL-2 laboratory safety training
- ✅ Good Laboratory Practice (GLP) training
- ✅ ALCOA+ data integrity training

### Materials for Training
- 1× Protocol 12 reagent kit (Step 2.5 reagents)
- Positive control: Circular plasmid DNA (3 kb, 100 ng/μL)
- Negative control: Nuclease-free water
- Mock plasma sample (pig plasma spiked with PCV2 synthetic DNA)

---

## Training Schedule (8 Hours)

### Morning Session (4 hours): Theory + Demonstration

#### 09:00-10:00 | Module 1: Theoretical Background (1h)
**Learning Objective**: Understand why Step 2.5 is necessary

**Topics Covered**:
1. **PMDA 91-pathogen requirement**
   - 厚労省別添2指定病原体リスト
   - PPA >95% requirement for ALL pathogens
   - Special Management pathogens (PCV2, PCV3)

2. **DNA genome structure diversity**
   - Linear dsDNA (most DNA viruses)
   - Circular dsDNA (Adenovirus, SV40)
   - Circular ssDNA (Circovirus, TTV)
   - Linear ssDNA (Parvovirus)

3. **SQK-LSK114 ligation kit limitations**
   - NEBNext Ultra II End-Repair requires free 3'/5' termini
   - T4 DNA Ligase: >95% efficiency on dsDNA, <5% on ssDNA
   - MinION motor proteins require dsDNA input

4. **Why circular DNA cannot be sequenced without linearization**
   ```
   Circular DNA (e.g., PCV2 1.7 kb):

        ╔════════════════╗
        ║  No free ends  ║  ← Cannot attach adapters
        ╚════════════════╝

   After DNase I linearization:

        5'═══════════════3'
        3'═══════════════5'  ← Free ends available
   ```

5. **Why ssDNA cannot be sequenced without second-strand synthesis**
   ```
   Single-stranded DNA (e.g., PPV 5.0 kb):

        5'─────────────────3'  ← T4 Ligase <5% efficiency

   After Klenow Fragment synthesis:

        5'═══════════════3'
        3'═══════════════5'  ← T4 Ligase >95% efficiency
   ```

**Assessment**: 5-question quiz (passing score: 4/5)

---

#### 10:00-11:00 | Module 2: Step 2.5 Protocol Walk-Through (1h)
**Learning Objective**: Understand each sub-step in detail

**Sub-step 2.5.1: Circular DNA Linearization (30 min total)**

**Principle**: DNase I creates random single-strand nicks in circular DNA, converting closed circles to linear molecules.

**Critical Parameters**:
- **DNase I concentration**: 0.005 U (ultra-low, 1/2000 of standard PCR cleanup)
- **Temperature**: 37°C (optimal enzyme activity)
- **Time**: 5 minutes (limited digestion, not complete degradation)
- **Inactivation**: 75°C, 10 minutes (irreversible heat denaturation)

**Step-by-Step**:
1. Prepare linearization master mix (on ice):
   - DNA sample (20 μL from CpG depletion)
   - DNase I Buffer (10×): 3 μL
   - Nuclease-free water: 6.5 μL
   - DNase I (0.005 U): 0.5 μL
   - **Total volume**: 30 μL

2. Mix gently (pipette 5× without bubbles)

3. Incubate 37°C, 5 minutes (PCR machine or heat block)

4. Heat inactivation 75°C, 10 minutes

5. Cool to room temperature (2 min on bench)

6. AMPure XP purification (0.8× ratio):
   - Add 24 μL beads (0.8× of 30 μL)
   - Mix 10× by pipetting
   - Incubate RT, 5 min
   - Magnet 5 min until clear
   - Remove supernatant (save for troubleshooting)
   - Wash 2× with 200 μL 80% EtOH
   - Air dry 5 min (until matte finish, no cracks)
   - Elute in 21 μL nuclease-free water
   - Final volume: 20 μL linearized DNA

**Expected Result**:
- Qubit: 70-90% recovery (if input was 10 ng, expect 7-9 ng)
- Bioanalyzer: Slight size shift (circular → linear mobility change)

---

**Sub-step 2.5.2: ssDNA → dsDNA Conversion (2h total)**

**Principle**: Random hexamer primers bind ssDNA templates, Klenow Fragment synthesizes complementary strand.

**Critical Parameters**:
- **Klenow Fragment (exo-)**: 3'→5' exonuclease-deficient (prevents template degradation)
- **Random hexamers**: 50 μM stock, 5 μM final (saturating concentration)
- **Denaturation**: 95°C, 3 min (denature any secondary structure)
- **Annealing**: 25°C, 5 min (hexamer binding)
- **Extension**: 16°C, 60 min (slow, processivity ~1 kb/min)

**Step-by-Step**:

**Phase 1: Denaturation (5 min)**
1. Transfer 20 μL linearized DNA to PCR tube
2. Denature: 95°C, 3 min (PCR machine)
3. Immediately transfer to ice (snap-cool, 2 min)
   - **Critical**: Prevents re-annealing of complementary regions

**Phase 2: Annealing (10 min)**
4. Prepare annealing master mix (on ice):
   - Denatured DNA: 20 μL
   - Random Hexamers (50 μM): 2 μL
   - **Total**: 22 μL

5. Mix gently, spin down

6. Anneal: 25°C, 5 min (heat block or PCR machine)
   - **Result**: Hexamers bind to ssDNA template every ~100-200 bp

**Phase 3: Extension (90 min)**
7. Prepare extension master mix (keep at RT):
   - Annealed DNA + hexamers: 22 μL
   - NEBNext Ultra II Reaction Buffer (10×): 3 μL
   - dNTP mix (10 mM each): 3 μL
   - Nuclease-free water: 1.5 μL
   - Klenow Fragment (exo-, 5 U/μL): 0.5 μL
   - **Total**: 30 μL

8. Mix gently (pipette 5× without bubbles)

9. Incubate 16°C, 60 min (PCR machine)
   - **Why 16°C?** Balances primer stability (higher Tm) with Klenow activity

10. Heat inactivation 75°C, 10 min

11. Cool to room temperature (2 min on bench)

**Phase 4: Purification (30 min)**
12. AMPure XP purification (1.8× ratio):
    - Add 54 μL beads (1.8× of 30 μL)
    - Mix 10× by pipetting
    - Incubate RT, 5 min
    - Magnet 5 min until clear
    - Remove supernatant (save for troubleshooting)
    - Wash 2× with 200 μL 80% EtOH
    - Air dry 5 min (until matte finish)
    - Elute in 21 μL nuclease-free water
    - Final volume: 20 μL dsDNA

**Expected Result**:
- Qubit: 80-95% recovery (if input was 8 ng, expect 6.4-7.6 ng)
- Bioanalyzer: Peak size unchanged (size remains same, but now dsDNA)
- **Key indicator**: Bioanalyzer concentration approximately 2× Qubit (dsDNA intercalates more dye)

---

#### 11:00-12:00 | Module 3: Quality Control & Troubleshooting (1h)
**Learning Objective**: Interpret QC results and troubleshoot failures

**QC Checkpoints**:

**Checkpoint 1: After Sub-step 2.5.1 (Linearization)**
- **Qubit dsDNA HS**: 70-90% recovery
  - ✅ PASS: 7-9 ng (if input was 10 ng)
  - ⚠️ WARNING: 5-7 ng (marginal recovery, may proceed)
  - ❌ FAIL: <5 ng (50% loss, troubleshoot before continuing)

- **Bioanalyzer DNA HS**:
  - ✅ PASS: Single peak, no smearing, size ~expected
  - ❌ FAIL: Smear or multiple peaks (over-digestion)

**Checkpoint 2: After Sub-step 2.5.2 (Second-Strand Synthesis)**
- **Qubit dsDNA HS**: 80-95% recovery
  - ✅ PASS: 6.4-7.6 ng (if input was 8 ng)
  - ❌ FAIL: <5 ng (re-do synthesis)

- **Bioanalyzer DNA HS**:
  - ✅ PASS: Bioanalyzer concentration ~2× Qubit (dsDNA characteristic)
  - ❌ FAIL: Bioanalyzer = Qubit (still ssDNA, incomplete synthesis)

---

**Common Issues & Solutions**:

| Issue | Cause | Solution |
|-------|-------|----------|
| **<50% recovery after linearization** | DNase I over-digestion | Reduce DNase I to 0.003 U or time to 3 min |
| **Bioanalyzer shows smear** | DNase I over-digestion | Use fresh DNase I, reduce concentration |
| **No linearization (circular DNA remains)** | DNase I inactive | Check DNase I storage (-20°C), use fresh aliquot |
| **Bioanalyzer = Qubit after synthesis** | Incomplete second-strand synthesis | Increase Klenow to 1 U or extend to 90 min |
| **<50% recovery after synthesis** | Klenow exonuclease contamination | Use Klenow (exo-) only, not regular Klenow |
| **No product after synthesis** | Random hexamers degraded | Use fresh random hexamers, store -20°C |

---

#### 12:00-13:00 | Lunch Break

---

### Afternoon Session (4 hours): Hands-On Practice

#### 13:00-15:30 | Module 4: Hands-On Practice Session 1 (2.5h)
**Learning Objective**: Perform Step 2.5 with positive control

**Sample**: Circular plasmid DNA (3 kb, pUC19 or similar, 10 ng input)

**Activities**:
1. **Pre-lab setup (15 min)**:
   - Calculate reagent volumes
   - Prepare ice bucket
   - Set up PCR machine programs
   - Label tubes

2. **Sub-step 2.5.1 execution (45 min)**:
   - Prepare linearization master mix
   - Incubate 37°C, 5 min
   - Heat inactivation 75°C, 10 min
   - AMPure XP purification (0.8×)
   - Qubit measurement
   - Record results

3. **Sub-step 2.5.2 execution (90 min)**:
   - Denaturation (95°C, 3 min → ice)
   - Annealing with random hexamers (25°C, 5 min)
   - Extension with Klenow (16°C, 60 min)
   - Heat inactivation (75°C, 10 min)
   - AMPure XP purification (1.8×)
   - Qubit measurement
   - Record results

**Expected Results**:
- Post-linearization: 7-9 ng (70-90% recovery)
- Post-synthesis: 6-8 ng (80-95% recovery)

**Trainer Observation Checklist**:
- ☐ Correct pipetting technique (no bubbles)
- ☐ Proper ice usage (reagents kept cold)
- ☐ Accurate timing (no delays between steps)
- ☐ Bead handling (complete EtOH removal, proper drying)
- ☐ Documentation (ALCOA+ compliant)

---

#### 15:30-15:45 | Break (15 min)

---

#### 15:45-17:00 | Module 5: Hands-On Practice Session 2 (1.25h)
**Learning Objective**: Perform Step 2.5 with mock plasma sample

**Sample**: Pig plasma spiked with PCV2 synthetic circular ssDNA (500 copies/mL, post-CpG depletion)

**Activities**:
1. **Independent execution** (trainee performs alone, trainer observes):
   - Complete Sub-step 2.5.1 (linearization)
   - Complete Sub-step 2.5.2 (second-strand synthesis)
   - QC measurements (Qubit)

2. **Results interpretation**:
   - Compare pre/post-Step 2.5 concentrations
   - Assess recovery percentages
   - Determine PASS/FAIL based on QC criteria

**Success Criteria for Competency**:
- ✅ Both QC checkpoints PASS
- ✅ Recovery >70% for linearization
- ✅ Recovery >80% for second-strand synthesis
- ✅ Documentation complete and accurate
- ✅ No major procedural errors

---

#### 17:00-17:30 | Module 6: Competency Assessment (30 min)
**Learning Objective**: Demonstrate mastery of Step 2.5

**Assessment Components**:

**1. Written Test (15 min, 10 questions)**
- Why can't circular DNA be sequenced directly? (1 pt)
- What is the DNase I concentration in Sub-step 2.5.1? (1 pt)
- Why use Klenow Fragment (exo-) instead of regular Klenow? (1 pt)
- What is the expected recovery after linearization? (1 pt)
- What is the expected recovery after second-strand synthesis? (1 pt)
- Which 4 PMDA pathogens require Step 2.5? (1 pt)
- What QC metric indicates successful dsDNA synthesis? (1 pt)
- Troubleshooting: Bioanalyzer shows smear after linearization. Cause? (1 pt)
- Troubleshooting: Bioanalyzer = Qubit after synthesis. Cause? (1 pt)
- Why is Step 2.5 critical for PMDA compliance? (1 pt)

**Passing score**: 8/10

**2. Practical Competency Observation (15 min)**
Trainer observes trainee perform one complete cycle:
- ☐ Correct reagent preparation (no calculation errors)
- ☐ Proper pipetting technique
- ☐ Accurate timing and temperature control
- ☐ Correct AMPure XP bead handling
- ☐ Proper documentation

**Passing criteria**: All 5 items checked

---

## Competency Certification

**Trainee Name**: _______________________________
**Training Date**: _______________________________
**Trainer Name**: _______________________________

### Assessment Results

**Written Test Score**: ______ / 10 (Passing: ≥8)
**Practical Observation**: ☐ Pass  ☐ Fail

### Certification Decision

☐ **CERTIFIED**: Authorized to perform Step 2.5 independently
☐ **REMEDIATION REQUIRED**: Additional training needed in areas:
  - ☐ Theoretical understanding
  - ☐ Pipetting technique
  - ☐ Timing/temperature control
  - ☐ QC interpretation
  - ☐ Troubleshooting
  - ☐ Documentation

**Trainer Signature**: _______________________________
**Date**: _______________________________

**Trainee Signature**: _______________________________
**Date**: _______________________________

---

## Appendix A: Pre-Training Checklist

**Equipment Preparation** (Complete 1 day before training):
- ☐ PCR machine programmed with Step 2.5 protocols
- ☐ Qubit calibrated with standards
- ☐ Magnetic rack available (8-tube or 96-well)
- ☐ Pipettes calibrated (P2, P10, P20, P200)
- ☐ Ice machine functional

**Reagent Preparation**:
- ☐ DNase I (RNase-free, NEB M0303) thawed, aliquoted (0.1 U aliquots)
- ☐ Klenow Fragment (exo-, NEB M0212) thawed, on ice
- ☐ Random Hexamers (NEB S1230) thawed, aliquoted (10 μL aliquots)
- ☐ dNTP Set (NEB N0446) thawed, on ice
- ☐ AMPure XP beads equilibrated to room temperature (30 min)
- ☐ 80% EtOH freshly prepared (50 mL per trainee)
- ☐ Nuclease-free water (multiple tubes)

**Sample Preparation**:
- ☐ Positive control: Circular plasmid DNA (pUC19, 10 ng/μL, 50 μL)
- ☐ Mock sample: PCV2-spiked plasma (post-CpG depletion, 50 μL)
- ☐ Negative control: Nuclease-free water (50 μL)

**Documentation**:
- ☐ Training manual printed (1 copy per trainee)
- ☐ Blank worksheets prepared
- ☐ Competency assessment forms printed

---

## Appendix B: Step 2.5 Quick Reference Card

**For Laminated Bench Card** (double-sided):

### SIDE 1: Sub-step 2.5.1 - Linearization (30 min)

```
┌─────────────────────────────────────────────┐
│  STEP 2.5.1: CIRCULAR DNA LINEARIZATION     │
└─────────────────────────────────────────────┘

INPUT: 20 μL DNA (post-CpG depletion)

MASTER MIX (30 μL total):
├─ DNA sample: 20 μL
├─ DNase I Buffer (10×): 3 μL
├─ Nuclease-free H₂O: 6.5 μL
└─ DNase I (0.005 U): 0.5 μL

PROTOCOL:
1. Mix on ice (gentle, 5× pipette)
2. Incubate 37°C, 5 min
3. Heat inactivation 75°C, 10 min
4. Cool to RT (2 min bench)

PURIFICATION (AMPure XP 0.8×):
1. Add 24 μL beads, mix 10×
2. RT 5 min, magnet 5 min
3. Remove supernatant (SAVE!)
4. Wash 2× with 200 μL 80% EtOH
5. Air dry 5 min (matte finish)
6. Elute 21 μL H₂O, collect 20 μL

QC CHECKPOINT:
✓ Qubit: 70-90% recovery
✓ Bioanalyzer: Single peak, no smear

OUTPUT: 20 μL linearized DNA → Sub-step 2.5.2
```

### SIDE 2: Sub-step 2.5.2 - Second-Strand Synthesis (2h)

```
┌─────────────────────────────────────────────┐
│  STEP 2.5.2: ssDNA → dsDNA CONVERSION       │
└─────────────────────────────────────────────┘

INPUT: 20 μL linearized DNA

PHASE 1: DENATURATION
1. 95°C, 3 min (PCR machine)
2. Ice immediately, 2 min

PHASE 2: ANNEALING (22 μL total)
├─ Denatured DNA: 20 μL
└─ Random Hexamers (50 μM): 2 μL
→ 25°C, 5 min

PHASE 3: EXTENSION (30 μL total)
├─ Annealed mix: 22 μL
├─ Reaction Buffer (10×): 3 μL
├─ dNTP mix (10 mM): 3 μL
├─ Nuclease-free H₂O: 1.5 μL
└─ Klenow (exo-): 0.5 μL
→ 16°C, 60 min
→ 75°C, 10 min (inactivate)

PURIFICATION (AMPure XP 1.8×):
1. Add 54 μL beads, mix 10×
2. RT 5 min, magnet 5 min
3. Remove supernatant (SAVE!)
4. Wash 2× with 200 μL 80% EtOH
5. Air dry 5 min (matte finish)
6. Elute 21 μL H₂O, collect 20 μL

QC CHECKPOINT:
✓ Qubit: 80-95% recovery
✓ Bioanalyzer ~2× Qubit (dsDNA indicator)

OUTPUT: 20 μL dsDNA → DNA Library Prep (Step 3)
```

---

## Appendix C: Troubleshooting Decision Tree

```
┌─────────────────────────────────────────────┐
│  QC FAIL: What to do?                       │
└─────────────────────────────────────────────┘

After Sub-step 2.5.1 (Linearization):
├─ <50% recovery?
│  ├─ Bioanalyzer shows smear? → DNase I over-digestion
│  │  └─ ACTION: Reduce DNase I to 0.003 U OR time to 3 min
│  └─ Bioanalyzer shows intact circular? → DNase I inactive
│     └─ ACTION: Use fresh DNase I aliquot, check storage
│
└─ 50-70% recovery? → Marginal, may proceed
   └─ Monitor next steps closely

After Sub-step 2.5.2 (Second-Strand Synthesis):
├─ <50% recovery?
│  └─ ACTION: Repeat synthesis with fresh reagents
│
├─ Bioanalyzer = Qubit (not ~2×)?
│  └─ Incomplete synthesis → ACTION: Extend time to 90 min
│     OR increase Klenow to 1 U
│
└─ 50-80% recovery? → Marginal, may proceed to library prep
   └─ Adjust library prep input volume if needed
```

---

## Appendix D: ALCOA+ Documentation Template

**Step 2.5 Execution Record**

**Sample ID**: _______________  **Date**: _______________  **Operator**: _______________

### Sub-step 2.5.1: Linearization

| Parameter | Value | Operator Initials | Time |
|-----------|-------|-------------------|------|
| Input DNA (ng) | | | |
| DNase I lot # | | | |
| Incubation start (37°C) | | | |
| Inactivation start (75°C) | | | |
| Post-purification Qubit (ng) | | | |
| Recovery (%) | | | |
| QC Result (Pass/Fail) | | | |

### Sub-step 2.5.2: Second-Strand Synthesis

| Parameter | Value | Operator Initials | Time |
|-----------|-------|-------------------|------|
| Input DNA (ng) | | | |
| Klenow lot # | | | |
| Random hexamer lot # | | | |
| Denaturation start (95°C) | | | |
| Extension start (16°C) | | | |
| Post-purification Qubit (ng) | | | |
| Recovery (%) | | | |
| QC Result (Pass/Fail) | | | |

**Overall Step 2.5 Status**: ☐ PASS  ☐ FAIL  ☐ REPEAT

**Deviations/Notes**: _________________________________________________________________

**Operator Signature**: _______________________________  **Date**: _______________

**Reviewer Signature**: _______________________________  **Date**: _______________

---

## Document Control

**Version**: 1.0
**Date**: 2025-11-13
**Author**: MinION Pipeline Development Team
**Approved by**: [Pending]
**Next Review**: 2026-05-13 (6 months post-deployment)

**Revision History**:
| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | 2025-11-13 | Initial version | Pipeline Team |

---

**End of Training Guide**
