# Quick Fill-In Template for Grant Application

Use this template to gather all the information you need to fill in the `[TO BE FILLED]` placeholders in the main grant application.

## Section 1: Principal Investigator Information

```
Full Name: _________________________________
Title/Position: _________________________________
Department: _________________________________
Institution: _________________________________
Email: _________________________________
Phone: _________________________________
ORCID (if applicable): _________________________________
```

## Section 2: Co-Investigators (Optional)

```
Co-Investigator #1:
Name: _________________________________
Institution: _________________________________
Expertise: _________________________________

Co-Investigator #2:
Name: _________________________________
Institution: _________________________________
Expertise: _________________________________
```

## Section 3: Graduate Students

```
Student #1:
Name: _________________________________
Program: (PhD / Master's)
Focus Area: _________________________________

Student #2:
Name: _________________________________
Program: (PhD / Master's)
Focus Area: _________________________________
```

## Section 4: Institutional Information

```
University/Institution Full Name: Meiji University (明治大学)
Department Full Name: _________________________________
Address: _________________________________
City, State, ZIP: _________________________________
Country: Japan

Sponsored Programs Office:
Contact Name: _________________________________
Email: _________________________________
Phone: _________________________________
```

## Section 5: PI Qualifications (Section 6 of proposal)

**Copy-paste this into the proposal Section 6.2:**

```markdown
Dr. [YOUR NAME] has [X]+ years of experience in [YOUR FIELD]. Their lab focuses on [YOUR RESEARCH FOCUS], with expertise in:

1. **[Area 1]:** [Brief description + 1-2 key publications]
2. **[Area 2]:** [Brief description + 1-2 key publications]
3. **[Area 3]:** [Brief description + 1-2 key publications]
4. **[Area 4]:** [Brief description + 1-2 key publications]

Key Publications (selected from [X]+ total):
- [Author] et al. (20XX). "[Title]." [Journal].
- [Author] et al. (20XX). "[Title]." [Journal].
- [Author] et al. (20XX). "[Title]." [Journal].
- [Author] et al. (20XX). "[Title]." [Journal].
- [Author] et al. (20XX). "[Title]." [Journal].

Funding History:
- [Agency] [Mechanism]: "[Project Title]" ($[Amount], [Years])
- [Agency] [Mechanism]: "[Project Title]" ($[Amount], [Years])

NVIDIA Collaborations:
- [Describe any previous NVIDIA collaborations, OR write:]
  "This represents our first formal collaboration with NVIDIA, though we have extensively used NVIDIA GPUs (Tesla K80, V100, A100) in our research computing infrastructure."
```

## Section 6: Your Most Relevant Publications

List your 5-10 most relevant publications for the References section:

```
1. _________________________________
2. _________________________________
3. _________________________________
4. _________________________________
5. _________________________________
6. _________________________________
7. _________________________________
8. _________________________________
9. _________________________________
10. _________________________________
```

## Section 7: Institutional Resources (Section 6.4)

```markdown
[UNIVERSITY NAME] provides:

1. **Research Computing:**
   - High-Performance Computing cluster with [X]+ NVIDIA GPUs ([Types: V100, A100, etc.])
   - [X] TB shared storage for genomics data
   - [Speed] network connectivity

2. **Biosafety Infrastructure:**
   - BSL-[1/2/3] laboratory for sample processing
   - MinION sequencing facility ([X]× devices)
   - [Relevant research programs]

3. **Administrative Support:**
   - Sponsored Programs Office (grant management)
   - Technology Transfer Office (IP, commercialization)
   - Research Data Management services

4. **Institutional Commitment:**
   - Letter of support from [Name, Title] (attached)
   - Cost-share commitment: $[Amount] for [personnel/equipment/etc.]
   - Space allocation: [X] sq ft laboratory
```

## Section 8: Institutional Cost-Share (if applicable)

```
[UNIVERSITY NAME] will provide cost-share for:

1. Personnel ([%] of [X] graduate student RA): $[Amount]/year × [Years] = $[Total]
2. Sequencing consumables (MinION flow cells): $[Amount]
3. Laboratory space and utilities: $[Amount]
4. Data storage (local NAS for backups): $[Amount]

Total institutional cost-share: $[TOTAL]
```

## Section 9: Collaborators

```
External Collaborator #1:
Organization: _________________________________
Contact Person: _________________________________
Role in Project: _________________________________

External Collaborator #2:
Organization: _________________________________
Contact Person: _________________________________
Role in Project: _________________________________
```

## Section 10: Diversity and Inclusion Statement (Section 8.4)

```markdown
We are committed to promoting diversity in computational biology:

1. **Recruitment:** Actively recruit students from underrepresented groups through partnerships with [INSTITUTION's specific programs: e.g., McNair Scholars, LSAMP, etc.]

2. **International Collaboration:** Partner with universities in [specific countries/regions] to broaden participation.

3. **Accessibility:** All tutorials and materials will be:
   - Free and open-source
   - Available in [languages]
   - Designed for learners with limited computational resources

4. **Mentorship:** PI has mentored [X]+ diverse students, with [Y]% advancing to PhD programs or industry positions.
```

## Submission Checklist

Once you've filled in all the above information:

- [ ] **Replace all `[TO BE FILLED]` in main proposal** with info from this template
- [ ] **Section 1 (Cover Sheet):** PI and institutional info complete
- [ ] **Section 6.2 (PI Qualifications):** CV summary, publications, funding added
- [ ] **Section 6.3 (Team):** Co-investigators and students listed
- [ ] **Section 6.4 (Resources):** Institutional support detailed
- [ ] **Section 7.4 (Cost-Share):** Budget breakdown provided (if applicable)
- [ ] **Section 8.4 (Diversity):** Diversity statement customized
- [ ] **Section 9 (References):** Your publications added to reference list

## Additional Documents to Prepare

### Document 1: PI CV (2 pages max)

**Format:** NSF or NIH style

**Sections to include:**
- Professional Preparation (degrees, institutions, years)
- Appointments (current and past positions)
- Products (5 publications most relevant to this proposal)
- Synergistic Activities (teaching, service, outreach)
- Collaborators & Co-Editors (past 48 months)

**Where to get template:**
- NSF: https://www.nsf.gov/bfa/dias/policy/biosketch.jsp
- NIH: https://grants.nih.gov/grants/forms/biosketch.htm

### Document 2: Letter of Institutional Support

**Request from:** Department Head or Dean

**Should include:**
- Statement of support for your research
- Commitment to resources (space, equipment, personnel)
- Cost-share commitment (if any)
- Alignment with institutional priorities
- Signature on official letterhead

**Timeline:** Request 2-3 weeks before submission deadline (they're busy!)

### Document 3: Budget Spreadsheet (Optional)

If NVIDIA requests a separate budget file, create Excel with:

| Item | Unit | Quantity | Cost/Unit | Total | Notes |
|------|------|----------|-----------|-------|-------|
| A100 GPU Hours - Basecalling | hours | 700 | [Rate] | [Total] | 350 samples × 2h |
| A100 GPU Hours - Pathogen | hours | 850 | [Rate] | [Total] | 350 samples × 2.5h |
| A100 GPU Hours - AI Training | hours | 500 | [Rate] | [Total] | Model training |
| A100 GPU Hours - Benchmarking | hours | 250 | [Rate] | [Total] | Testing |
| A100 GPU Hours - Buffer | hours | 200 | [Rate] | [Total] | 8% contingency |
| **TOTAL A100 HOURS** | | **2,500** | | **[TOTAL]** | |
| DGX Spark Systems | units | 2 | $3,999 | $7,998 | Hardware |
| **GRAND TOTAL** | | | | **[TOTAL]** | |

## Quick Reference: Where Info Goes in Proposal

| Your Info | Goes in Proposal Section |
|-----------|-------------------------|
| PI name, title, email | Section 1 (Cover Sheet) |
| Co-investigators | Section 6.3 (Team Composition) |
| PI CV summary | Section 6.2 (PI Qualifications) |
| Your publications | Section 6.2 + Section 9 (References) |
| Institutional resources | Section 6.4 (Institutional Resources) |
| Cost-share details | Section 7.4 (Budget Narrative) |
| Diversity statement | Section 8.4 (Diversity & Inclusion) |
| Collaborators | Section 6.3 (Team) + Appendix D |

## Time Estimate

| Task | Estimated Time |
|------|----------------|
| Fill this template | 30-60 minutes |
| Write PI qualifications section | 1-2 hours |
| Request & receive support letter | 1-2 weeks (calendar time) |
| Prepare CV (if not already available) | 1-2 hours |
| Add your publications to references | 30 minutes |
| Final review of completed proposal | 1-2 hours |
| **TOTAL EFFORT** | **5-8 hours** |
| **TOTAL CALENDAR TIME** | **1-2 weeks** |

## Submission

**Final step:** Upload to NVIDIA portal

**URL:** https://academicgrants.nvidia.com/academicgrantprogram/s/Application

**You'll need:**
1. Main proposal PDF (export from Markdown or Word)
2. PI CV PDF (2 pages)
3. Letter of support PDF (signed)
4. Budget spreadsheet (if requested)

**Deadline:** December 31, 2025 (11:59 PM your local time, typically)

---

**Pro Tip:** Fill this out FIRST, then do find-replace in the main proposal document. This saves time and reduces errors!

**Need Help?** Contact your institutional Sponsored Programs Office - they review grant applications for free and provide feedback!
