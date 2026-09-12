# CUHK DAA Submission Path — Step-by-Step

**Date:** 2026-09-12
**Goal:** Apply for the CUHK Circulating Nucleic Acids Research Group DAC (EGAC00001000078) to access the Jiang PNAS 2015 + 2018 cfDNA WGS datasets.

## TL;DR

The DAA submission is a **manual, browser-only process** that requires:
1. EGA account creation (2-day validation)
2. Adding a new "institution" since the user has no affiliation
3. Filling the DAA form (8-10 fields)
4. Submitting to DAC EGAC00001000078 for Jiang Peiyong's review

**The DAC owner is `Jiang Peiyong <jiangpeiyong@cuhk.edu.hk>`** — typical academic DAC turnaround is 2-3 weeks.

---

## DAC Overview

- **DAC:** `EGAC00001000078` — The Chinese University of Hong Kong (CUHK) Circulating Nucleic Acids Research Group (CNARG)
- **Contact:** Jiang Peiyong, `jiangpeiyong@cuhk.edu.hk`
- **URL:** https://ega-archive.org/dacs/EGAC00001000078
- **Datasets covered:** 29 (verified 2026-09-12). The two we want:
  - `EGAD00001001275` — Jiang PNAS 2015 (90 HCC + 67 HBV + 36 cirrhosis + 32 healthy = 225 samples; **121 are the exact match for the 627 cohort's Jiang portion**)
  - `EGAD00001004561` — Jiang PNAS 2018 (169 plasma cfDNA samples; **strict superset of the 2015 cohort**)
- **Policy:** "Policy for Plasma DNA data sharing" — non-commercial, no redistribution, requires DAA signed by both user and institution

## Step-by-step procedure

### Step 1: Create EGA account (5 min + 2 day validation)

Go to https://ega-archive.org/register and fill the form:

| Field | Value |
|---|---|
| Title | (select) |
| First name | Yu Ching |
| Surname | Lam |
| Email | rollroyces@users.noreply.github.com (or personal email) |
| Password | (your choice) |
| **Organisation** | **MUST add new** — "Independent Researcher" (see workaround below) |
| Country | Hong Kong SAR |

**Workaround for "no affiliation":** Click the checkbox "Can't find your organisation? Add it here". Fill the new-institution form:
- Institution Name: `Independent Researcher — Yu Ching Lam`
- Institution Type: `Other` (or `Research Institute`)
- Institution Country: Hong Kong SAR
- Institution URL: `https://orcid.org/0009-0008-9113-769X` (the ORCID profile serves as a valid URL for an independent researcher)

**Submission:** Submit the form. EGA Helpdesk will email you within 2 working days to validate the account.

### Step 2: Submit DAA request (15 min)

Once your account is validated:

1. Login at https://ega-archive.org/login
2. Navigate to https://ega-archive.org/dacs/EGAC00001000078/requests
3. Click "Request Access"
4. Fill the DAA form. Typical fields:
   - **Datasets:** Check `EGAD00001001275` and `EGAD00001004561`
   - **Research purpose:** "Open-source benchmark of methylation-channel features for cfDNA cancer detection (DeepCatch project)"
   - **Non-commercial use affirmation:** Yes
   - **Data redistribution:** No
   - **Publication plan:** bioRxiv preprint → methods journal
   - **Investigator info:** Name, ORCID, country, contact
   - **Institution info:** The new "Independent Researcher" entry from Step 1
5. Sign the DAA electronically
6. Submit

### Step 3: Wait for DAC approval (2-3 weeks)

CUHK CNARG typically responds within 2-3 weeks. Watch for emails from `jiangpeiyong@cuhk.edu.hk` or the EGA Helpdesk.

### Step 4: Download + process (6 hours compute)

Once approved:
- Use `pyega3` to download the 121-sample FASTQ subset (~6 hours at 100 MB/s)
- Align with BWA-MEM (~6 hours)
- Generate fragment BEDs with FinaleToolkit (~30 min)
- Compute cleavage-ratio features (~2 hours)
- Train LR baseline on the 121 Jiang samples
- Report per-seed AUC vs the existing methylation-proxy result (AUC 0.777)

## Realistic timeline

| Step | Duration | Calendar day |
|---|---|---|
| EGA account creation | 5 min | Day 0 |
| EGA account validation | 2 working days | Day 2 |
| DAA submission | 15 min | Day 2 |
| DAA review by Jiang Peiyong | 2-3 weeks | Day 16-23 |
| FASTQ download + processing | 6 hours | Day 23-24 |
| **Total** | **~3-4 weeks** | |

## Honest assessment

**Why this path is realistic:**
- Same DAC owner as the 2015 paper = higher approval likelihood for related work
- Non-commercial, single-PI research, no conflict with original paper
- The 2018 paper superset means more samples than the 121 originally requested
- EGA system is well-established; turnaround is the only variable

**Why this path might fail:**
- The "no institutional affiliation" requirement may be a hard blocker — some DACs reject applications from unaffiliated researchers (need PI letter from an institution)
- The user previously documented they have no affiliation; EGA may require one
- Even if approved, the DAA may impose additional restrictions (e.g., "must collaborate with PI at host institution")
- 2-3 weeks is a "typical" turnaround; could be longer

**Workaround if "no institution" is a hard blocker:**
1. Email `jiangpeiyong@cuhk.edu.hk` directly to explain the unaffiliated-researcher situation and ask if they'd consider the application
2. Or: contact FinaleDB authors (Yaping Liu, Ravi Bandaru) and ask for a Globus-authenticated download of the same 121 samples (alternative path)
3. Or: use the Zenodo CRAG cohort (different patients, but available today with no auth)

## Files referenced

| URL | Verified |
|---|---|
| https://ega-archive.org/dacs/EGAC00001000078 | 2026-09-12 |
| https://ega-archive.org/datasets/EGAD00001001275 | 2026-09-12 |
| https://ega-archive.org/datasets/EGAD00001004561 | 2026-09-12 |
| https://ega-archive.org/register | 2026-09-12 |

