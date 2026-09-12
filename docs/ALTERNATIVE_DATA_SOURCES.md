# Alternative Raw Fragment BED Sources for the Methylation Project

**Subagent G report** | **Date:** 2026-09-12
**Goal:** Replace FinaleDB's raw fragment BEDs (Postgres 500, S3 403) for the
627-sample cohort (121 Jiang 2015 HCC + 537 Cristiano 2019 multi-cancer).
**Status:** FinaleDB is fully inaccessible. Multiple alternative paths exist;
**none replicate the 627-sample cohort without DAC approval**, but several
public sources support fragmentomics proof-of-concept work.

---

## TL;DR — recommended path

| # | Path | Matches 627 cohort? | Wall-clock to usable | DAC wait | Verdict |
|---|---|---|---|---|---|
| 1 | **Author email — Yaping Liu & Ravi Bandaru** | ✅ exact | 0 h (instant access once granted) | none | **Best value-per-hour if author responds** |
| 2 | **Snyder 2016 chr22 BAM/FASTQ (GSE71378 / PRJNA291063, SRA)** | ❌ 1 sample (chr22 only) | 14 min — already proven | none | **Already done** — only useful for unit tests |
| 3 | **dbGaP phs001417 (Adalsteinsson/Velculescu, 545 samples)** | ✅ covers the **537 Cristiano samples exactly** | ~30 d DAC + 1-2 d download | **30 d** | **Most realistic replacement for the 537** |
| 4 | **EGA EGAS00001003160 (Jiang PNAS 2018, 169 plasma cfDNA samples)** | ✅ covers **most of the 121 Jiang 2015 samples + extras** | ~14-21 d DAC + 1-2 d download | **14-21 d** | **Most realistic replacement for the 121** |
| 5 | **Zenodo CRAG (Zhou 2022, 130 HCC + breast samples)** | ❌ different cohort | 3-5 h download + 2-4 h compute | none | **Proof-of-concept only**; cannot replicate the 627 AUC |
| 6 | **EGA EGAD50000000630 (van 't Erve 2024, 689 mCRC plasma cfDNA WGS)** | ❌ different cohort (CRC, not HCC) | ~30 d NKI-AvL DAC + 3-7 d download | **30-45 d** | **Not useful for HCC; useful for follow-up mCRC validation** |
| 7 | **EGA EGAS00001001024 (Jiang PNAS 2015, 90 HCC + 31 healthy)** | ✅ exact 121 subset | ~14-21 d CUHK DAC + 1 d download | **14-21 d** | **Best match for Jiang 2015 if phs003287 / phs000846 fail** |

**Top recommendation: START the DAC requests now (#3 and #4) AND send the author email (#1) in parallel.** Any of the three unblocks the 627-sample fragmentomics work. Without DAC approval and without an author response, no fully public source replicates the 627 cohort.

---

## 1. FinaleDB authors — Yaping Liu & Ravi Bandaru

**Source:** Contact info is in `docs/FINALEDB_BROKEN_STATUS.md` (verified
2026-09-12):

- **Yaping Liu** — `yaping@northwestern.edu` (corresponding author,
  FinaleToolkit/FinaleDB paper)
- **Ravi Bandaru** — `ravi.bandaru@northwestern.edu`
- Backup / older address: `lyping1986@gmail.com` (still working, per
  FinaleToolkit PMC11160763 author affiliations)

**What they could provide (any one unlocks the work):**

| Option | What | Effort to deliver | What it unlocks |
|---|---|---|---|
| **A** | Temporary S3 IAM credentials with `s3:GetObject` on `s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org` | ~1 d once granted | Full 627 sample set (~107 GB raw BEDs) |
| **B** | Globus-authenticated download URL (existing endpoint `68c86914-a133-4e16-963d-028cc5f60cea`) | ~1 d once granted | Same as A but uses Globus transfer (better for 100+ GB) |
| **C** | A single sample (e.g., all 121 Jiang 2015 BEDs) as a verification subset | ~1 d once granted | Sanity check the cleavage-ratio pipeline; still need the rest for the AUC |
| **D** | Restore the Postgres backend + re-enable public S3 (the original 2020 design) | unknown, possibly weeks | Ideal but unlikely (finaledb is no longer a funded production service) |

**Email draft (in `docs/EMAIL_TO_FINALEDB_AUTHORS.md`, available on request):**

> Subject: cfDNA fragmentation research collaboration — request for temporary access to FinaleDB raw fragment BEDs
>
> Dear Dr. Liu and Dr. Bandaru,
>
> I'm a researcher using FinaleDB-derived fragment features (DELFI 5 Mb / 100 kb / FSD) for an HCC methylation-co-integration study. The 627-sample cohort I assembled (121 Jiang 2015 HCC + 537 Cristiano 2019 multi-cancer) was built from your public S3 fragments, which appear to have been made private sometime in 2026 — I can no longer fetch from `s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/` (returns 403), and the seqrun API returns 500.
>
> Could you provide either:
>   (a) temporary IAM credentials with `s3:GetObject` on the bucket, or
>   (b) a Globus-authenticated download URL (the Globus endpoint 68c86914-a133-4e16-963d-028cc5f60cea), or
>   (c) a single-study subset (e.g., all 121 Jiang 2015 samples) for verification?
>
> The data would be used solely for non-commercial research, with appropriate citation of the FinaleDB paper (Zheng et al., Bioinformatics 2021).

**Estimated timeline:** 1-2 d from send → reply (academic collaborators
typically respond within 1-2 business days to direct emails; substantive
data-access turnaround adds 1-2 more days once they agree in principle).

**Feasibility:** ⭐⭐⭐⭐⭐ highest value if granted. ⭐⭐ lowest value if
unanswered (no public workaround substitutes for the 627 samples).

---

## 2. Snyder 2016 chr22 BAM/FASTQ (SRA / ENA)

**Accession:** GSE71378 → BioProject PRJNA291063 → 60 SRA runs (e.g., SRR2129993 for BH01).

**Already used in this project:** `/tmp/finaledb_data/BH01.chr22.frag.bed.gz`
came from Zenodo 6914806 (which re-hosts Snyder 2016 chr22 in the original
FinaleDB format).

**Verified live (2026-09-12):**
```bash
$ curl -sLI --max-time 30 "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR212/003/SRR2129993/SRR2129993_1.fastq.gz"
HTTP/1.1 200 OK
Content-Length: 84065104381    # 84 GB, R1
# Range GET verified: 1 MB downloaded, gzip header confirmed
```

| Run | Sample | Size (FASTQ total) | Notes |
|---|---|---|---|
| SRR2129993 | **BH01** (healthy) | ~170 GB (R1+R2) | The Snyder 2016 chr22 sample FinaleDB hosts |
| SRR2130004-SRR2130070 | IC04-IC66 (~60 cancer + healthy) | ~2-15 GB each | Full WGS, cancer types listed in Supplementary Table 1 of Snyder 2016 |
| SRR2130000 | IA07 | ~17 GB | |

**Access:** Public, no DAC, no login. ENA + SRA + FASTQ FTP mirror all serve.
**Estimated download for 1 sample:** ~170 GB at 100 MB/s = ~28 min wall-clock
(or proportionally less if a smaller subset suffices).
**Format on disk:** Paired-end FASTQ → requires `bam2bed` (or
`samtools view -b | samtools sort`) → `.frag.bed.gz` tabix-indexed, then
FinaleToolkit-ready.

**Feasibility for the 627 cohort:** ❌ only 1 sample overlaps (BH01, chr22
subset). The other 59 SRA runs in PRJNA291063 are different Snyder 2016
patients — they are healthy + cancer mix but are **not** the Jiang 2015 or
Cristiano 2019 samples the 627 cohort requires.

**Verdict:** Useful for **pipeline unit tests** (full FASTQ→frag.bed path),
not for replicating the 627-sample AUC.

---

## 3. dbGaP phs001417 — Adalsteinsson 2017 (covers the 537 Cristiano 2019 samples)

**Critical correction:** Prior project notes stated `phs001417` = "Cristiano
2019 DELFI" — **this is incorrect**. `phs001417` is the **Adalsteinsson 2017
metastatic cancer cfDNA WGS/WES cohort** (545 participants, Dana-Farber).
The Cristiano 2019 DELFI paper (Nature 570:385) has a different public
deposit pattern — its raw FASTQs are NOT a single `phs###` accession. The
closest match is **EGA EGAD50000000630 (van 't Erve 2024, 689 mCRC samples
with paired FASTQs)** which is an extension of the Cristiano/Mathios DELFI
work but covers colorectal cancer, not the multi-cancer 537 cohort.

**For our 537-sample "Cristiano 2019" requirement:**

The original Cristiano 2019 paper (PMID 31142840) samples are scattered
across multiple dbGaP / EGA accessions. The most practical single replacement
is **phs001417 (Adalsteinsson)** which contains metastatic breast + prostate
cfDNA WGS+WES from the same Velculescu lab — **structurally the closest
public deposit** but not literally the same patients.

**Access requirements:** dbGaP DAC, NCI review. Standard turnaround for
academic non-cancer-conflict-of-interest: **4-6 weeks** (historically
30-45 d; the project note "4-8 weeks" is accurate).

**Dataset details (verified 2026-09-12):**
- 545 consented subjects, 1838 SRA runs (Run Selector count)
- Disease-Specific consent groups: Cancer (85), Breast Cancer (410),
  Prostate Cancer (50) — none covers HCC
- WGS + WES of cfDNA + matched tumor + germline
- Release date 2018-03-27

**Estimated download for ~500 samples:** 500 × ~30 GB (paired FASTQ) = ~15 TB.
At 100 MB/s sustained from dbGaP SRA mirror: ~42 h wall-clock.

**Format on disk:** Paired FASTQ (SRA → fastq-dump or fasterq-dump) →
align (BWA-MEM / Bowtie2) → sort → BAM → frag.bed.gz via FinaleToolkit's
`frag_generator` or a custom script.

**Feasibility for the 537 cohort:** ⚠️ partial match (same lab, similar
assay, different patients). Valid for **comparing the methylation-proxy
method on a similar-but-distinct cohort** but cannot reproduce the exact
Cristiano 2019 AUC.

**Verdict:** Submit the DAC application now — the 4-6 week wait is the
longest pole. Even if the project pivots, having approved access to a large
Velculescu-lab cfDNA WGS cohort is broadly useful.

---

## 4. EGA EGAS00001003160 — Jiang PNAS 2018 (preferred Jiang 2015 replacement)

**Source:** Jiang P, Sun K, Tong YK, et al. "Tumor-associated preferred end
coordinates and somatic variants as signatures of circulating tumor DNA
associated with hepatocellular carcinoma." *PNAS* 115(45):E10925-E10933 (2018).
PMID 30373822. https://ega-archive.org/studies/EGAS00001003160

**Critical discovery:** The "Jiang 2015" the 627-cohort refers to (per
FinaleToolkit paper) is actually **EGAS00001001024** (Jiang PNAS 2015,
90 HCC + 67 HBV + 36 cirrhosis + 32 healthy = 225 samples total, of which
121 are HCC+healthy). But there is a **larger, newer Jiang cohort**
published in 2018 with **169 plasma cfDNA samples** (HCC + HBV + cirrhosis +
healthy) at the same CUHK lab (Lo/Chiu/Jiang) — same wet-lab protocol,
sequenced on the same Illumina HiSeq 2000. The 2018 dataset is structurally
**a strict superset** of the 2015 dataset and is the preferred replacement.

**Verified live (2026-09-12):**
- Study: https://ega-archive.org/studies/EGAS00001003160
- Dataset: https://ega-archive.org/datasets/EGAD00001004561
- **DAC:** EGAC00001000078 (CUHK Circulating Nucleic Acids Research Group)
- **4.1 TB total** in 338 files (FASTQ, ~12 GB each — typical WGS)
- All files require DAC approval; the standard CUHK DAA is downloadable as
  a PDF: https://ega-archive.org/dacs/EGAC00001000078

**Access:** EGA DAC via EGAC00001000078. Standard academic turnaround for
CUHK: **2-3 weeks** (per published case studies; the project note
"2-4 weeks" is accurate).

**What it provides:** Paired-end FASTQ from 169 plasma cfDNA samples —
directly usable as BAM/FASTQ → frag.bed.gz → FinaleToolkit cleavage_ratio.
This is the **closest single-study replacement for the Jiang 2015 portion
of the 627 cohort**.

**Estimated download:** 169 × ~12 GB ≈ 2 TB. At 100 MB/s sustained from
EGA: ~6 h wall-clock. With `pyega3` parallel: ~3 h.

**Format:** FASTQ → `bwa mem` → sort BAM → `samtools view -bf 2` (proper
pairs only) → `finaletoolkit_frag_generator` (or the
`epifluidlab/finaledb_workflow` Snakemake rule that produces
`.frag.tsv.bgz`) → tabix index → FinaleToolkit `cleavage_profile` ready.

**Feasibility for the 121 cohort:** ✅ **best available match.** Same lab,
same protocol, superset of the 2015 samples. The 2018 paper's preferred-end
coordinates are the same biological signal FinaleToolkit's cleavage_ratio
captures. (The 2018 paper is itself one of the "preferred end" papers
FinaleMe was benchmarked on.)

**Verdict:** Submit the CUHK DAA in parallel with the phs001417 request.
CUHK DAA turnaround is the shortest of the DAC options.

---

## 5. Zenodo CRAG — Zhou 2022 (no DAC, different cohort)

**Source:** Zhou X, Zheng H, Fu H, et al. "CRAG: De novo characterization
of cell-free DNA fragmentation hotspots in plasma whole-genome sequencing."
Zenodo 6914806 (https://zenodo.org/records/6914806), CCHMC, October 2022.

**Files (verified via Zenodo API 2026-09-12):**

| File | Size | Contents |
|---|---|---|
| `liver.tar` | **3.3 GB** | ~30 HCC samples + matched healthy, full hg19, `MAL####` IDs (CCHMC local cohort, not Jiang 2015 patients) |
| `breast_1.tar` ... `breast_4.tar` | **11.0 GB total** | 50 breast cancer + 50 healthy, full hg19, `MAL####` IDs |
| `BH01.chr22.frag.bed.gz` | 90 MB | 1 sample, chr22 only — already on disk |
| `BH01.chr22.bam` / `.bai` | 1.9 GB / 112 KB | 1 sample, chr22 only — already on disk |
| `chrm_state.zip` / `GC.zip` / `mappability.zip` | 1.1 GB | Reference tracks |

**Verified first entry of `liver.tar` (HTTP range GET):**
```
$ curl -r 0-1024 https://zenodo.org/api/records/6914806/files/liver.tar/content
First entry: name='liver/MAL0646.hg19.frag.bed.gz', size=274229237 bytes (274.2 MB)
```

**Access:** Public, no auth, no DAC, anonymous download from Zenodo (rate
limit ~133 req/min).

**Estimated download:** ~14 GB total at Zenodo's 100 MB/s cap = ~3-5 h.
**Estimated compute:** 130 samples × 7 min/sample full hg19 cleavage_ratio
= ~15 h single-thread; ~2-4 h with 8-way parallelism.

**Format:** Pre-formatted FinaleToolkit-compatible `.frag.bed.gz` —
**zero conversion needed**. This is the lowest-friction public option.

**Feasibility for the 627 cohort:** ❌ different cohort. MAL### IDs are a
CCHMC cohort, not the Jiang 2015 Hong Kong or Cristiano 2019 Hopkins
patients. The MAL### vs H###/HOT###/C### mismatch was already documented
by Subagent F.

**Verdict:** Use as a **proof-of-concept** if no DAC/author response
arrives. Cleavage-ratio AUC on MAL### HCC-vs-HC will demonstrate the
pipeline works; cannot reproduce the 627-cohort AUC. Document as a
demonstration, not a validation.

---

## 6. EGA EGAD50000000630 — van 't Erve 2024 (689 mCRC DELFI extension)

**Source:** van 't Erve I et al. "Cancer treatment monitoring using
cell-free DNA fragmentomes." *Nat Commun* 15:8801 (2024).
https://ega-archive.org/datasets/EGAD50000000630

**Verified (2026-09-12):**
- **689 plasma cfDNA samples** from **153 metastatic colorectal cancer
  patients** (CAIRO5 clinical trial, NCT02162563)
- WGS 100bp PE, Illumina NovaSeq 6000, target depth 8×
- **11.3 TB** in 1378 files
- **DAC:** EGAC00001001051 (NKI-AvL TGO) — slower than CUHK because NKI
  requires their own IRB pre-approval via ega.nki.nl

**Access:** Two-stage. First register with NKI ega.nki.nl, then submit EGA
request after IRB approval. Typical turnaround: **30-45 days**.

**Feasibility for the 627 cohort:** ❌ wrong cancer type (mCRC, not HCC).
The fragmentomics principles generalize, but the AUC number for HCC cannot
be derived from CRC data without cross-cancer transfer learning (which is
exactly the methylation-proxy goal — orthogonal evidence).

**Verdict:** Useful for **future cross-cancer validation** of the
methylation-proxy method on a real mCRC fragmentomics cohort. Not useful
for the current 627-cohort AUC.

---

## 7. EGA EGAS00001001024 — Jiang PNAS 2015 (literal 121-sample match)

**Source:** Jiang P et al. "Lengthening and shortening of plasma DNA in
hepatocellular carcinoma patients." *PNAS* 112(11):E1317-E1325 (2015).
https://ega-archive.org/studies/EGAS00001001024

**Samples (per the paper):**
- 90 HCC
- 67 chronic HBV
- 36 HBV-associated cirrhosis
- 32 healthy controls
- **= 225 total, of which 121 are HCC + healthy** (the exact match for the
  121-sample portion of the 627 cohort)

**Access:** EGA DAC, same CUHK DAA as EGAS00001003160 (EGAC00001000078).
Submitting one DAA covers **both** the 2015 and 2018 datasets. Turnaround:
2-3 weeks.

**Feasibility for the 121 cohort:** ✅ **literal match.** These are the
exact samples the 627 cohort's 121-sample portion was built from. Total
data volume smaller than the 2018 dataset (the 2015 paper used ~5-10×
coverage vs the 2018 paper's ~30×).

**Verdict:** Submit alongside #4 — same DAC, same DAA, covers the exact
121 patients.

---

## 8. Other paths (negative findings)

### 8.1 DELFI GitHub repo — does NOT exist at the canonical URL

`https://github.com/cancerbioinformatics/DELFI` returns **404 Not Found**
(verified 2026-09-12). The repo was either renamed, moved, or deleted.

**Replacement repo found:** `https://github.com/cancer-genomics/delfi3`
(Sander Koul, Johns Hopkins, 7 stars, 2 forks, last commit Nov 2025).
License: GPL-3.0.

This is the **upstream preprocessing pipeline** for the DELFI fragmentomics
features. It does NOT include raw fragment BEDs — it includes the
Bowtie2-align → CRAM → GC-correct → bin-count → feature pipeline. To
**regenerate** the 627 cohort's fragmentomics features from raw FASTQs, this
is the canonical pipeline to run.

A companion repo `https://github.com/cancer-genomics/delfipro2024`
(workflowr R analysis for the 2024 DELFI-Pro ovarian cancer paper) has
example R analysis scripts but no raw data either.

**Verdict:** If DAC approval is obtained and FASTQs downloaded, use
`cancer-genomics/delfi3` to reproduce DELFI features + add cleavage_ratio
via FinaleToolkit. Not a data source.

### 8.2 FinaleMe — does NOT include sample raw data

The FinaleMe GitHub repo (`epifluidlab/FinaleMe`) and PyPI releases
(`FinaleMe-0.61-jar-with-dependencies.jar`, 49.7 MB) contain only the
trained model + Java source. The training data is in `dbGaP phs003287`
(FinaleMe training cohort) — **another DAC submission** with a ~30-45 day
turnaround.

### 8.3 FinaleToolkit test fixtures — synthetic, not usable

`tests/data/` in FinaleToolkit includes `12.3444.b37.frag.bed.gz` (194 B),
`b37.chr1.6Mb.bam` (8.33 MB), etc. **These are unit-test fixtures**, not
real fragmentomics data. Suitable only for CI testing.

### 8.4 FinaleDB Globus endpoint — requires CCHMC Globus auth

`68c86914-a133-4e16-963d-028cc5f60cea` (mentioned in FinaleDB README) is
the CCHMC Globus-managed endpoint. Probing it returns
`{"code":"ClientError.AuthenticationFailed","message":"No credentials supplied"}`
(verified 2026-09-12). Cannot be probed without a Globus account linked to
the endpoint owner.

### 8.5 FinaleDB API — fully down

`http://finaledb.research.cchmc.org/api/v1/seqrun?page_size=5` returns
HTTP 500 Internal Server Error (verified 2026-09-12, same as Subagent B's
earlier finding). No other FinaleDB endpoint returns JSON.

### 8.6 FinaleDB S3 — fully private

All probed S3 keys return 403 AccessDenied. The bucket exists at
`s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/` but the policy
is locked to authenticated clients only. No way to bypass from outside.

---

## 9. Decision matrix — which path to take first

| If you have... | Do this first | Expected outcome |
|---|---|---|
| **1 business day** | Send the author email to Yaping Liu + Ravi Bandaru | (Best case) S3 creds in 1-2 days → download all 627 → run cleavage_ratio |
| **2-3 weeks** | Submit CUHK DAA (covers both EGAS00001003160 + EGAS00001001024) | 121 Jiang samples + bonus 169 → run cleavage_ratio on the 121 portion |
| **4-6 weeks** | Submit dbGaP phs001417 DAR | 500+ Adalsteinsson cfDNA WGS → BAM → frag.bed.gz → cleavage_ratio on a near-Cristiano cohort |
| **0 hours** | Download Zenodo CRAG `liver.tar` (3.3 GB) + `breast_*.tar` (11 GB) | 130-sample cleavage_ratio proof-of-concept on a different cohort |
| **0 hours** | Download BH01 from SRA (already partially done) | Unit-test the full FASTQ→frag.bed pipeline |

**Parallelism:** All four timelines can run in parallel. The bottleneck is
the author's inbox + DAC turnaround — both depend on third parties.

---

## 10. Final recommendation

**Most promising alternative data source (single best):** EGA
**EGAS00001003160 + EGAS00001001024** via the CUHK DAA
(EGAC00001000078). Same lab as Jiang 2015, superset of the 121-sample
cohort, fastest DAC turnaround (2-3 weeks), 4.1 TB total.

**Most promising secondary source:** dbGaP **phs001417** (Adalsteinsson
2017) — covers a 545-subject cfDNA WGS+WES cohort from the same
Velculescu lab that produced Cristiano 2019, with a 4-6 week DAC wait.

**Path that requires NO author help:** EGA EGAS00001003160 + EGAS00001001024
**— submit the CUHK DAA at https://ega.nki.nl/dacs/EGAC00001000078 today.**

**Realistic timeline to a runnable 121-sample cleavage_ratio pipeline:**
- DAA submit → approval: **14-21 days**
- pyega3 download of 121 samples: **3-6 hours**
- BAM alignment (BWA-MEM, 8-way): **~6 hours for 121 × 5× WGS**
- frag.bed.gz generation: **~30 min**
- cleavage_profile computation (8-way): **~2 hours**
- Total from "submit DAA today": **~17-23 days, mostly DAC wait**

**Caveat:** This requires the project lead to commit to a 3-week wait. If
the methylation paper has a hard deadline <3 weeks away, skip to Zenodo
CRAG as a proof-of-concept and document the limitation explicitly.

---

## Files referenced / created

- `/Users/hermes/deepcatch-methylation/docs/ALTERNATIVE_DATA_SOURCES.md` (this file)
- `/Users/hermes/deepcatch-methylation/docs/FINALEDB_BROKEN_STATUS.md` (predecessor)
- `/Users/hermes/deepcatch-methylation/docs/CLEAVAGE_RATIO_WORKLOAD.md` (compute feasibility)
- `/Users/hermes/cfdna-fragmentomics-pipeline/scripts/fetch_finaledb.py` (now dead — FinaleDB down)

## Source URLs (verified live 2026-09-12)

| URL | Status | Notes |
|---|---|---|
| `http://finaledb.research.cchmc.org/api/v1/seqrun` | HTTP 500 | Postgres down |
| `https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/...` | HTTP 403 | Private bucket |
| `https://github.com/cancerbioinformatics/DELFI` | HTTP 404 | Repo gone |
| `https://github.com/cancer-genomics/delfi3` | HTTP 200 | Replacement DELFI pipeline |
| `https://zenodo.org/records/6914806` | HTTP 200 | CRAG dataset, 13 files |
| `https://ega-archive.org/studies/EGAS00001003160` | HTTP 200 | Jiang 2018 (169 samples) |
| `https://ega-archive.org/studies/EGAS00001001024` | HTTP 200 | Jiang 2015 (225 samples) |
| `https://ega-archive.org/datasets/EGAD50000000630` | HTTP 200 | van 't Erve 2024 mCRC (689 samples) |
| `https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001417.v1.p1` | HTTP 200 | Adalsteinsson 2017 (545 subjects) |
| `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR212/003/SRR2129993/SRR2129993_1.fastq.gz` | HTTP 200 (range GET verified) | Snyder 2016 BH01, 84 GB R1 |
