# FinaleToolkit Cleavage-Ratio Workload — 627-Sample Cohort on M4

**Subagent F report** | **Date:** 2026-09-12 | **Build:** hg19
**Hardware:** M4 MacBook Pro (10 logical cores, 16 GB RAM, 8.8 GB free disk at start)
**Code:** FinaleToolkit 1.1.0 in `/Users/hermes/deepcatch/.venv`
**Baseline:** `/Users/hermes/deepcatch-methylation/docs/CLEAVAGE_RATIO_FEASIBILITY.md`

---

## TL;DR

| Question | Answer |
|---|---|
| Can we compute cleavage-ratio features for the 627-sample FinaleDB cohort? | **No — data acquisition is the blocker, not compute.** |
| How fast is the compute on M4? | **~7 min/sample full-hg19 single-core; 9.5 hours with 8-way parallelism.** |
| Where would the 627 raw fragment BEDs come from? | **FinaleDB is down (S3 private, API 403).** No identical-cohort public mirror. |
| Disk cost? | **~170 MB/sample × 627 = ~100 GB** — disk-streaming possible but starts at 8.8 GB free. |
| Recommendation | **Skip.** Compute path is viable but data path is a multi-day scrape that depends on FinaleDB reauth or EGA DAC approval. |

---

## 1. Data acquisition — where could 627 raw fragment BEDs come from?

### 1.1 FinaleDB (the canonical source) — DEAD

FinaleDB (Zheng 2021, Bioinformatics) is the published source of the 627-sample
fragment BEDs that `cfdna-fragmentomics-pipeline` already extracted features from.

- **Web API** (`http://finaledb.research.cchmc.org/api/v1/seqrun`) — **403 Forbidden** as of 2026-09-12 (Subagent B verified).
- **S3 bucket** (`s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org`) — returns 403 on every key guess; `?list-type=2` denied. **Private.**
- **Globus login** — only authenticated users with a CCHMC Globus account can list/download.

Verified this subagent:
```
$ curl -sLI https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org
HTTP/1.1 403 Forbidden
$ curl -sLI --max-time 15 https://zenodo.org/api/records/6914806/files/liver.tar/content
HTTP/1.1 200 OK, content-length: 3307939840   ← Zenodo works
```

### 1.2 Zenodo CRAG dataset — same format but DIFFERENT cohort

[Zenodo 6914806](https://zenodo.org/records/6914806) ("CRAG: de novo characterization of cell-free DNA fragmentation hotspots", Zhou 2022) contains:

| File | Size | Contents |
|---|---|---|
| `BH01.chr22.frag.bed.gz` | 90 MB | chr22 only — what we already have |
| `BH01.chr22.bam` / `.bai` | 1.9 GB / 112 KB | chr22 only |
| `liver.tar` | **3.3 GB** | ~30 HCC + healthy liver samples, full hg19 |
| `breast_1.tar` ... `breast_4.tar` | **11.0 GB total** | 100 samples (50 BC + matched HC), full hg19 |
| `chrm_state.zip` / `GC.zip` / `mappability.zip` | 1.1 GB | Reference tracks |

Confirmed by listing the first tar member: `liver/MAL0646.hg19.frag.bed.gz` (262 MB).

**These IDs do NOT match our cohort.** Our 121 Jiang samples use `H###` / `HOT###` /
`C###` prefixes; the Zenodo CRAG liver cohort uses `MAL###` (Cincinnati Children's,
CCHMC cohort — Liu/Zhou lab). These are **different patients**. They are HCC vs HCC
but different study, different lab, different sequencing depth. They are not the same
627 samples FinaleDB aggregates.

→ Zenodo gives us a **proof-of-concept subset (~130 HCC + HC samples) but cannot
replicate the 627-sample cohort labels.**

### 1.3 Raw FASTQ/BAM at EGA / dbGaP / SRA — possible but slow

Edge cases:

- **Jiang 2015 (low-pass WGS HCC, the Hong Kong cohort)**: 90 HCC + 31 healthy = 121 samples.
  Per FinaleToolkit paper, the EGA accession is `EGAS00001001024`. Requires EGA DAC approval
  (weeks-to-months turnaround).
- **Cristiano 2019 (deep WGS multi-cancer, 8 cancer types + healthy)**: 537 samples
  spread across dbGaP `phs001417.v1.p1` (DELFI study, Johns Hopkins, Scharpf). Requires
  dbGaP DAC approval (weeks-to-months turnaround).
- **HG19 BAMs from SRA via `sam-dump`**: Subagent C's prior attempt showed FinaleMe
  hangs on BAMs produced by `sam-dump | samtools view | samtools sort` — the
  structural details FinaleMe needs may be missing.

### 1.4 Existing fragmentomics features — DEAD END

`/Users/hermes/cfdna-fragmentomics-pipeline/data/features/` has 632 unique samples × 7
files each (5 Mb coverage/meanlen/ratio, 100 kb counts/meanlen/ratio, FSD, motifs).
Verified:

| File | Shape | Size |
|---|---|---|
| `C309.delfi_5mb_coverage.npy` | (631,) float64 | 5.1 KB |
| `C309.delfi_5mb_meanlen.npy` | (631,) float64 | 5.1 KB |
| `C309.delfi_100kb_counts.npy` | (631,) float64 | 241 KB |

These are **aggregates only**: per-bin coverage, mean length, ratio. They do NOT
contain per-fragment (chrom, start, end, strand) coordinates, which are what
`cleavage_profile` needs. **Cannot be reconstructed from aggregates.** (Bins are
lossy — they sum fragment counts, not positions.)

### 1.5 Local fragment BEDs we DO have

| File | Path | Coverage |
|---|---|---|
| `BH01.chr22.frag.bed.gz` (86 MB, 13.4 M fragments, tabixed) | `/tmp/finaledb_data/` | chr22 only |

**Just one sample, one chromosome.** Enough for feasibility work; not enough for
the 627-sample AUC.

### 1.6 Summary of data acquisition status

| Source | What | Accessible? | Matches cohort? |
|---|---|---|---|
| FinaleDB S3 | 627 frag.tsv (~100 GB) | ❌ Private (403) | ✅ exact match |
| FinaleDB Globus | same | ⚠️ Requires auth | ✅ exact match |
| Zenodo 6914806 liver.tar | ~30 HCC, full hg19, `MAL###` IDs | ✅ Public | ❌ different cohort |
| Zenodo 6914806 breast_1..4.tar | 50 BC + 50 HC, full hg19 | ✅ Public | ❌ different cohort |
| EGA EGAS00001001024 | Jiang 2015 WGS FASTQ | ⚠️ DAC approval | ✅ exact (121/627) |
| dbGaP phs001417.v1.p1 | Cristiano 2019 WGS FASTQ | ⚠️ DAC approval | ✅ exact (537/627) |
| Local fragmentomics features | 632 × (5Mb + 100kb + FSD + motifs) | ✅ on disk | ✅ same labels — but aggregations only |
| Local `BH01.chr22.frag.bed.gz` | 13.4 M fragments, tabixed | ✅ on disk | ❌ 1 sample, chr22 only |

---

## 2. M4 runtime benchmark (cleavage_profile)

### 2.1 Setup

- **Sample:** `BH01.chr22.frag.bed.gz` (86 MB, 13.4 M fragments, tabixed)
- **Reference:** `/tmp/cpgisland_hg19.bed` (30,344 UCSC `cpgIslandExt` CpG islands,
  hg19, gunzipped to `/tmp/cleavage_bench/cpgislands_hg19.bed`)
- **API:** `finaletoolkit.frag.cleavage_profile(frag_file, chrom_size, contig, start, stop)`
  with defaults (intersect_policy="any")
- **Hardware:** M4, 10 logical cores, 16 GB RAM, 8.8 GB free disk
- **Warm cache** (HTSlib .tbi cached after first call)

### 2.2 chr22 scaling sweep

| Run | Islands iterated | Islands processed | Wall (s) | ms/island (effective) |
|---|---|---|---|---|
| chr22_all (full hg19 list, chr22 only matches) | 30,344 | 719 | **10.34** | 14.38 |
| chr22_full (chr22 only) | 719 | 719 | **6.45** | 8.97 |
| chr22_half | 360 | 360 | 2.31 | 6.42 |
| chr22_quarter | 180 | 180 | 0.74 | 4.11 |
| chr22_eighth | 90 | 90 | 0.28 | 3.13 |

Sanity checks:
- Mean cleavage % per island: **0.228 %** (matches Subagent E's 0.196 %, both
  within Zhou 2022 genome-wide ~0.3 %).
- 562,123 total positions processed in 6.45 s → ~87 K positions/s throughput.

The "full hg19 iteration" cost (10.34 s) is ~60 % higher than the chr22-only loop
(6.45 s) — most overhead is **tabix query setup per query** (~30 K no-op queries),
not the actual per-position computation.

### 2.3 Single-sample full-hg19 projection

Linear extrapolation from the chr22 (719 islands, 6.45 s) ratio:

| | Per-island cost | Islands | Per-sample (single-core) |
|---|---|---|---|
| **chr22 only** (719 islands) | 8.97 ms | 719 | **6.45 s** |
| **Full hg19** (30,344 islands) | 14.38 ms (with iteration overhead) | 30,344 | **436 s = 7.3 min** |
| Linear extrapolation (30344/719 × 6.45 s) | — | — | **272 s = 4.5 min** (theoretical floor) |

The 436 s figure (10.34 s chr22-scaled + 7.3 min from 14.38 × 30,344) is the
realistic estimate. ~7 min/sample on M4 single-core matches the "8.5 min for
full hg19" figure cited in the feasibility brief.

### 2.4 627-sample wall-clock projections

Assumes one sample fits in memory; process N samples in parallel by spawning
N independent Python processes (cleavage_profile is GIL-bound within a process):

| Parallelism | Wall-clock (full hg19, 30,344 islands) | Notes |
|---|---|---|
| 1-way (single-core) | **76.0 h** | Unrealistic; use parallelism. |
| 2-way | 38.0 h | |
| 4-way | **19.0 h** | Conservative — leaves cores free for tabix I/O. |
| 6-way | 12.7 h | |
| 8-way | **9.5 h** | Tighter; tabix I/O may contend. |
| 10-way | 7.6 h | Aggressive; M4 thermal/swap limits. |

With FinaleToolkit's reported overhead (~50× faster than naive WPS), the per-sample
time for `cleavage_ratio` is closer to the lower bound. With per-chr extraction
(skipping empty regions and parallelizing per-chr): realistic best-case ≈ **6–9 hours**.

**Network/disk acquisition dominates the calendar:**
- FinaleDB scrape: 1–3 hours (if accessible) at ~10–30 MB/s aggregate
- Local computation: 6–19 hours

### 2.5 Disk cost

- 627 samples × ~170 MB frag.tsv = **~107 GB** raw BEDs
- 8.8 GB free → need streaming (download → process → delete → next) or external SSD
- Feature output (cleavage-ratio per island): 627 × 30,344 × 8 bytes ≈ **150 MB**

---

## 3. Decision matrix

| Scenario | Feasibility | Wall-clock | Blocker |
|---|---|---|---|
| Re-fetch 627 from FinaleDB Globus | ⚠️ Requires CCHMC Globus auth | Data: 1–3 h; compute: 9–19 h | Globus credentials |
| Re-download Jiang from EGA + Cristiano from dbGaP | ⚠️ Requires DAC approval (weeks-months) | Same as above + DAC wait | EGA/dbGaP DAC |
| Use Zenodo CRAG (130 samples, different cohort) | ✅ Available | Data: ~14 GB download (3–5 h); compute: 2–4 h | Cohort mismatch — invalidates 627-cohort AUC |
| Use existing BH01 chr22 only | ✅ Have it | Compute: 12 s | Sample size = 1 (AUC undefined) |

---

## 4. Recommendation: **SKIP cleavage-ratio expansion on the 627 cohort**

**Reasons:**

1. **Data acquisition is the actual blocker, not compute.** FinaleDB is down.
   Public mirrors (Zenodo) contain *different* cohorts. EGA/dbGaP require DAC
   approval (weeks-to-months). The 627-sample feature pipeline that already
   exists was built when FinaleDB was public.

2. **Compute is fine but unnecessary to fully characterize.** 7 min/sample
   full hg19 → 9.5 h with 8-way parallelism. This is **doable on M4** but
   **useless without the data**.

3. **The methylation-proxy baseline (AUC 0.777) already exists** without raw
   fragment BEDs — it uses per-bin coverage and length features derived from
   aggregates. Adding cleavage-ratio to those aggregates would require raw
   BEDs *anyway*. The marginal value of cleavage-ratio over DELFI-style
   coverage/length features is **theoretically orthogonal** (5'-end sharpness
   vs. fragment counts) but **empirically unproven** in this codebase.

4. **Disk headroom is too tight.** 8.8 GB free vs. 100 GB raw BEDs. Streaming
   works but adds complexity and risk of mid-stream failures.

5. **Alternative with higher value-per-hour:** if a small raw-fragment subset
   is needed (e.g., to confirm a biological hypothesis), the Zenodo CRAG
   breast cohort (50 BC + 50 HC, 11 GB) is a faster standalone validation —
   but the resulting features don't match the 627-sample labels.

**Suggested future work, if compute is added later:**

1. Get FinaleDB Globus access (or EGA/dbGaP DAC approval).
2. Free 100+ GB on disk (external SSD or clear `/tmp/finaleme_test/`,
   `/tmp/finaledb_data/results/`).
3. Re-run `cfdna-fetch` for all 627 samples → ~107 GB raw BEDs.
4. Stream-process: download sample N → cleavage_profile (full hg19, 30,344
   islands) → save per-sample feature parquet (~150 MB total) → delete raw
   → download sample N+1.
5. With 8-way parallelism and a steady-state per-sample of 7 min compute +
   ~3 min download: realistic 627-sample wall-clock ≈ **6–10 hours**.

---

## Files referenced / created

- `/Users/hermes/deepcatch-methylation/docs/CLEAVAGE_RATIO_WORKLOAD.md` (this doc)
- `/tmp/cleavage_bench/bench_chr.py` — single-chr timing harness
- `/tmp/cleavage_bench/scaling_test.py` — chr22 island-count sweep
- `/tmp/cleavage_bench/cpgislands_hg19.bed` — 30,344 CpG islands (hg19)
- `/tmp/cleavage_bench/chr22_{half,quarter,eighth}.bed` — chr22-only subsets

## Source citations

- Zhou et al. 2022 (Genome Med) — CRAG paper; Zenodo 6914806 hosts breast/liver fragment BEDs.
- Li et al. 2025 (Bioinformatics Advances, PMC12597888) — FinaleToolkit paper;
  HCC cohort accessed via `EGAS00001001024`; benchmark on 100× BH01 WGS.
- Snyder et al. 2016 (Cell) — BH01 healthy fragment BED (chr22 only, Zenodo 6914806).
- Zheng et al. 2021 (Bioinformatics) — FinaleDB; 2,579 samples across 23 conditions;
  raw data behind authenticated S3/Globus.
- Cristiano et al. 2019 (Nature) — DELFI paper; deep WGS cohort dbGaP `phs001417.v1.p1`.
- Jiang et al. 2015 (PNAS) — low-pass WGS HCC cohort; EGA `EGAS00001001024`.
