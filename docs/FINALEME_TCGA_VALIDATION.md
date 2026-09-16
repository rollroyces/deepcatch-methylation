# FinaleMe vs TCGA-LIHC HM450 Zero-Shot Methylation Validation

**Date:** 2026-09-16
**Sample:** BH01 (Snyder 2016 healthy plasma cfDNA; chr22 only)
**Reference:** TCGA-LIHC solid tissue normal, Illumina HumanMethylation450 (SeSAMe β)
**Script:** `scripts/finaleme_validation/run_finaleme_tcga_validation.py`
**Output:** `results/finaleme_tcga_validation.json`

---

## TL;DR

| Item | Value |
|------|-------|
| FinaleMe chr22 CpGs scored (healthy + cancer) | **489,370** each |
| Healthy-model mean β | **71.97 %** |
| Cancer-model mean β | **69.31 %** (Δ = **−2.66 pp**) |
| Global KS statistic (healthy vs cancer) | **0.115** (p ≈ 0) |
| HM450 chr22 probes in Illumina manifest | **8,552** (8,502 hg19, 50 hg18) |
| HM450 chr22 probes liftOver-validated to hg19 | **8,551** (1 hg18 liftover failure) |
| TCGA-LIHC normals loaded | **6** |
| CpG-triple overlap (TCGA + FinaleMe healthy + cancer) | **7,883** (1,036 had NA TCGA β) |
| Spearman ρ (FinaleMe healthy vs TCGA mean β) | **0.808** (n = 6,847) |
| Spearman ρ (FinaleMe cancer vs TCGA mean β) | **0.806** (n = 6,847) |
| Pearson r (FinaleMe healthy vs TCGA mean β) | **0.821** |

**Headline:** the FinaleMe pretrained HMM predictions agree strongly in **rank order** with TCGA-LIHC normal liver tissue β-values across 6,847 chr22 CpGs that overlap the Illumina HM450 array. The agreement is highest inside CpG-island shores (ρ ≈ 0.81) and inside CpG islands themselves (ρ ≈ 0.71), where the canonical "CpG islands stay unmethylated" pattern dominates. FinaleMe does **not** match TCGA exactly — the β scale differs by a few percentage points in mean (Δβ ≈ −5 pp FinaleMe healthy − TCGA, ≈ −7 pp cancer − TCGA), but the global chromatin pattern is preserved.

This is the first time FinaleMe pretrained HMM predictions have been benchmarked against a public tissue methylation array. **The biological-direction checks pass (CpG-island hypomethylation, global cancer hypomethylation, KS significance in every region).** A rigorous plasma-vs-tissue study needs paired samples and is out of scope here.

---

## 1. What was measured

### 1.1 Per-region β-value distribution on BH01 chr22 (FinaleMe only)

UCSC hg19 annotations:
- **CpG islands:** 719 chr22 intervals, 562 kb total
- **refGene genes:** 1,669 chr22 transcripts (TSS ± 2 kb promoter windows + gene-body windows)
- **Liftover chain:** `hg18ToHg19.over.chain.gz` (UCSC)

Region labels (priority order: Island > Shore > Shelf > Promoter > GeneBody > Intergenic):

| Region | n CpGs | healthy mean β (%) | cancer mean β (%) | Δβ (pp) | KS stat | KS p |
|--------|--------|---------------------|--------------------|---------|---------|------|
| Overall | 489,370 | 71.97 | 69.31 | −2.66 | 0.115 | < 1e-300 |
| Island | 50,682 | 19.81 | 19.65 | −0.16 | 0.030 | 5.98e-20 |
| N_Shore | 16,070 | 62.54 | 59.11 | −3.43 | 0.114 | 1.27e-91 |
| S_Shore | 16,571 | 61.57 | 60.79 | −0.78 | 0.102 | 5.09e-76 |
| N_Shelf | 19,765 | 77.82 | 75.13 | −2.69 | 0.127 | 1.70e-139 |
| S_Shelf | 21,896 | 78.37 | 76.29 | −2.08 | 0.122 | 7.35e-143 |
| Promoter (TSS±2 kb, non-CGI) | 16,068 | 68.52 | 65.72 | −2.80 | 0.124 | 7.32e-108 |
| GeneBody (non-CGI, non-promoter) | 220,039 | 81.21 | 78.41 | −2.80 | 0.128 | < 1e-300 |
| Intergenic | 128,279 | 77.69 | 74.05 | −3.64 | 0.125 | < 1e-300 |

**Biological reading:**

- **CpG islands stay unmethylated** in both models: median β = 0 % in islands (well-known biology, validated).
- **Global cancer hypomethylation** visible everywhere except the islands themselves: Δβ ranges −0.8 to −3.6 pp across non-island regions, consistent with the cancer methylome literature (Ehrlich 2002, Hansen et al. 2011).
- **Islands show the smallest KS shift** (0.030 vs ≈ 0.12 elsewhere), exactly as expected — islands are the most protected from cancer-associated hypomethylation.
- **Median healthy β = 91.3 %** reflects that FinaleMe's posterior emits binary-like scores (methy_count/total_count × 100); the median is pulled toward 100 by fully-methylated CpGs and the mean is pulled down by partially-methylated ones.

### 1.2 Per-region Spearman correlation with TCGA-LIHC HM450 normals

| Region | n CpGs (complete data) | ρ healthy vs TCGA | ρ cancer vs TCGA |
|--------|-------------------------|--------------------|--------------------|
| Island | 3,009 | **0.711** | **0.696** |
| N_Shore | 607 | **0.809** | **0.802** |
| S_Shore | 664 | **0.764** | **0.743** |
| N_Shelf | 491 | 0.482 | 0.497 |
| S_Shelf | 567 | 0.397 | 0.417 |
| Promoter | 748 | 0.471 | 0.476 |
| GeneBody | 615 | 0.359 | 0.385 |
| Intergenic | 146 | 0.375 | 0.364 |
| **Overall (all 6,847)** | — | **0.808** | **0.806** |

**Reading:**

- The strongest rank-order agreement is in **island shores** (ρ ≈ 0.81), the regions with the most dynamic β range (intermediate methylation). CpG-island-shore probes are exactly where DNA methylation is most biologically informative (illuminaHM450 design biased towards them).
- **Islands** also correlate strongly (ρ ≈ 0.71) — both methods agree that island CpGs are unmethylated in healthy tissue.
- **Shelves / gene-body / intergenic** have lower ρ because these regions have a narrow high-β band where ranks compress — Spearman is rank-based and a probe at β = 0.85 vs 0.95 has high measurement noise but the same biological "methylated" state.

### 1.3 Δβ (FinaleMe − TCGA) summary

| Statistic | FinaleMe healthy − TCGA (pp) | FinaleMe cancer − TCGA (pp) |
|-----------|-------------------------------|-------------------------------|
| Mean | −5.03 | −6.86 |
| Median | −2.82 | −3.24 |

Both models predict **slightly lower β** than TCGA-LIHC normals report. That is consistent with **two known biological effects**:

1. **Plasma cfDNA dilution** — FinaleMe is scoring fragments from plasma, where tumour-derived and other-tissue DNA may dilute the hepatocyte signature; cfDNA β is shifted relative to tissue β.
2. **Cell-type mixture** — TCGA "solid tissue normal" liver is mostly hepatocytes; plasma cfDNA from a healthy donor is dominated by peripheral blood mononuclear cells (PBMCs), which have a different methylation baseline at many CpGs.

The cancer model is shifted further (more hypomethylation), matching the cancer-vs-tissue direction even though the BH01 sample is healthy plasma.

---

## 2. What was NOT measured (and what would be needed)

1. **No paired plasma ↔ tissue validation.** BH01 is published Snyder 2016 cfDNA from a healthy donor; no matched tissue is available. The 0.808 ρ is therefore **a cross-tissue rank-correlation, not a plasma methylation benchmark.**
2. **No technical replicates.** FinaleMe chr22 was decoded once (n = 1). Sample-level correlations are noisy.
3. **No multichromosome coverage.** This study is chr22-only. FinaleMe full-genome Step 3 OOM is a known blocker (see `docs/FINALEME_PRETRAINED_INCOMPATIBILITY.md`); extending to all autosomes is out of scope.
4. **No matched-cancer comparator.** The TCGA-LIHC "normal" β is from non-tumour liver tissue adjacent to HCC tumours. FinaleMe's cancer model was trained on ENCODE/Blueprint cancer WGBS, which is much broader. A head-to-head cancer-vs-cancer comparison (FinaleMe cancer vs TCGA-LIHC tumour) would be the natural next step and would directly test the cancer-vs-healthy Δβ signal at the per-CpG level.
5. **No probe-level masking** of cross-reactive probes (the HM450 manifest flags ~5 % of probes as unreliable). The current triple-overlap (7,883) does not filter on Illumina's `Random_Loci` or `Methyl27_Loci` flags.
6. **No cell-type deconvolution** — TCGA "solid tissue normal" liver contains hepatocytes + immune infiltrate; plasma cfDNA contains PBMC + liver + other tissue signatures. Without deconvolution the comparison cannot reach quantitative agreement, only qualitative rank-order agreement.

**To make this validation rigorous** you would need:
- FinaleMe decoded β for ≥ 20 paired plasma samples (case + control, ideally from a published HCC cfDNA cohort).
- A second methylation platform on the same plasma (HM450 or EPICv2) for direct per-probe correlation.
- Cell-type deconvolution of the tissue methylation truth (e.g. reference-based: Hepatic / PBMC / Immune signatures).
- Multi-chromosome coverage (not chr22-only) — needs the Step 3 OOM fix.

---

## 3. Methods (honest)

### 3.1 Inputs

- **FinaleMe output (BH01 chr22):**
  `results/finaleme_decode/BH01_chr22_{healthy,cancer}_decoded.bed.gz`
  Schema: `#chr / start / end / methy_perc_predict / methy_count_predict / total_count_predict / methy_perc_obs / methy_count_obs / total_count_obs`. We use `methy_perc_predict` as the predicted β (in %).
- **UCSC CpG islands:** `data/reference/cpgIslandExt.txt.gz` (full hg19 table → filtered to chr22 = 719 intervals). Sub-annotation: N_Shore / S_Shore = 0–2 kb from island; N_Shelf / S_Shelf = 2–4 kb from island (Illumina HM450 convention).
- **UCSC refGene:** `data/reference/refGene.txt.gz` → 1,669 chr22 transcripts. Promoter = TSS ± 2 kb; GeneBody = txStart..end (minus promoter / CGI labels).
- **HM450 manifest:** Illumina v1.2 CSV (192 MB, downloaded once and cached at `/tmp/HM450_manifest.csv`). 8,552 chr22 probes (8,502 in build 37 = hg19; 50 in build 36 = hg18). LiftOver hg18 → hg19 via `pyliftover` chain for the 50 build-36 probes; 1 liftover failure.
- **TCGA-LIHC HM450 normals:** 6 SeSAMe-processed β-value TSVs at `data/raw/tcga_lihc_subset/{fid}.tsv` (486,427 probes each; 1,036 chr22 probes had NA in ≥ 1 sample, dropped from correlation).

### 3.2 Coordinate harmonisation (the bug I had to fix)

The Illumina MAPINFO column is **1-based**. FinaleMe BED `start` is **0-based**. A naive match of MAPINFO = FinaleMe `start` gave only 99 of 8,552 matches. After subtracting 1 from MAPINFO (MAPINFO − 1 = 0-based), 7,883 matches. Of the 50 build-36 (hg18) probes, 1 failed liftover; the remaining 49 were lifted and matched directly.

### 3.3 Statistical tests

- **Per-region β distribution:** mean, median, std, 25th/75th percentile (pandas `.describe()`).
- **Healthy-vs-cancer shift:** two-sample Kolmogorov–Smirnov test (`scipy.stats.ks_2samp`).
- **FinaleMe vs TCGA-LIHC:** Spearman rank correlation (`scipy.stats.spearmanr`) — primary metric. Pearson r reported for completeness. NaN β dropped per probe (6,847 of 7,883 complete).

### 3.4 Runtime

Wall-clock for the full pipeline (re-running end-to-end on M4): **≈ 4 s** (no FinaleMe decoding, only the annotation + comparison step).

---

## 4. Files

| Path | Purpose |
|------|---------|
| `scripts/finaleme_validation/run_finaleme_tcga_validation.py` | Main script — load FinaleMe, annotate, KS test, correlation. |
| `results/finaleme_tcga_validation.json` | Machine-readable results (schema_version 1.0). |
| `data/reference/cpgIslandExt.txt.gz` | UCSC hg19 CpG islands (downloaded once). |
| `data/reference/cpg_islands_chr22.bed` | chr22-only island BED6 (derived). |
| `data/reference/refGene.txt.gz` | UCSC hg19 refGene (downloaded once). |
| `data/reference/hg18ToHg19.over.chain.gz` | LiftOver chain for 50 hg18 probes (UCSC). |
| `/tmp/HM450_manifest.csv` | Illumina HM450 v1.2 manifest (cached, 192 MB). |
| `test/test_finaleme_tcga_validation.py` | pytest guard rails. |
| `docs/FINALEME_TCGA_VALIDATION.md` | This document. |

---

## 5. Honest bottom line

- **Per-region BH01 chr22 statistics are real and biologically sensible:** islands unmethylated, non-island regions ~70–80 % methylated, cancer uniformly hypomethylated by ~2.7 pp (less in islands — exactly as the literature predicts).
- **Spearman ρ ≈ 0.81 between FinaleMe healthy and TCGA-LIHC normals** is strong rank-order agreement at the per-CpG level for n = 6,847 chr22 probes — this is a real cross-platform signal, not noise.
- **This is not a paired-plasma validation.** The right next step is decoding FinaleMe on ≥ 1 cfDNA cohort with matched-tissue or matched-HM450 truth; until then, the Spearman ρ is a sanity check, not a regulatory-grade metric.
- **Step 3 OOM (full-genome FinaleMe) is the actual blocker** for a multi-chromosome version of this validation; chr22 alone is too narrow to make clinical claims.