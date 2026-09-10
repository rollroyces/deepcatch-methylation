# Fragmentomics vs Methylation — Head-to-Head Status

**Date:** 2026-09-10
**Task:** Option C — establish a fragmentomics baseline on TCGA tissue data for a direct head-to-head against the existing TCGA-LIHC methylation baseline (Phase 0 AUC 0.972).

---

## Bottom line

**A direct head-to-head at the sample level is NOT possible with current local data.** The two pipelines were developed on disjoint datasets with no shared samples:

| Modality | Local cohort | n | Sample type | Assay |
|---|---|---|---|---|
| Fragmentomics | FinaleDB (Cristiano 2019 + Jiang 2015) | 627 | Plasma cfDNA | Low-pass WGS |
| Methylation (Phase 0) | TCGA-LIHC | 18 | Bulk tissue | Illumina 450K array |

The cross-cohort comparison that IS available is informative but indirect — see "Proxy comparison" below.

---

## Why a direct head-to-head is blocked

### 1. No TCGA-LIHC samples in the fragmentomics cache

```
$ ls data/features/ | grep -i lihc      → (empty)
$ ls data/features/ | grep -i tcga      → (empty)
$ grep -i "lihc" data/features/labels_cross_study.tsv   → (empty)
```

The fragmentomics features directory contains **5,048 files across 658 samples** from two studies only:

| Study | n_cancer | n_healthy | n_total | Sample type |
|---|---|---|---|---|
| Cristiano 2019 | 275 | 262 | 537 | Plasma cfDNA, low-pass WGS |
| Jiang 2015 | 89 | 32 | 121 | Plasma cfDNA, low-pass WGS |

No TCGA tissue samples are present, so the fragmentomics pipeline cannot be run on the 18 TCGA-LIHC samples that have methylation data.

### 2. The fragmentomics pipeline technically could run on TCGA WXS — but it's not the same signal

The fragmentomics feature set (5Mb short/long coverage, 5Mb short/long mean fragment length, 100kb ratio, nucleosome features) was designed for cfDNA WGS. TCGA tissue samples are mostly **WXS (whole exome)**, not WGS:

- WXS covers ~1-2% of the genome (exome capture, ~60 Mb of coding regions).
- 5Mb genomic bins have highly non-uniform coverage on WXS — most bins have near-zero reads.
- Fragment-length distributions from WXS differ from cfDNA WGS because:
  - Tissue DNA is from intact cells, not apoptotic fragments.
  - WXS captures short fragments (post sonication) preferentially.
  - No plasma-derived cfDNA fragmentation biology applies.

A "fragmentomics AUC on TCGA-LIHC" would measure mostly **how well the model distinguishes tumor-tissue vs normal-tissue exome coverage variability** — not cfDNA cancer signal. This is a fundamentally different biological question.

### 3. TCGA methylation vs FinaleDB cfDNA are different biological substrates

Even if the same model architecture were applied to both:

- **Methylation signal strength**: bulk tissue tumor purity is typically 60-90%; cfDNA tumor fraction is typically 0.1-5%. Tissue methylation is intrinsically a stronger cancer signal.
- **Assay platform**: 450K array probes are pre-selected for cancer-relevant CpG islands; cfDNA WGS has uniform coverage and no probe selection.
- **Feature dimensionality**: methylation uses ~175K features after filter; fragmentomics uses ~5,000 bins. Different overfitting profiles.

---

## Proxy comparison (cross-cohort, with caveats)

Even though a sample-level comparison is impossible, the existing results can be tabulated to give a high-level signal comparison. **All numbers are 5-fold pooled-OOF L2-LR; this is the only thing genuinely comparable across the two pipelines.**

| Metric | Fragmentomics | Methylation (Phase 0) |
|---|---|---|
| **Cohort** | FinaleDB cross-study (Cristiano + Jiang) | TCGA-LIHC |
| **n samples** | 627 | 18 |
| **n cancer / healthy** | 364 / 263 (≈ 58% / 42%) | 12 / 6 (67% / 33%) |
| **Sample type** | Plasma cfDNA | Bulk tissue |
| **Assay** | Low-pass WGS | Illumina 450K array |
| **Tumor fraction** | ~0.1-5% (ctDNA) | ~60-90% (tissue purity) |
| **Feature count** | ~5,000 (5Mb + 100kb bins) | 175,893 (CpG after var filter) |
| **Model** | L2-LR, C=1000 | L2-LR, C=1.0 |
| **CV protocol** | 5-fold StratifiedKFold × 5 seeds | 5-fold StratifiedKFold × 5 seeds |
| **Headline AUC** | **0.978 ± 0.001** | **0.972 ± 0.009** |
| **Per-seed AUCs (fragmentomics, C=1000 sweep)** | [0.9796, 0.9785, 0.9761, 0.9753, 0.9781] | [0.958, 0.972, 0.986, 0.972, 0.972] |
| **AUC range across model ablations** | 0.963-0.978 (5 configs tested) | single config only (Phase 0) |
| **Source** | `cfdna-fragmentomics-pipeline/results/lr_reg_sweep.json` (best_C=1000, AUC=0.9782) | `deepcatch-methylation/results/phase0_baseline.json` |

### Interpretation

- Both modalities reach **near-0.97 AUC** on their respective cohorts.
- The fragmentomics number is **more impressive** for the screening context because cfDNA has intrinsically lower signal-to-noise than bulk tissue.
- The methylation number is a **sanity check** that the pipeline correctly reads TCGA β-values; tissue methylation is known to be near-perfectly separable at array resolution.
- The std difference (0.001 vs 0.009) reflects sample-count differences (627 vs 18), not model quality.

### Honest statement of what this comparison does NOT prove

- ❌ It does NOT prove one modality is better for cancer screening.
- ❌ It does NOT prove the two are equivalent.
- ❌ It does NOT address whether they capture complementary signal (which is the actual Phase 2 fusion question).
- ❌ It does NOT validate either pipeline on the same samples.

---

## What WOULD be a true head-to-head

### Option A — FinaleMe imputation on FinaleDB cfDNA WGS (recommended)

**Approach:** Run FinaleMe on the same 627 cfDNA WGS samples in the fragmentomics pipeline. Get imputed β-values for the identical cohort. Run the methylation baseline on those imputed β-values. Now both modalities have an AUC on the same 627 samples.

| Step | Time | Notes |
|---|---|---|
| Install Java 21 + Maven | ~30 min | One-time |
| Clone/build FinaleMe | ~1 hour | MIT-licensed, ~10k LoC Java |
| Download reference panels | ~30 min | ~5-10 GB ENCODE/Roadmap |
| Run FinaleMe on 627 cfDNA WGS | ~1-2 days | Per-sample inference, multi-thread |
| Run methylation baseline on imputed β | ~10 min | Same code as Phase 0 |

**Expected methylation AUC:** ~0.80-0.90 (per FinaleMe paper, cancer-vs-healthy cfDNA auROC ≈ 0.91 — this is documented lower than tissue methylation because the imputation noise + low tumor fraction reduce signal).

**Expected fragmentomics AUC:** ~0.978 (already measured).

**Phase 2 fusion prediction:** With both modalities at ~0.85-0.98 on the same samples, fusion could plausibly reach 0.99+ — but only the fusion ablation on identical samples answers the actual scientific question.

### Option B — Published-consensus table (quick proxy)

Compile peer-reviewed auROC numbers from:
- **Cristiano 2019** (DELFI): auROC ~0.94 for multi-cancer detection on their cohort.
- **Stackhouse FinaleMe 2024**: auROC ~0.91 for cancer-vs-healthy on cfDNA WGS.
- **Moss 2018** (GSE122126): auROC 0.85-0.95 on cfDNA 450K array.
- **Our local fragmentomics**: 0.978 on pooled 627-sample FinaleDB cohort.
- **Our local methylation**: 0.972 on TCGA-LIHC tissue (Phase 0).

Quick literature pull (~1-2 hours). Lower confidence because studies use different cohorts. Treats cross-study comparability as a literature problem rather than solving it.

### Option C — Synthetic-data validation only

Generate matched synthetic fragmentomics + methylation features (plant tumor signal at varying fractions in baseline healthy distributions). Validates the methodology but does not produce a real-data number.

**Not recommended as a substitute for A or B** — synthetic numbers cannot be compared to real-data baselines.

---

## Recommended path forward

**Do Option A.** It's the only path to a true head-to-head on identical samples, and the only path that answers the fusion question (do fragmentomics + methylation capture complementary signal?). Until then, all cross-cohort comparisons are illustrative only.

If Option A is too heavy to start now, **do Option B** as a quick stopgap — get a published-consensus table into `docs/FRAGMENTOMICS_HEAD_TO_HEAD.md` and cite the per-study cohort sizes and assays.

**Do not** invest in Option C; synthetic validation does not substitute for either of the above.

---

## Files written

- `docs/FRAGMENTOMICS_HEAD_TO_HEAD.md` (this file)
- `results/fragmentomics_comparison.json` — structured comparison data
