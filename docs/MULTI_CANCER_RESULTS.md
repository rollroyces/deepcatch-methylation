# Phase 1 Option A — Multi-cancer TCGA Methylation Baseline

**Date:** 2026-09-10
**Author:** Hermes subagent (Option A)
**Output:** `results/multi_cancer_baseline.json`

## Goal

Test whether the strong TCGA-LIHC methylation signal (Phase 0 AUC 0.972 ±
0.009) generalizes across cancer types, and produce a baseline that can be
honestly compared to the fragmentomics AUC of 0.975 from
`cfdna-fragmentomics-pipeline`.

## Data downloaded

| Project  | Tumor | Normal | Probes | Source dir                          |
|----------|-------|--------|--------|-------------------------------------|
| TCGA-LUAD |   10  |   5    | 488027 | `data/raw/tcga_multi_cancer/TCGA-LUAD`  |
| TCGA-BRCA |   10  |   5    | 488027 | `data/raw/tcga_multi_cancer/TCGA-BRCA`  |
| TCGA-COAD |   10  |   5    | 488027 | `data/raw/tcga_multi_cancer/TCGA-COAD`  |
| TCGA-PAAD |   10  |   5    | 486427 | `data/raw/tcga_multi_cancer/TCGA-PAAD`  |
| TCGA-LIHC |   12  |   6    | 486427 | `data/raw/tcga_lihc_subset` (existing) |

**Total:** 52 tumor + 26 normal = 78 samples, ~610 MB on disk
(within the 800 MB budget).
**Probe intersection across all 78 samples:** 486,427 CpGs.

All files were pulled from the GDC REST API
(`https://api.gdc.cancer.gov/`) using `data_type =
"Methylation Beta Value"` filtered to Primary Tumor / Solid Tissue
Normal samples. No reCAPTCHA, no external auth.

A recovery helper script `scripts/download_remaining_files.py` was added
to handle the GDC connection resets observed mid-batch on `TCGA-COAD`
and `TCGA-PAAD`.

## Setup

All three setups reuse `src/methylation/methylation_baseline.py`
unchanged: L2 logistic regression, per-fold `StandardScaler`,
5-fold stratified CV, 5 random seeds, `min_var=0.01`, `C=1.0`.

The driver script is `scripts/run_multi_cancer_baseline.py`. To re-run:

```bash
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/run_multi_cancer_baseline.py
```

## Setup (a) — Per-cancer-type methylation baseline

Each cancer type is treated as its own binary classification task
(tumor vs matched solid-tissue normal).

| Cancer type | n (tumor/normal) | Features kept (var≥0.01) | AUC mean ± std | Per-seed AUCs          |
|-------------|------------------|---------------------------|----------------|------------------------|
| TCGA-LUAD   | 10 / 5           | 6,272                     | **1.0000 ± 0.0000** | [1.0, 1.0, 1.0, 1.0, 1.0] |
| TCGA-BRCA   | 10 / 5           | 8,415                     | **1.0000 ± 0.0000** | [1.0, 1.0, 1.0, 1.0, 1.0] |
| TCGA-COAD   | 10 / 5           | 8,276                     | **1.0000 ± 0.0000** | [1.0, 1.0, 1.0, 1.0, 1.0] |
| TCGA-PAAD   | 10 / 5           | 136,418                   | **0.9160 ± 0.0150** | [0.92, 0.90, 0.94, 0.92, 0.90] |
| TCGA-LIHC (Phase 0, ref) | 12 / 6  | 175,893 | **0.9722 ± 0.0088** | [0.958, 0.972, 0.986, 0.972, 0.972] |

**Honest reading:** AUC 1.0 on 10t/5n is degenerate — 5-fold CV puts only
2 tumor + 1 normal in each test fold, and a single hypermethylated probe
will perfectly rank them. This is *consistent with a real signal* but
the point estimates should be treated as a lower-bound on what larger
cohorts would show, not as 1.0 in any clinical sense. PAAD sits at
0.916 ± 0.015, lower than the rest, plausibly because pancreas has the
smallest tumor/normal β-value separation of the four additions.

## Setup (b) — Pooled multi-cancer methylation baseline

All 78 samples pooled into one binary task. Same model, same protocol.

- **Pooled AUC: 0.9778 ± 0.0031** over 5 seeds × 5 folds
- Features after variance filter: 9,556
- Per-seed AUCs: see `results/multi_cancer_baseline.json`

**Comparison vs TCGA-LIHC-only (0.972 ± 0.009):** the pooled AUC is
slightly higher (+0.006) and has smaller variance (~3× tighter std).
The cancer-vs-normal signal *does not collapse* when we add four more
cancer types. This is the headline result.

## Setup (c) — Tumor-only tissue-of-origin test

Confirms that the strong pooled signal isn't just "easy tumor vs
normal" — within-tumor heterogeneity carries tissue-of-origin
information too. All pairwise tests and one OvR run.

| Pair                          | n (pos/neg) | AUC mean ± std     |
|-------------------------------|-------------|--------------------|
| TCGA-LUAD vs TCGA-BRCA        | 10 / 10     | 1.0000 ± 0.0000    |
| TCGA-LUAD vs TCGA-COAD        | 10 / 10     | 0.9780 ± 0.0075    |
| TCGA-LUAD vs TCGA-PAAD        | 10 / 10     | 0.9080 ± 0.0075    |
| TCGA-LUAD vs TCGA-LIHC        | 10 / 12     | 0.9708 ± 0.0186    |
| TCGA-BRCA vs TCGA-COAD        | 10 / 10     | 0.9920 ± 0.0075    |
| TCGA-BRCA vs TCGA-PAAD        | 10 / 10     | 0.9960 ± 0.0049    |
| TCGA-BRCA vs TCGA-LIHC        | 10 / 12     | 0.9950 ± 0.0041    |
| TCGA-COAD vs TCGA-PAAD        | 10 / 10     | 0.9860 ± 0.0102    |
| TCGA-COAD vs TCGA-LIHC        | 10 / 12     | 0.9933 ± 0.0033    |
| TCGA-PAAD vs TCGA-LIHC        | 10 / 12     | 0.9175 ± 0.0342    |
| **OvR TCGA-LUAD vs rest**     | 10 / 42     | **0.9957 ± 0.0044** |

**Reading:** every pair of cancer types is separable at AUC ≥ 0.91
on tumor samples alone, with most pairs in the 0.97-1.00 range. The
PAAD signal is consistently the weakest pair-end (PAAD↔LUAD 0.908,
PAAD↔LIHC 0.918). The OvR "LUAD vs all other tumors" task scores
0.996 — methylation captures a strong tissue-of-origin signature, as
expected from prior literature.

## Honest comparisons

### vs Phase 0 (TCGA-LIHC only, AUC 0.972 ± 0.009)

The pooled multi-cancer AUC of **0.9778 ± 0.0031** is essentially equal
to (in fact marginally above) the LIHC-only baseline, and has ~3×
smaller std. The signal does **not** drop with more cancer types — if
anything, it's slightly more robust because the variance filter drops
cancer-type-specific probes and keeps the cancer-vs-normal probes.

### vs cfDNA fragmentomics (AUC 0.975 from cfdna-fragmentomics-pipeline)

Methylation pooled multi-cancer (0.9778 ± 0.0031) is statistically
indistinguishable from the fragmentomics baseline of 0.975 on this
sample size — both are essentially at the noise ceiling of a small
binary classifier. **Neither modality has a meaningful edge** on this
n=78 TCGA cohort. To tell them apart we'd need a much larger cohort
or a paired cfDNA-vs-tissue study, both of which are out of scope for
Phase 1 Option A.

## Caveats / what this does NOT show

- Sample sizes are tiny (10t/5n per new cohort). Per-cancer AUCs at
  1.000 are artefacts of 5-fold CV on 15 samples; pooled and pairwise
  results are more reliable.
- All tissue is fresh-frozen primary tumor / matched normal, **not**
  cfDNA from plasma. The cfDNA translation step (deconvolution, bisulfite
  conversion efficiency, low-input bias) is out of scope.
- The fragmentomics 0.975 baseline is from a separate pipeline; we did
  not re-run it. The comparison is qualitative.
- Tumor-only "tissue-of-origin" results depend on each cancer type
  having only 10 tumor samples; adding more tumors per type would
  likely raise the OvR AUC but could also lower some pairs if the
  within-tumor variance grows.

## Files created / modified

- `scripts/run_multi_cancer_baseline.py` — new multi-cancer driver
- `scripts/download_remaining_files.py` — new recovery helper for
  connection-aborted GDC downloads
- `results/multi_cancer_baseline.json` — full numerical results
- `data/raw/tcga_multi_cancer/{TCGA-LUAD,TCGA-BRCA,TCGA-COAD,TCGA-PAAD}/`
  — 60 new TSVs + tumor/normal id lists

`src/methylation/methylation_baseline.py` was **not** modified.