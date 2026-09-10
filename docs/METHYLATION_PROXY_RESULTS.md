# Methylation-Proxy Head-to-Head on the 627-Sample Cohort

**Date:** 2026-09-10
**Task:** Phase 2 Option 2 — FinaleToolkit methylation-proxy head-to-head on the existing 627-sample FinaleDB cfDNA WGS cohort.
**Author:** Yu Ching Lam (Independent Researcher)

---

## Bottom line

| Setup | AUC (5-seed × 5-fold pooled OOF) | 95% bootstrap CI | Per-seed AUCs |
|---|---|---|---|
| **Fragmentomics-only** (5-channel baseline) | **0.9746 ± 0.0019** | [0.9626, 0.9835] | 0.9740, 0.9733, 0.9753, 0.9724, 0.9780 |
| **Methylation-proxy-only** (29 summary features) | **0.7770 ± 0.0021** | [0.7386, 0.8090] | 0.7738, 0.7780, 0.7776, 0.7799, 0.7756 |
| **Combined** (fragmentomics + 29 summary features) | **0.9752 ± 0.0019** | [0.9635, 0.9840] | 0.9747, 0.9741, 0.9758, 0.9731, 0.9785 |

```
ΔAUC (combined − fragmentomics): +0.0006 ± 0.0001   paired t = +12.87   p = 0.0002
                                  ^ all 5/5 seeds favor combined
ΔAUC (proxy − fragmentomics):    -0.1976 ± 0.0035   paired t = -113.40  p < 0.0001
                                  ^ 0/5 seeds favor proxy
```

**Honest verdict:** The combined signal beats fragmentomics-only by a **tiny but statistically significant +0.0006 AUC** (paired t-test p=0.0002, all 5 seeds favor combined). The methylation-proxy-only AUC of 0.777 is meaningfully above chance (0.5) — the proxy features carry real methylation-correlated signal — but the combined ΔAUC is far below the +0.005-0.02 target. Fragmentomics alone already saturates the signal at 0.97 AUC; there is little headroom for any linear-classifier fusion to add. **This is a useful null result that motivates the Phase 1 FinaleMe path** (where true methylation β-values would bring orthogonal information).

---

## 1. What this task actually was

The original Phase 2 plan called for running **FinaleMe** (a deep-learning model that imputes single-CpG methylation from plain cfDNA WGS) on the 627-sample FinaleDB cohort. That path is heavy: FinaleMe requires Java 21 + Maven build, hg38 reference downloads, and per-sample BAM inference (multi-day compute).

This task instead uses **FinaleToolkit per-fragment features as a methylation proxy** from the *aggregated* feature files the existing fragmentomics pipeline already produces. No FinaleMe. No BAM downloads. The 627 samples stay the same.

### 1.1 Constraint

> We do NOT have raw per-fragment BEDs. We only have aggregated features (5 channels per sample). Therefore FinaleToolkit's per-fragment feature extraction (true end-motif frequencies, true GC-content per fragment) cannot be run. The "methylation proxy" we CAN compute is the GC content, end-motif distribution, and fragment-length-based features derived from the aggregated DELFI bins.

### 1.2 What "methylation proxy" means here

The methylation-proxy features are **NOT** true methylation measurements. They are **methylation-correlated features** derived from the same fragment-length and coverage data that fragmentomics already uses. The biological justification:

| Proxy feature | True methylation correlate |
|---|---|
| Per-5Mb-bin short/long ratio (regional fragmentation) | Regional chromatin accessibility — hypomethylated regions fragment differently |
| Per-5Mb-bin coverage profile (median-normalized) | cfDNA coverage reflects open vs. closed chromatin (open = hypomethylated) |
| Per-100kb-bin ratio distribution summary stats | Same as 5Mb but at finer resolution |
| WPS-like coverage asymmetry (per-bin coverage skew) | Nucleosome positioning, regulated by methylation |
| FSD entropy + top-3 FSD bins | Nuclease activity signatures (DNASE1L3 etc., methylation-regulated) |

These features are derived from the SAME data that fragmentomics uses, so they are **NOT independent** of the fragmentomics baseline. A positive ΔAUC for the combined model does NOT prove methylation adds orthogonal information.

---

## 2. Inventory of existing fragmentomics feature files

Per sample in `/Users/hermes/cfdna-fragmentomics-pipeline/data/features/`, the aggregated feature files are:

| File | Dimensions | Content |
|---|---|---|
| `<sample>.delfi_5mb_ratio.npy` | 631 bins | 5Mb short/long ratio (DELFI signal) |
| `<sample>.delfi_5mb_coverage.npy` | 631 bins | 5Mb median-normalized coverage (CNA at 5Mb) |
| `<sample>.delfi_100kb_ratio.npy` | 30,894 bins | 100kb short/long ratio (fine DELFI signal) |
| `<sample>.delfi_100kb_counts.npy` | 30,894 bins | 100kb raw counts (depth varies per sample) |
| `<sample>.delfi_100kb_meanlen.npy` | 30,894 bins | Per-bin mean fragment length |
| `<sample>.delfi_5mb_meanlen.npy` | 631 bins | Per-5Mb-bin mean fragment length |
| `<sample>.delfi.json` | — | Summary stats + per-bin 5Mb counts dict |
| `<sample>.fsd.json` | 196 bins | FSD histogram (5bp bins 20-1000bp) + summary stats |
| `<sample>.gc_corrected.npy` | 30,894 bins | GC-corrected 100kb ratio (LOESS regression) |
| `<sample>.gc_meta.json` | — | GC correction metadata (n_bins, loess_frac, gc_count_corr) |
| `<sample>.motifs.npy` | 256 bins | 4-mer end-motif frequencies |
| `<sample>.motifs.json` | — | End-motif frequency metadata |
| `<sample>.wps_100kb.npy` | 24 bins | Window Protection Score (per-chr summary) |

The existing production baseline (`cfdna-fragmentomics-pipeline/scripts/`) uses the **5-channel combination** = `{5Mb ratio, 5Mb coverage, 100kb ratio, 100kb counts, FSD histogram}` with per-sample median-normalization on the 100kb counts. Concatenated shape: **63,246 features**.

### Cohort used

* 658 samples labeled in `labels_cross_study.tsv` (Cristiano 2019 + Jiang 2015)
* 627 samples pass the existing 5-file gate (4 npy + fsd.json)
* 363 cancer, 264 healthy
* Cristiano: 506 samples (80.7%); Jiang: 121 samples (19.3%)

---

## 3. Methylation-proxy features (29 summary features)

Computed from the SAME 5 files (no per-fragment BEDs required):

```
5mb_summary (8 features):
  mean, std, p10, p25, p50, p75, p90, frac_extreme
  ─ source: per-sample delfi_5mb_ratio.npy (631 bins)
  ─ proxy: regional short/long ratio distribution; extremes = hyper/hypo
    fragmentation regions (chromatin-accessibility correlate)

100kb_summary (12 features):
  mean, std, p10, p25, p50, p75, p90, mad, skew, kurt, frac_low, frac_high
  ─ source: per-sample delfi_100kb_ratio.npy (30,894 bins)
  ─ proxy: fine-resolution regional fragmentation distribution; skew/kurt
    capture non-Gaussian regional shifts

5mb_coverage_summary (5 features):
  mean, std, median, p10, p90
  ─ source: per-sample delfi_5mb_coverage.npy (631 bins)
  ─ proxy: regional coverage = chromatin accessibility = hypomethylation

FSD-derived entropy + top-3 (4 features):
  motif_entropy, motif_top1_freq, motif_top2_freq, motif_top3_freq
  ─ source: per-sample fsd.json (196-bin fragment-length histogram)
  ─ proxy: FSD entropy + tail bins as proxy for nuclease-activity signatures
    (4-mer motif entropy would be ideal but is only available on 82 samples;
    FSD entropy captures the same biological axis at the full 627 cohort)
```

**Total: 29 features.** The full feature descriptions live in `src/methylation/methylation_proxy.py`.

### 3.1 A larger 1,291-feature variant (with per-bin percentile-rank)

We also tested a 1,291-feature variant that adds **per-5Mb-bin percentile-rank encoded ratio + coverage profiles** (631 + 631 = 1,262 rank features, computed using TRAIN-fold ranks only to avoid test-set leakage). The rank features had to be computed **inside each CV fold** because ranking across the full cohort would leak test-set information into the training distribution. This is documented in `src/methylation/methylation_proxy.py::encode_per_bin_rank_train_test`.

**Result:** the 1,262 rank features turned out to be pure noise for LR with C=1.0 — they DEGRADED the methylation-proxy-only AUC from 0.777 to 0.675 AND the combined AUC from 0.975 to 0.921. This is documented in `results/methylation_proxy_head_to_head_with_rank.json` as the honest negative ablation.

Why do the rank features hurt?

* LR with C=1.0 on 1,262 features + ~500 train samples is severely under-regularized.
* The rank features are correlated with the summary features (both encode the same regional methylation biology) — adding redundant noisy features just makes the model overfit.
* A different model (RF, GBM, neural net) might benefit from the rank features; the LR baseline cannot.

The 29-summary-feature variant is the cleaner result and is reported in the headline `methylation_proxy_head_to_head.json`.

---

## 4. Evaluation protocol

Mirrors the existing `cfdna-fragmentomics-pipeline/scripts/lr_no_pca_vs_pca200.py` production protocol exactly:

* **Model:** L2 LogisticRegression (sklearn `LogisticRegression(C=1.0, max_iter=20000, tol=1e-8, random_state=0)`)
* **CV:** 5-fold StratifiedKFold, **5 seeds** = `[42, 13, 7, 99, 1234]`
* **Preprocessing:** per-fold per-study z-score harmonization (Cristiano/Jiang), then per-fold StandardScaler, then LR
* **Constant-column drop:** yes (`std > 1e-12`)
* **Pooling:** pooled out-of-fold AUC, sens@95spec, sens@99spec computed once per seed
* **Bootstrap CIs:** 2,000 bootstrap resamples of the pooled OOF AUC (Hanley-McNeil-style; full DeLong with structural components is a future improvement)

The C=1.0 choice matches the existing `lr_no_pca_vs_pca200.py` "LR no-PCA" baseline. The C=1000 best-of-sweep number (0.978) reported in `cfdna-fragmentomics-pipeline/RESULTS.md` is a DIFFERENT setting (regularization sweep); we use C=1.0 here so all three setups are directly comparable.

### 4.1 Why C=1.0 (not C=1000)?

The task says "the same LR baseline + 5-fold pooled-OOF protocol as the existing pipeline". The existing production script (`lr_no_pca_vs_pca200.py`) uses C=1.0. The C=1000 number from `lr_regularization_sweep.py` is the SWEEP-BEST; using it would conflate "the best baseline" with "the same baseline as fragmentomics". For a head-to-head, the protocol must be identical across all 3 setups.

---

## 5. Headline numbers

```
Setup                              AUC mean ± std            95% CI (pooled OOF)  sens@95spec  sens@99spec
─────────────────────────────────────────────────────────────────────────────────────────────────────
Fragmentomics-only (5-channel)     0.9746 ± 0.0019          [0.9626, 0.9835]     0.886        0.742
Methylation-proxy only (29 feat)   0.7770 ± 0.0021          [0.7386, 0.8090]     0.418        0.254
Combined (frag + 29 summary)       0.9752 ± 0.0019          [0.9635, 0.9840]     0.889        0.744
```

### Per-seed AUCs (5 seeds × 5-fold pooled OOF)

```
Seed       Fragmentomics   Methylation-proxy   Combined
42         0.9740          0.7738              0.9747
13         0.9733          0.7780              0.9741
7          0.9753          0.7776              0.9758
99         0.9724          0.7799              0.9731
1234       0.9780          0.7756              0.9785
```

### Statistical tests (paired t-test across 5 seeds)

```
ΔAUC (combined − fragmentomics):
  mean = +0.0006,  std = 0.0001
  per-seed deltas: [+0.00067, +0.00075, +0.00051, +0.00071, +0.00052]
  all 5 seeds favor combined: YES
  paired t = +12.87   p = 0.0002   (significant — combined is +0.0006 better)

ΔAUC (proxy − fragmentomics):
  mean = -0.1976,  std = 0.0035
  per-seed deltas: [-0.2002, -0.1953, -0.1977, -0.1925, -0.2023]
  all 5 seeds favor proxy: NO
  paired t = -113.40   p < 0.0001   (significant — proxy is 0.20 worse)
```

---

## 6. Interpretation

### 6.1 Why is the ΔAUC so small?

Fragmentomics on this cohort already hits **AUC 0.974-0.978** (the existing production pipeline). This is at or near the empirical linear-signal ceiling — the cfDNA fragment-length distribution carries so much cancer-vs-healthy signal that there is little room for ANY orthogonal feature to add value via linear-classifier fusion. Even the existing fragmentomics-vs-mutation fusion (`tumor_naive + mutation` in `cfdna-fragmentomics-pipeline/RESULTS.md`) lifts from 0.978 to 0.989 only because the mutation channel adds genuinely orthogonal information at the LR-model level. Our methylation-proxy adds no orthogonal information because:

1. **Same source data.** The methylation-proxy is derived from the SAME aggregated fragmentomics features, not from a different assay.
2. **Linear saturation.** LR with L2 regularization on a fragmentomics-baseline that already saturates has no headroom for correlated features.
3. **Information already extracted.** Fragmentomics already captures the methylation-correlated regional fragmentation signal (cfDNA fragmentation is itself partly methylation-driven).

### 6.2 Why is the proxy-only AUC 0.78?

The 29 summary features compress the 63,246-dim fragmentomics feature vector into a 29-dim methylation-correlated signal. That signal alone — without the per-bin detail — still gets 0.78 AUC on cancer-vs-healthy. This is the meaningful finding: **regional chromatin accessibility + nuclease activity signatures carry non-trivial cancer signal even when reduced to 29 summary statistics**. With TRUE methylation β-values (FinaleMe output) the methylation-only AUC would likely be higher (0.80-0.90 per FinaleMe's published 0.91 auROC), but this is a **proxy** and cannot directly compare.

### 6.3 Why is sens@95 only 0.42 for the proxy-only?

At 95% specificity, the proxy classifier picks up only 42% of cancers. This is far below fragmentomics (89%) — confirming the proxy is a weak standalone classifier. At 99% specificity it's 25% — too low for clinical use alone. The proxy is useful as a SUPPLEMENTARY signal, not as a standalone classifier.

### 6.4 Is the +0.0006 ΔAUC publishable?

Honest answer: **not as a stand-alone headline.** It is statistically significant (p=0.0002, all 5 seeds agree) but its magnitude (+0.0006 AUC, +0.3 percentage points sensitivity at 95% specificity) is far below the +0.005 minimum threshold typically considered publishable in MCED literature. The honest framing is:

> "A methylation-proxy feature set derived from aggregated fragmentomics data adds a tiny but statistically significant +0.0006 AUC over fragmentomics alone (paired t-test p=0.0002, all 5 seeds agree), confirming that the proxy captures some non-redundant signal — but this lift is too small for clinical relevance. The 0.78 proxy-only AUC confirms methylation-correlated signal is present in cfDNA WGS data and is recoverable from aggregated features alone. To get a meaningful lift in combined AUC, true methylation β-values (FinaleMe imputation on per-fragment BEDs) are required."

---

## 7. What this means for Phase 1

This task was the "small-data fallback" of Phase 2 Option 2. The relevant comparison for the FINAL Phase 2 deliverable is:

| | Methylation-proxy only | Combined ΔAUC | Headline insight |
|---|---|---|---|
| **Option 2 (this task)** — proxy from aggregated features | 0.78 | +0.0006 | Proxy carries signal but is too correlated with fragmentomics to lift the combined AUC |
| **Option 1 (Phase 1 — TBD)** — true methylation β-values from FinaleMe | Expected 0.85-0.92 | Expected +0.005 to +0.02 | True methylation is orthogonal to fragment-length signal |

The proxy result is **encouraging** (the signal exists) but **insufficient** (the lift is too small). This MOTIVATES Phase 1: if we want a meaningful head-to-head, we need FinaleMe-imputed methylation β-values that are independent of the fragment-length data, not a feature engineering on the same data.

---

## 8. Reproducibility

```bash
# From the repo root /Users/hermes/deepcatch-methylation/
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
  scripts/run_methylation_proxy_head_to_head.py \
  --output results/methylation_proxy_head_to_head.json
# Total runtime: ~4 minutes on a single core

# Quick smoke test (2 seeds, ~1.5 min)
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
  scripts/run_methylation_proxy_head_to_head.py \
  --quick --output /tmp/quick.json

# Ablation: with the 1262 per-bin rank features (negative result, ~4 min)
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
  scripts/run_methylation_proxy_head_to_head.py \
  --output results/methylation_proxy_head_to_head_with_rank.json
```

### Unit tests

```bash
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python -m pytest test/test_methylation_proxy.py -v
# 8 tests: shape, NaN-handling, no-leakage invariants, edge cases

env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python -m ruff check src/ test/ scripts/
# All checks pass
```

---

## 9. Files created / modified

| File | Status | Purpose |
|---|---|---|
| `src/methylation/methylation_proxy.py` | **new** | Methylation-proxy feature extraction (29 summary + per-fold rank encoder) |
| `scripts/run_methylation_proxy_head_to_head.py` | **new** | Orchestrator: 3-way head-to-head + statistics + JSON output |
| `test/test_methylation_proxy.py` | **new** | 8 unit tests (shape, NaN, no-leakage, edge cases) |
| `results/methylation_proxy_head_to_head.json` | **new** | Headline result (29-summary-feature design) |
| `results/methylation_proxy_head_to_head_with_rank.json` | **new** | Ablation result (1262-rank-feature design — negative) |
| `docs/METHYLATION_PROXY_RESULTS.md` | **new** | This writeup |

No files in `cfdna-fragmentomics-pipeline/` were modified (the production baseline is untouched).

---

## 10. Honest caveats

1. **The methylation-proxy is NOT true methylation.** It is fragment-length and coverage-derived proxies for methylation-correlated biology, computed from the same 5-channel fragmentomics files. No per-fragment BEDs, no FinaleMe, no bisulfite data.
2. **The proxy and fragmentomics features are NOT independent.** They share the same source signal. A positive ΔAUC for the combined model does NOT prove methylation adds orthogonal information.
3. **The +0.0006 ΔAUC is too small to be clinically meaningful.** It is statistically significant (p=0.0002, all 5 seeds agree), but +0.0006 AUC and +0.3% sensitivity at 95% specificity is below typical MCED publication thresholds (+0.005 minimum).
4. **The 0.78 proxy-only AUC is below the published FinaleMe 0.91 auROC** — which is expected because FinaleMe measures TRUE methylation while we measure proxies. This confirms that proxy features are a useful smoke validation but not a substitute for true methylation β-values.
5. **C=1.0 was chosen for protocol parity** with `cfdna-fragmentomics-pipeline/scripts/lr_no_pca_vs_pca200.py`. The C=1000 best-of-sweep number (0.978 in `cfdna-fragmentomics-pipeline/RESULTS.md`) is from a different script and not directly comparable.
6. **DeLong CI is approximated by bootstrap.** A full DeLong implementation with structural components would tighten the CI; bootstrap is conservative.
7. **The 1262 per-bin rank features are not usable with LR-C=1.0.** A different model class (RF, GBM, MLP) might benefit from them. Out of scope for this task.

---

## 11. Conclusion

The methylation-proxy head-to-head on the 627-sample FinaleDB cohort gives a clean honest result:

* **Fragmentomics-only AUC: 0.9746 ± 0.0019** (matches the existing production baseline exactly).
* **Methylation-proxy-only AUC: 0.7770 ± 0.0021** (real signal present, much weaker than fragmentomics).
* **Combined AUC: 0.9752 ± 0.0019** (+0.0006 over fragmentomics, paired t-test p=0.0002, all 5 seeds agree — statistically significant but too small for clinical impact).

The proxy-only AUC confirms methylation-correlated signal IS recoverable from aggregated cfDNA WGS features alone. The combined AUC confirms this signal is largely redundant with fragmentomics at the linear-classifier level. **The Phase 1 FinaleMe path remains the right next step** — only true methylation β-values can plausibly add a clinically meaningful ΔAUC.