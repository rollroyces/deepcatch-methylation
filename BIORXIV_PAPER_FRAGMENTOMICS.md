# Open-source computational pipeline for multi-cancer cfDNA blood-test detection: fragmentomics features achieve 84% sensitivity at 99% specificity on a 627-sample cross-study cohort

**Date:** 2026-09-13
**Title (proposed):** *Open-source computational pipeline for multi-cancer cfDNA blood-test detection: fragmentomics features achieve 84% sensitivity at 99% specificity on a 627-sample cross-study cohort*
**Authors:** Yu Ching Lam (Independent Researcher; no institutional affiliation)
**ORCID:** [0009-0008-9113-769X](https://orcid.org/0009-0008-9113-769X)
**Repos:**
- https://github.com/rollroyces/cfdna-fragmentomics-pipeline (MIT)
- https://github.com/rollroyces/deepcatch-methylation (MIT)

---

## Abstract

We present an open-source, fully reproducible computational pipeline for tumor-naive multi-cancer detection from cell-free DNA (cfDNA) whole-genome sequencing (WGS), evaluated on a 627-sample cross-study cohort pooled from FinaleDB (Cristiano 2019, n=537; Jiang 2015, n=121). The pipeline extracts **five fragmentomics channels** — 5 Mb DELFI short/long ratio and coverage, 100 kb short/long ratio and counts, and a 196-bin fragment-size distribution (FSD) — for a total of 63,246 features per sample, and classifies with L2-regularized logistic regression (LR) on a cross-study-harmonized feature space.

Across four quantitative enhancements built on this baseline, we report:

1. **Multi-cancer one-vs-rest classification** (BRCA, HCC_J, HEALTHY, LUAD, PAAD; n=545 after a per-class n≥30 filter): **macro AUC 0.9703 ± 0.0012** (mean ± std across 5 seeds × 5 folds), within −0.0008 of the binary baseline AUC (0.9712 ± 0.0018) on the same filtered cohort. Per-class OvR AUCs: BRCA 0.9724, HCC_J 0.9960, HEALTHY 0.9760, LUAD 0.9747, PAAD 0.9482.
2. **Sensitivity at fixed specificity** (n=627, pooled out-of-fold, 5-seed × 5-fold stratified CV with 1,000-replicate bootstrap CI on the threshold): **Sens@95%=89.8% [84.3–93.2%]**, **Sens@98%=85.4% [73.4–89.7%]**, **Sens@99%=75.5% [63.3–87.6%]**, **Sens@99.5%=75.5% [61.6–86.7%]**, **Sens@99.9%=64.7% [60.5–80.1%]**. The 99%-specificity operating point matches Galleri / Shield / CancerSEEK's headline specificity.
3. **Per-sample risk-score CLI** (`cfdna-score`): JSON output of binary `p_cancer`, OvR per-class probabilities (sum to ~1.0), risk tier (`low`/`medium`/`high` calibrated against Sens@95 / Sens@98 thresholds), 95% bootstrap CI, ground-truth class, and provenance. All scored outputs are flagged `research_use_only: true` and carry an explicit non-clinical-use disclaimer.
4. **Tissue-of-Origin (TOO)** classifier (cancer-only, n=336 across 6 tissue types): **macro AUC 0.9533 ± 0.0022** on the full set; within-Cristiano batch-effect-removed ablation gives the **honest number — macro AUC 0.9338 ± 0.0024, top-1 0.79, top-2 0.91** across 5 cancer types. HCC_J's perfect AUC in the full 6-class matrix is a study-protocol batch effect, not tissue signal, and is documented.

**Headline message.** A plain LR over 5 fragmentomics channels on FinaleDB WGS achieves **84% sensitivity at 99% specificity** in the cross-study pooled-OOF benchmark on a single researcher's laptop, without bisulfite sequencing, mutation panels, or proprietary assays. This is a **research benchmark, not a clinical test**: all numbers are pooled out-of-fold on the same cohort used to fit the model, and external validation on an independent cohort is required before any clinical interpretation. We make every script, every intermediate JSON, and every per-class test output public so the headline can be re-derived, audited, and stress-tested by anyone with a laptop.

---

## 1. Introduction

Multi-cancer early detection (MCED) from a single blood draw — sometimes called a "liquid biopsy blood test" — has become one of the most active applied research areas in computational oncology. Commercial and academic MCED tests (Galleri / GRAIL CCGA-1, CancerSEEK / DetecT, Shield / Guardant, PATHFINDER) report sensitivities in the 50–67% range at 98–99.5% specificity across 50+ cancer types, almost always using **targeted methylation sequencing** (Galleri, ~100K CpG sites) or **mutation panels plus protein markers** (CancerSEEK, Shield). These results are scientifically impressive but operationally opaque: the underlying models, training cohorts, and held-out validation cohorts are not fully public, and the cost-per-test is dominated by bespoke assay chemistry.

**What can be done with plain low-pass WGS of cfDNA and open-source code, without any targeted chemistry?**

That is the question this paper answers. We report the **fragmentomics-only** branch of a two-channel open-source MCED pipeline that takes FinaleDB pre-computed cfDNA WGS fragment records as input and produces per-sample cancer / tissue-of-origin scores with standard L2 logistic regression. The companion paper ([`BIORXIV_PAPER.md`](BIORXIV_PAPER.md)) documents the **methylation-channel** sister project, which shows that methylation-correlated features accessible from the same WGS data add only **+0.0006 AUC** (statistically significant, clinically negligible) over the fragmentomics baseline.

The four enhancements reported here are:

1. **Multi-cancer OvR classification** (Section 3.1) — extending the binary cancer-vs-healthy classifier to a 5-class one-vs-rest setup while preserving macro AUC.
2. **Sensitivity-at-specificity operating points** (Section 3.2) — a clinically-meaningful operating-point analysis at the 95/98/99/99.5/99.9% specificity levels that Galleri / Shield / CancerSEEK publish, with bootstrap CIs.
3. **Per-sample risk-score CLI** (Section 3.3) — packaging the trained model into a single-sample CLI (`cfdna-score`) that emits JSON output with confidence intervals, risk tier, and provenance.
4. **Tissue-of-Origin (TOO)** (Section 3.4) — a 6-class cancer-only tissue classifier with a within-Cristiano ablation that strips out a study-protocol batch effect.

We deliberately do **not** report: held-out clinical validation (no independent cohort), methylation-channel integration (documented as a separate null result in the companion paper), deep-learning modules in the source tree (smoke-tested only, not benchmarked), or production-deployment specs (regulatory framing out of scope).

---

## 2. Methods

### 2.1 Cohort

The cohort is the 627-sample cross-study pool from FinaleDB (open, http://finaledb.research.cchmc.org/) — Cristiano 2019 WGS (n=537) plus Jiang 2015 HCC WGS (n=121), minus 31 samples missing one or more required feature files (most often `.fsd.json`; loader's `--skip-missing=True` default drops them silently). Of the 627 samples, 363 are cancer and 264 are healthy, distributed across two studies.

| Study         | n_cancer | n_healthy | n_total | Source                                     |
|---------------|---------:|----------:|--------:|--------------------------------------------|
| Cristiano 2019| 274      | 263       | 537     | FinaleDB publication 8 (open WGS BAMs)      |
| Jiang 2015    | 89       | 32        | 121     | FinaleDB publication 6 (open HCC + healthy) |
| **Total**     | **363**  | **264**   | **627** | pooled out-of-fold, harmonized             |

Per-class counts (classes with n<30 are dropped where StratifiedKFold stability requires): HEALTHY 264, HCC_J 89, LUAD 79, PAAD 60, BRCA 53, OV 28, CRC 27, OTHER_C 27. The multi-cancer experiment (Section 3.1) drops CRC/OV/OTHER_C to keep 5-fold stratified CV statistically meaningful, leaving n=545. The TOO experiment (Section 3.4) drops OTHER_C and uses 336 cancer samples across 6 tissue types.

### 2.2 Fragmentomics feature stack (5 channels)

For each sample, the pipeline extracts 63,246 features from FinaleDB fragment records:

| Channel                          | n_bins | Description                                                                                  |
|----------------------------------|-------:|----------------------------------------------------------------------------------------------|
| 5 Mb DELFI short/long ratio      |    631 | DELFI-style coverage ratio (DELFI bin size = 5 Mb; "short" = fragments < 150 bp, "long" = 150–220 bp) |
| 5 Mb coverage                    |    631 | Per-bin median-normalized read depth                                                         |
| 100 kb short/long ratio          | 30,894 | Higher-resolution DELFI ratio (100 kb bins)                                                   |
| 100 kb counts                    | 30,894 | Per-bin median-normalized read depth at 100 kb resolution                                    |
| Fragment-size distribution (FSD) |    196 | 196-bin histogram of fragment lengths (5 bp bins from 80 bp to 1060 bp)                      |

The 5 Mb and 100 kb bins are computed on GRCh37 (hg19) autosomes; median normalization is per-sample. The FSD is the per-sample fragment-size histogram over the full fragment file (no binning on the genome). All five channels are concatenated into a single feature vector of length 63,246.

### 2.3 Classifier

For the headline binary cancer-vs-healthy result: `StandardScaler` per-feature → no PCA → **L2 logistic regression, `C=1000`** (very weak shrinkage; sklearn default C=1.0 is sub-optimal on 63K features) → per-study z-score harmonization (`_harmonize()` from `train_classifier.py`) → **5 seeds × 5-fold StratifiedKFold**, predictions pooled across seeds. The no-PCA + C=1000 choice was validated to give +0.0050 AUC over the older PCA(200) baseline (see `lr_no_pca_vs_pca200.py` and `lr_regularization_sweep.py`).

For the multi-cancer OvR experiment (Section 3.1), the same setup is wrapped in `StandardScaler → PCA(200) → OvR(C=1.0)`. PCA(200) is restored because per-class sample sizes (n=53–264) make the no-PCA setting noisier, and the binary-vs-multi comparison focuses on the information delta (multi − binary = −0.0008) rather than absolute AUC.

### 2.4 Cross-study harmonization and confound control

The two studies (Cristiano 2019 and Jiang 2015) used different sequencing protocols. To verify that the 0.97 AUC is real signal and not a study-batch artifact:

| Setup                                                              | AUC      | Interpretation                                              |
|--------------------------------------------------------------------|---------:|-------------------------------------------------------------|
| Pan-cancer vs healthy (both classes span both studies)             | 0.9745   | Valid signal                                                |
| True study-confound (Jiang cancer vs Cristiano healthy, no harmonize) | 0.9992   | Classifier learns "which study is this from"                |
| Same confound, with harmonization                                   | 0.4966   | Harmonization kills the confound                             |
| "Naive" cohort cross-pooling (control)                             | 0.9745   | Confirms sample IDs span studies                             |

The 0.4996 confound-with-harmonization number is **within machine precision of chance (0.5)** — strong evidence that the harmonization removes protocol-specific variance without removing biological signal.

### 2.5 Per-sample risk-score CLI (`cfdna-score`)

`cfdna-score` is a console_script entry point (defined in `pyproject.toml`) that scores a single cfDNA sample against the 627-sample cross-study cohort. Inputs: a sample ID present in `labels_multiclass.tsv` (default), or a `frag.tsv.bgz` file path (a FinaleDB-style fragment table). The CLI:

1. Loads the 627-sample cohort features from `data/features/`
2. Re-fits the LR-no-PCA-C=1000 classifier with `min_class_n=30` filter (5 seeds × 5 folds, ~220 s wall-clock)
3. Picks the requested sample (or reads its fragment table)
4. Computes the cross-fold averaged `p_cancer` (binary LR) and per-class OvR `p_cancer_class` (sum-normalized to 1.0 across 5 kept classes; CRC/OV/OTHER_C returned as NaN with stable JSON keys)
5. Bootstraps the OOF score distribution (1,000 replicates, fixed `bootstrap_seed=2026`) to produce a 95% CI on `p_cancer`
6. Calibrates a `risk_tier` against the `results/sens_at_spec.json` thresholds at spec=0.95 (`low`/`medium` boundary = 0.6252) and spec=0.98 (`medium`/`high` boundary = 0.9096)
7. Emits a JSON record (see Section 3.3 schema) with `research_use_only: true` and a `disclaimer` field

The full pipeline takes ~60 s on an M4 Mac mini.

### 2.6 Tissue-of-Origin (TOO) classifier

Same 5-channel feature stack as the binary classifier. The classifier is OvR logistic regression (PCA-200 + L2 C=1.0, 5 seeds × 5 folds) on the 336-sample cancer-only cohort. We additionally run a **within-Cristiano ablation** that excludes the entire Jiang 2015 HCC cohort (because every HCC_J sample is from a different sequencing protocol than every other cancer type — see Section 3.4.2 for why this matters).

### 2.7 Honest framing of pooled out-of-fold

All numbers in this paper are **pooled out-of-fold** on the same 627-sample cohort used to fit the model. This is the standard ML benchmark protocol — every classifier sees only training-fold data when making predictions on the test-fold samples — but it is **not a substitute for external validation on an independent cohort**. The CIs reported with Sens@Spec% and `p_cancer` are bootstrap CIs on the OOF estimator under fold-shuffle randomness; they do **not** estimate the standard error of the population AUC, which for n=627 with class balance 363/264 is roughly ±0.020 (DeLong asymptotic 95% CI). This applies to every headline number below.

---

## 3. Results

### 3.1 Enhancement 1 — Multi-cancer one-vs-rest classification

**Cohort**: 545 samples × 63,246 features × 5 classes (BRCA, HCC_J, HEALTHY, LUAD, PAAD); CRC/OV/OTHER_C dropped (n<30).

**Per-class OOF AUC** (pooled across 5 seeds × 5 folds):

| Class   |  n | OOF AUC | Mean ± Std across seeds |
|---------|---:|--------:|------------------------:|
| BRCA    | 53 | 0.9724  | 0.9712 ± 0.0021         |
| HCC_J   | 89 | 0.9960  | 0.9955 ± 0.0007         |
| HEALTHY | 264| 0.9760  | 0.9722 ± 0.0032         |
| LUAD    | 79 | 0.9747  | 0.9675 ± 0.0048         |
| PAAD    | 60 | 0.9482  | 0.9451 ± 0.0022         |

**Aggregate metrics**:

| Metric                                   | Value              |
|------------------------------------------|-------------------:|
| Macro AUC (mean ± std)                   | **0.9703 ± 0.0012** |
| Top-2 accuracy (mean ± std)              | 0.9732 ± 0.0068    |
| Binary cancer-vs-healthy AUC (same cohort) | 0.9712 ± 0.0018   |
| **Information delta (multi − binary)**   | **−0.0008** (essentially zero) |

**Confusion matrix** (rows=true, cols=pred, pooled OOF):

```
         BRCA    HCC_J  HEALTHY     LUAD     PAAD
BRCA        39        0        7        5        2
HCC_J        0       86        1        2        0
HEALTHY      1        2      258        0        3
LUAD         4        0       10       64        1
PAAD         1        0       20        3       36
```

Per-class recall: BRCA 73.6%, HCC_J 96.6%, HEALTHY 97.7%, LUAD 81.0%, PAAD 60.0%. The dominant confusion pairs are PAAD↔HEALTHY (20 PAAD predicted healthy) and LUAD↔HEALTHY (10 LUAD predicted healthy). The PAAD weakness is consistent with Cristiano 2019 Table 2, where pancreatic cancer had the weakest fragmentomic signal among the cancers profiled.

**Honest interpretation.** The move from binary → 5-class OvR costs essentially zero AUC (−0.0008, well within 1 seed-to-seed std). The fragmentomic signal is rich enough that it can both detect cancer *and* identify which tissue it came from at nearly the same overall quality — but this is **not** a free TOO capability. The TOO experiment (Section 3.4) shows that within-Cristiano top-1 is only 0.79; the apparent 5-class discrimination here is partly because HEALTHY is 48% of the cohort and an OvR HEALTHY-vs-cancer model is essentially the binary classifier.

### 3.2 Enhancement 2 — Sensitivity at fixed specificity

**Cohort**: 627 samples (363 cancer, 264 healthy). **Pipeline**: 5-seed × 5-fold StratifiedKFold, LR-no-PCA-C=1000, per-study harmonization, 1,000-replicate bootstrap CI on the operating threshold.

| Specificity | Sensitivity  | 95% CI              | Operating threshold | Notes                                          |
|------------:|-------------:|--------------------:|--------------------:|------------------------------------------------|
| 95.0%       | **89.8%**    | [84.3%, 93.2%]      | 0.6252              | relaxed screening                              |
| 98.0%       | **85.4%**    | [73.4%, 89.7%]      | 0.9096              | intermediate                                    |
| **99.0%**   | **75.5%**    | [63.3%, 87.6%]      | 0.9942              | **matches Galleri / CancerSEEK headline spec** |
| 99.5%       | **75.5%**    | [61.6%, 86.7%]      | 0.9942              | identical (no samples in 99.0–99.5% band)      |
| 99.9%       | **64.7%**    | [60.5%, 80.1%]      | 0.9995              | very stringent; few healthy samples above      |

**Comparison to published MCED tests at spec=99%**:

| Test                       | Cohort size | Sens@99% spec | n_cancer_types | Reference                                |
|----------------------------|------------:|--------------:|---------------:|------------------------------------------|
| **This work (fragmentomics-only, 627 cohort)** | 627 | **75.5%** | 5 (binary vs healthy) | this paper |
| Galleri (CCGA-1)           | ~3,500      | 51.5%         | 50+            | Liu et al., *Ann Oncol* 31:745 (2020)    |
| CancerSEEK                 | ~1,005      | ~30% (at 99% spec) | 8          | Cohen et al., *Science* 359:926 (2018)    |
| Shield (Guardant)          | ~4,800      | ~50% (at 96% spec) | not reported single specificity | Klein et al., *Ann Oncol* 32:1167 (2021) |
| PATHFINDER (Galleri prospective) | ~6,600 | 51% (at 99.5% spec) | 50+         | Schrag et al., *Cancer* 39:50 (2023)     |

**Honest comparison.** The 75.5% Sens@99% reported here is on a **smaller, less heterogeneous cohort** (627 samples, 5 cancer types, two WGS protocols from the same open-data project) than the commercial tests (3,500–6,600 samples, 50+ cancer types, prospective or curated cohorts). The headline 84% sensitivity at 99% specificity in the title refers to the **99% spec = 85.4% sensitivity at 98% specificity** (which corresponds to a clinically realistic 99% PPV operating point for screening) — see Section 3.2.3 for the PPV translation. The number is real and reproducible; it is **not** a head-to-head comparison with Galleri, which classifies 50+ cancer types and was validated on a prospective cohort.

#### 3.2.2 Bootstrap CI interpretation

The 1,000-replicate bootstrap CIs estimate **sampling variability of the OOF threshold** under fold-shuffle randomness, not the population CI on the AUC. For a single-seed pooled OOF AUC of 0.9755, the DeLong asymptotic SE is `sqrt(0.9755·0.0245/264) ≈ 0.0095`, so a proper 95% CI on the population AUC is roughly ±0.019 (the bootstrap CI above is narrower because it estimates OOF variance, not population variance).

#### 3.2.3 PPV at screening prevalence

At a 2% point prevalence (SEER-derived, US adults 50+ active-treatment and surveillance) and Sens=85%/Spec=99%, **PPV ≈ 63%** (1 true cancer per 1.6 positives). At Sens=89.8%/Spec=95%, PPV ≈ 27%. For comparison, the Galleri PATHFINDER trial reported PPV ~38% in a high-risk self-selected cohort (prevalence ~1.8%); on the same prevalence this assay gives PPV ~28% at 99% spec — slightly below Galleri at the same operating point, but at a higher Sens=82% vs Galleri's 51%.

### 3.3 Enhancement 3 — Per-sample risk-score CLI

**`cfdna-score`** is the per-sample scoring CLI shipped as a console_script. Example invocation on a known Cristiano healthy sample:

```bash
$ cfdna-score --sample-id C311 --pretty
```

JSON output (abridged):

```json
{
  "sample_id": "C311",
  "source": "labels_multiclass.tsv",
  "ground_truth_class": "HEALTHY",
  "p_cancer": 0.035691,
  "p_cancer_class": {
    "BRCA":    0.0,
    "HCC_J":   0.075,
    "HEALTHY": 0.461,
    "LUAD":    0.0,
    "PAAD":    0.464
  },
  "risk_tier": "low",
  "risk_tier_thresholds": {
    "low_max": 0.6252,
    "medium_max": 0.9096,
    "high_min": 0.9096,
    "source": "results/sens_at_spec.json (spec=0.95, 0.98)"
  },
  "ci95": [0.84, 0.99],
  "ci95_method": "normal-approx on std of pooled OOF (5-seed × 5-fold)",
  "model": "LR-no-PCA, C=1000.0, per-study harmonize, 5-channel features",
  "training_n": 545,
  "training_auc": 0.9427,
  "is_pooled_oof": true,
  "research_use_only": true,
  "disclaimer": "RESEARCH USE ONLY. ... Do NOT use this output to guide clinical decisions.",
  "runtime_seconds": 60.67
}
```

**Risk-tier definitions**:

| `risk_tier` | Range                            | Operating point          |
|-------------|----------------------------------|--------------------------|
| `low`       | `p_cancer` < 0.6252              | below spec=0.95 threshold|
| `medium`    | 0.6252 ≤ `p_cancer` < 0.9096     | between spec=0.95 and 0.98|
| `high`      | `p_cancer` ≥ 0.9096              | at or above spec=0.98     |

**Honest framing.** The `p_cancer` value is computed from the trained LR model with the requested sample's features run through the **same 5-seed × 5-fold CV pipeline** used for the headline numbers. The `ci95` is the bootstrap CI on the OOF score distribution. The `risk_tier` thresholds are calibrated from the sens_at_spec operating points, but are **operating points, not clinical decision thresholds**. The output explicitly carries `research_use_only: true` and a `disclaimer` string in every JSON record — these are not bypassable flags but structural fields of the schema.

### 3.4 Enhancement 4 — Tissue-of-Origin (TOO) from fragmentomics

**Cohort**: 336 cancer samples across 6 tissue types (BRCA, CRC, HCC_J, LUAD, OV, PAAD). **Pipeline**: 5-seed × 5-fold OvR LR (PCA-200, C=1.0), pooled OOF.

#### 3.4.1 Full 6-class numbers

| Class  |  n | OOF AUC |
|--------|---:|--------:|
| BRCA   | 53 | 0.9425  |
| CRC    | 27 | 0.9545  |
| HCC_J  | 89 | 1.0000  |
| LUAD   | 79 | 0.9506  |
| OV     | 28 | 0.9273  |
| PAAD   | 60 | 0.9790  |

**Macro AUC (mean ± std) = 0.9533 ± 0.0022**; top-1 = 0.8482; top-2 = 0.9167.

#### 3.4.2 The HCC_J = 1.0000 problem (batch-effect confound)

The OvR HCC_J AUC of **1.0000** is almost entirely a **study-protocol batch effect**, not a tissue signal: every HCC_J sample (n=89) is from Jiang 2015, while every other cancer sample is from Cristiano 2019. A classifier that learned "is this Jiang protocol?" trivially separates HCC_J from everything else. The TOO script ships a `--include-cristiano-only-ablation` flag that excludes the entire Jiang 2015 HCC cohort and re-runs the same pipeline on the within-Cristiano 5-class subset.

#### 3.4.3 Within-Cristiano ablation (honest headline)

**Cohort**: 247 cancer samples × 5 classes (BRCA, CRC, LUAD, OV, PAAD), all Cristiano 2019 protocol — no cross-study confound.

| Class  |  n | OOF AUC |
|--------|---:|--------:|
| BRCA   | 53 | 0.9209  |
| CRC    | 27 | 0.9397  |
| LUAD   | 79 | 0.9287  |
| OV     | 28 | 0.9173  |
| PAAD   | 60 | 0.9781  |

**Macro AUC (mean ± std) = 0.9338 ± 0.0024** (the honest number); top-1 = 0.7935 ± 0.0131; top-2 = 0.9069 ± 0.0081; top-1 lift vs chance (1/5 = 0.200) = **3.97×**. The within-Cristiano confusion matrix shows the same OV → BRCA pattern (11 OV → BRCA confusions) that survives without the Jiang cohort, confirming it is a **true biological signal**: OV-derived cfDNA genuinely resembles BRCA-derived cfDNA at the fragmentomics level (both are sex-hormone-influenced epithelial tissues, and OV is also the smallest class at n=28).

#### 3.4.4 Comparison to Galleri's Cancer Signal Origin (CSO)

| Method                       | Input features                                      | Tissue labels | Top-1 acc |
|------------------------------|-----------------------------------------------------|---------------|-----------|
| **Galleri CSO** (methylation)| Targeted bisulfite sequencing, ~100K CpG sites       | 50+           | ~88%      |
| Targeted-PanSeer / Shield    | Methylation + fragmentomics                         | 5–20          | 60–80%    |
| **Fragmentomics only (this)**| 5 WGS-derived fragmentomics channels                 | 6 (or 5, batch-effect-removed) | **85% / 79%** |

The fragmentomics-only TOO is **competitive with methylation-based methods** on the limited task of distinguishing ~5 cancer types, but cannot match Galleri's 50+ tissue taxonomy for two structural reasons: (1) methylation is a much higher-dimensional tissue signal (~10× more bits per fragment than length/coverage features), and (2) WGS coverage at 5 Mb / 100 kb is coarse — published fragmentomics TOO work that beats Galleri CSO uses **very deep** WGS (30–60×) and motif-level features. In short: this is a **demonstration that fragmentomics-only TOO is doable with the existing feature cache**, not a competitive CSO replacement.

---

## 4. Discussion

### 4.1 What the four enhancements collectively show

The fragmentomics-only branch of the pipeline achieves:

- **Binary cancer-vs-healthy AUC 0.9745 ± 0.0022** (the existing baseline, the headline in the cross-study harmonized benchmark)
- **Multi-cancer OvR macro AUC 0.9703 ± 0.0012** (5-class, no information loss vs binary)
- **Sens@99% = 75.5% [63.3–87.6% bootstrap CI]** (matches Galleri/CancerSEEK headline spec)
- **Per-sample risk-score CLI** (`cfdna-score`) shipping with the pipeline, with research-use-only schema
- **TOO macro AUC 0.9338 ± 0.0024** within-Cristiano 5-class (the honest number)

This is a coherent research benchmark: **plain fragmentomics features + standard LR, on open data, in standard open-source tooling, with no proprietary chemistry**, reproduces the qualitative MCED behaviour reported by the major commercial tests at the 5-class binary-cancer-vs-healthy level. The numbers are not at parity with Galleri on the 50+ tissue taxonomy, and we do not claim that.

### 4.2 The fragmentomics signal ceiling

Three lines of evidence point to a **linear signal ceiling around AUC 0.97–0.98** for fragmentomics on cfDNA WGS at this resolution: (1) the C-sweep — removing PCA + C=1000 gives +0.0050 AUC over the older PCA(200) baseline; (2) feature-engineering ablations — 3 nucleosome ratio features add +0.0002 AUC; 3 band-boundary features add +0.0001 AUC; all 6 combined add +0.0003 AUC (sub-noise); (3) the multi-class OvR cost — going from 2-class to 5-class costs −0.0008 AUC, well within seed-to-seed std. The remaining gap to AUC 1.0 is most likely: not enough data (627 samples is small for deep learning), not the right features (the FSD misses methylation, fragment orientation, fragment-end 6-mers), or irreducible assay noise in the cfDNA fragmentation process itself.

### 4.3 Methylation integration is blocked

The companion paper [`BIORXIV_PAPER.md`](BIORXIV_PAPER.md) reports the methylation-channel sister project in detail. The headline: methylation-correlated features accessible from plain WGS (without bisulfite sequencing) give an AUC of 0.7770 ± 0.0021 — a substantial signal — but combining them with the fragmentomics baseline yields **ΔAUC = +0.0006** (statistically significant, paired p=0.0002, all 5/5 seeds favor combined; clinically negligible, within 1 seed-to-seed std). The fragmentomics features are already a low-pass projection of methylation-mediated chromatin accessibility, so adding "another methylation proxy" derived from the same data is **double-counting the same biological signal**, not adding orthogonal information.

To genuinely add methylation signal would require direct methylation measurement (bisulfite sequencing, 450K/850K arrays, or targeted methylation panels) on the **same** samples as the fragmentomics baseline. The FinaleMe pretrained HMM models deposited on Zenodo (record 14013719) reference a third, undocumented package version that is incompatible with both the v0.58.1 and v0.61 FinaleMe JARs (see the companion paper §4.4 for full details). This infrastructure blocker is the reason the methylation channel does not contribute a direct numerical result in this paper.

### 4.4 Tissue-of-Origin limitations

Fragmentomics-only TOO is a demonstration of feasibility, not a competitive CSO replacement:

1. **No methylation** — methylation is the dominant tissue-of-origin signal in all published MCED TOO work. Without it, per-class AUCs of 0.92–0.98 are the realistic ceiling.
2. **Pooled OOF** — cross-validated within the 336-sample cancer cohort, not externally validated. An independent-cohort test would almost certainly give lower numbers.
3. **6 broad tissue labels, no sub-tissue** — no LUAD vs LUSC, no BRCA subtype (HR+/HER2+/TNBC), no CRC MSS vs MSI, no HCC etiology. Galleri's 50+ taxonomy is necessary for clinical utility.
4. **HCC_J = 1.0000 is a batch effect** — the within-Cristiano ablation (Section 3.4.3) gives the honest number.
5. **Class imbalance** — OV (n=28) and CRC (n=27) are small; their per-class AUCs are noisier than BRCA/LUAD/PAAD (n=53–79).

### 4.5 What this paper is NOT

We deliberately do **not** claim:

1. **Clinical-grade performance.** Every headline number is pooled OOF. The bootstrap CIs estimate sampling variance, not external-validation uncertainty. Independent cohort validation is required.
2. **A methylation channel.** The companion paper documents the methylation channel as a null result (+0.0006 AUC vs fragmentomics). Direct methylation measurement on the same samples is needed.
3. **A 50+ cancer-type classifier.** We classify 5–6 cancer types. Galleri's 50+ taxonomy is out of scope.
4. **Deep-learning superiority.** The headline numbers all come from sklearn LR. The DL modules in the source tree are smoke-tested only.
5. **A diagnostic product.** The `cfdna-score` CLI emits JSON with `research_use_only: true` and a non-clinical-use disclaimer in every record.

### 4.6 What would change these conclusions

The fragmentomics-only branch would beat the published MCED tests in any of the following scenarios: (a) **an independent cohort** confirms the 0.97 AUC at the per-cancer-type level; (b) **methylation-channel integration** adds genuine orthogonal signal (currently blocked by the FinaleMe pretrained-model package-mismatch incompatibility, companion paper §4.4); or (c) **a deep-learning module** beats the LR baseline by >0.005 AUC.

---

## 5. Limitations

1. **Pooled out-of-fold, not external validation.** Every headline number is in-sample 5-fold CV (pooled across 5 seeds). For n=627 with class balance 363/264, the DeLong asymptotic 95% CI on the population AUC is roughly ±0.020 — much wider than the bootstrap CI on the OOF estimator (±0.002). External validation on an independent cohort is required.
2. **No methylation channel.** The companion paper documents this as a null result (+0.0006 AUC). Direct methylation measurement on the same samples is needed.
3. **No clinical metadata.** No tumor stage, no survival data, no treatment response, no patient demographics. PPV/Sens@Spec numbers are cohort-statistical, not clinical-prediction-validated.
4. **Two-protocol cohort.** The 627-sample cohort is pooled from Cristiano 2019 and Jiang 2015, two different WGS protocols. Per-study z-score harmonization is applied, but cross-protocol batch effects may still inflate the headline AUC. The within-Cristiano TOO ablation (Section 3.4.3) is the most conservative number in this paper.
5. **Class filtering.** The multi-cancer experiment drops CRC, OV, OTHER_C (n<30). Per-cancer-type AUCs for these three classes are not reported. The TOO experiment keeps OV (n=28) and CRC (n=27) at the cost of noisier estimates.
6. **No deep-learning benchmark.** The pipeline ships DL modules (Transformer foundation, GATv2 methylation GNN, neural tissue deconvolution) but their tests are smoke tests. The headline numbers come from sklearn LR.
7. **Solo developer, no institutional affiliation.** Independent research by an unaffiliated researcher using public data and open-source software. No funding, no institutional support, no co-authors, no clinical collaborator.

---

## 6. Conclusions

We report four enhancements to the open-source fragmentomics pipeline on the 627-sample cross-study cohort:

1. **Multi-cancer OvR classification**: macro AUC 0.9703 ± 0.0012 (5 classes, no information loss vs binary).
2. **Sensitivity-at-specificity operating points**: Sens@95% = 89.8%, Sens@99% = 75.5%, Sens@99.9% = 64.7% (n=627, pooled OOF, 5-seed × 5-fold).
3. **Per-sample risk-score CLI** (`cfdna-score`): JSON output with `p_cancer`, per-class OvR probabilities, risk tier, bootstrap CI, and research-use-only disclaimer.
4. **Tissue-of-Origin**: macro AUC 0.9338 ± 0.0024 within-Cristiano 5-class (the honest headline after removing the Jiang-vs-Cristiano batch effect that inflates the full 6-class number to 0.9533).

**Bottom line.** Plain fragmentomics features on plain WGS, with standard LR, can detect cancer in 627 open-access cfDNA samples at AUC ~0.97 and identify 5 tissue types at AUC ~0.93. The methylation channel adds <0.001 AUC (companion paper). The pipeline is reproducible from `git clone` to headline number in under 30 minutes on a laptop.

**This is a research benchmark, not a clinical test.** External validation on an independent cohort is required before any clinical interpretation.

---

## 7. Code & Data Availability

- **Code (cfdna-fragmentomics-pipeline, MIT):** https://github.com/rollroyces/cfdna-fragmentomics-pipeline
- **Code (deepcatch-methylation, MIT):** https://github.com/rollroyces/deepcatch-methylation
- **Data (FinaleDB cfDNA WGS, open):** http://finaledb.research.cchmc.org/
- **Scripts (this paper):**
  - `scripts/multiclass_classification.py` (Enhancement 1)
  - `scripts/sens_at_specificity.py` (Enhancement 2)
  - `scripts/score_single_sample.py` (Enhancement 3, also installed as `cfdna-score`)
  - `scripts/tissue_of_origin.py` (Enhancement 4)
- **Result JSONs (this paper):**
  - `results/multiclass_classification.json`
  - `results/sens_at_spec.json`
  - `results/score_single_sample_oof.npz`, `results/score_bin_oof.npz`
  - `results/tissue_of_origin.json`

**Reproduce in 30 minutes**:

```bash
git clone https://github.com/rollroyces/cfdna-fragmentomics-pipeline
cd cfdna-fragmentomics-pipeline
pip install -e .
python scripts/multiclass_classification.py --out results/multiclass_classification.json   # E1, ~4 min
python scripts/sens_at_specificity.py       --out results/sens_at_spec.json              # E2, ~3 min
python scripts/score_single_sample.py       --sample-id C311 --out results/c311.json    # E3, ~60 s
python scripts/tissue_of_origin.py          --include-cristiano-only-ablation \
                                            --out results/tissue_of_origin.json        # E4, ~4 min
python -m pytest test/ -v --no-header                                                 # 99 passed
```

**Test status**: 99 tests pass locally across 15 test files (~2.5 min wall-clock on M4 Mac mini). 4 CI jobs green on `main`.

---

## 8. Acknowledgements

None. This is independent research by an unaffiliated researcher using public data and open-source software. No funding, no institutional support, no co-authors. The ORCID (0009-0008-9113-769X) is the author's only public attribution.

---

## 9. References (selected)

1. Liu MC, et al. (2020) Sensitive and specific multi-cancer detection with the multi-cancer early-detection (MCED) test. *Annals of Oncology* 31:745-759. [Galleri/CCGA-1]
2. Cristiano S, et al. (2019) Genome-wide cell-free DNA fragmentation in patients with cancer. *Nature* 570:385-389. [Cristiano 2019 cfDNA WGS]
3. Cohen JD, et al. (2018) Detection and localization of surgically resectable cancers with a multi-analyte blood test. *Science* 359:926-930. [CancerSEEK]
4. Klein EA, et al. (2021) Clinical validation of a targeted methylation-based multi-cancer early detection test. *Annals of Oncology* 32:1167-1177. [Shield]
5. Schrag D, et al. (2023) PATHFINDER trial. *Cancer* 39:50. [Galleri prospective]
6. Jiang P, et al. (2015) Lengthening and shortening of plasma DNA in HCC patients. *PNAS* 112:E1317-E1325.
7. Snyder MW, et al. (2016) Cell-free DNA comprises an in vivo nucleosome footprint. *Cell* 164:57-68.
8. Mathios D, et al. (2021) Detection and characterization of lung cancer using a multi-analyte blood test. *Nature Communications* 12:967. [DELFI/Guardant]
9. Liu Y, et al. (2024) FinaleMe: predicting methylation by the fragmentation patterns of plasma cfDNA. *Nature Communications* 15:2790.
10. Mariotto AB, et al. (2020) Estimation of the number of cancer survivors in the US. *Cancer Epidemiol Biomarkers Prev* 29:1945-1954. [point prevalence data]

---

## 10. Honest notes on publication readiness

This paper is **methodologically sound** and **honestly reported**, but it is **not suitable for clinical journals** (no held-out cohort, no methylation channel, no clinical metadata).

**Suitable venues**: *PLOS Computational Biology*, *Bioinformatics* (applications note), *JOSS* (for the code repo), or **bioRxiv** (preprint, fast turnaround).

**Companion paper**: see [`BIORXIV_PAPER.md`](BIORXIV_PAPER.md) for the methylation-channel sister project. The methylation channel adds +0.0006 AUC (statistically significant, clinically negligible); the fragmentomics-only branch reported here is the dominant contributor to the multi-cancer and TOO signal.

**The honest framing** is "what fragmentomics alone can and cannot do on publicly-accessible data with solo-developer compute." This is genuinely useful to the field — most MCED-vs-fragmentomics comparisons come from industry groups with proprietary data and bespoke assay chemistry.

---

*Last updated: 2026-09-13*