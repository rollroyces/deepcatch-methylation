# Methylation Scaffolding Inventory — deepcatch & cfdna-fragmentomics-pipeline

**Date:** 2026-09-10
**Scope:** Read-only survey. No code changes made.
**Repos surveyed:**
- `/Users/hermes/deepcatch` (primary — full `src/methylation_gnn/` package)
- `/Users/hermes/cfdna-fragmentomics-pipeline` (no code; only docs reference methylation)

---

## TL;DR

The methylation capability lives entirely in **`deepcatch/src/methylation_gnn/`** (9 Python files, ~2,200 LOC plus a 127-line `CHANGES.md` and an 836-line `test_integration.py`). It is a **self-contained, documented, fully unit-tested scaffolding layer** for a GATv2-based pre-cancer "field-defect" detector.

**Key reality check:** the GNN branch is **not wired into any headline pipeline**:
- `src/fragmentomics/fusion_ablation.py` (the published fusion-ablation benchmark) does not import or reference `methylation_gnn`.
- `src/foundation/{config,data}.py` declare `"gnn": 1` as a placeholder 1-dim modality, but **no code path actually populates it** — no import of `src.methylation_gnn` anywhere outside the package itself.
- All 46 tests are **synthetic-only**; the only real methylation artefact in either repo is `results/gse185307_methylation_validation.json` (ONT cfDNA, 13 samples, LOOCV AUC = **0.357** — below chance, i.e. the existing methylation module *failed* validation).

**Production-readiness status: scaffolding-only, not production-ready.**

---

## 1. File inventory — `/Users/hermes/deepcatch/src/methylation_gnn/`

| File | Lines | Size (B) | Purpose | Tests? |
|---|---|---|---|---|
| `__init__.py` | 82 | 2,670 | Module docstring + exports of 17 public names (config, builder, model, trainer, inference, data utils, integration adapters) | Covered by `TestImports` (2 tests) |
| `config.py` | 346 | 12,403 | `GNNConfig` dataclass (validation, device auto-detect, param estimator, `to_dict`) + `NODE_FEATURE_SPEC` dict (20 feature definitions) + three presets (`DEFAULT_`, `PROTOTYPE_`, `PRODUCTION_GNN_CONFIG`) | 5 tests in `TestConfig` |
| `graph_builder.py` | 762 | 27,734 | `RegulatoryGraphBuilder`, `GenomicRegion`, `REGION_TYPES`/`EDGE_TYPES` constants; builds heterogeneous methylation graphs (5 region types × 5 edge types) with BED parsing + co-methylation/co-fragmentation edge estimation; graceful fallback when Hi-C missing | 8 tests in `TestGraphBuilder` |
| `data.py` | 619 | 22,310 | `ReferenceDataCatalog` of public dataset URLs (UCSC CpG, GENCODE, FANTOM5, 4DN Hi-C, ENCODE, TCGA, Roadmap) + `preprocess_methylation_betas`, `build_comethylation_matrix`, `parse_regions_from_bed`, `generate_cpg_regions`. **Download URLs only — no fetching implemented.** | 10 tests in `TestDataUtils` |
| `gnn_model.py` | 640 | 21,781 | `MethylationGNN` (GATv2Conv × 3 layers, 4 heads, dual head), `ReconstructionDecoder`, `AnomalyHead`, `FieldDefectLoss`, `compute_field_defect_score`, heterogeneous-edge aggregation with homogeneous fallback | 7 tests in `TestGNNModel` |
| `gnn_trainer.py` | 686 | 24,419 | `GNNTrainer`, `GNNTrainerPhase` enum, `TrainingHistory`; 2-phase training (self-supervised masked-node pretrain → joint reconstruction+anomaly finetune); MPS/CUDA/CPU, AMP, early stopping, checkpointing | 5 tests in `TestGNNTrainer` |
| `gnn_inference.py` | 366 | 11,071 | `GNNInference` (low-level: load checkpoint, predict, predict_details with node-level anomaly maps) + `MethylationGNNPredictor` (high-level: builder + checkpoint) | 4 tests in `TestGNNInference` |
| `integration.py` | 408 | 13,440 | `MethylationBranchAdapter` + `ModularArmsBuilder` + `extend_fusion_with_gnn`; bridges GNN output → `CrossAttentionFusion` as a 5th modality | 5 tests in `TestIntegration` |
| `test_integration.py` | 836 | 33,778 | Self-contained unittest suite (no pytest fixtures), with graceful skip when torch/PyG missing; **46 tests** across 8 test classes | — (the suite itself) |
| `CHANGES.md` | 127 | 5,900 | Detailed v2.1 changelog: rationale, architecture diagram, integration pseudocode, test coverage claims, **explicit list of remaining data-acquisition work** (UCSC, GENCODE, FANTOM5, GM12878 Hi-C, ENCODE, TCGA) | — |

**Total:** ~4,872 LOC (code) + 127 LOC (changelog) + 836 LOC (tests). Compiled `.pyc` files all from Python 3.14.

---

## 2. Per-file summaries

### `__init__.py` (82 lines)
Package docstring describes the GNN as Stage 1 (Capture) of the DeepCatch CET pipeline, citing 5 references (Rao 2014 Hi-C, Brody 2022 GATv2, Velickovic 2018 GAT, Schlichtkrull 2018 R-GCN, Hu 2020 self-supervised pretraining). Re-exports 17 public names grouped into config, graph builder, model, trainer, inference, reference data, and fusion adapters.

### `graph_builder.py` (762 lines)
Heterogeneous methylation-graph construction. `GenomicRegion` dataclass holds chrom/start/end/region_type/metadata. `RegulatoryGraphBuilder` consumes a methylation dict (β-values, CpG density, GC, coverage) and an optional cfDNA coverage array and emits a PyG `Data` object with `x` (20-dim node features) and `edge_index` (5 edge-type aware). Falls back to genomic-proximity edges when reference Hi-C is absent; falls back to fully synthetic regions when none are provided. Constants: `REGION_TYPES` (cpg_island/enhancer/promoter/ctcf/dhs) and `EDGE_TYPES` (physical_interaction/co_methylation/co_fragmentation/genomic_proximity/regulatory_domain).

### `data.py` (619 lines)
Reference-data metadata layer. `ReferenceDataCatalog` is a dataclass-of-records with URLs and sizes for UCSC CpG islands, GENCODE v44 promoters, FANTOM5 enhancers, 4DN Hi-C, ENCODE chromatin, TCGA methylation, Roadmap Epigenomics, GeneHancer. **All download functions return `(url, suggested_path)` tuples — the module deliberately does not perform network I/O.** Includes `preprocess_methylation_betas`, `build_comethylation_matrix`, `parse_regions_from_bed`, `generate_cpg_regions` for offline preprocessing.

### `gnn_model.py` (640 lines)
Core architecture: `MethylationGNN` = feature projection (20→64) → 3× GATv2Conv (64→128→256, 4 heads, BatchNorm + ReLU + 0.2 dropout) with heterogeneous-edge aggregation (one GATv2Conv per edge type, mean-pooled; falls back to one homogeneous GATv2Conv when `n_edge_types=1`) → dual head: `ReconstructionDecoder` (MLP 256→128→64→20 for masked-node prediction) + `AnomalyHead` (MLP 256→128→64→1 → sigmoid per-node anomaly). `FieldDefectLoss` combines reconstruction MSE + BCE anomaly + optional temporal consistency; `compute_field_defect_score` reduces per-node anomalies to a single graph-level scalar via weighted reconstruction error.

### `gnn_trainer.py` (686 lines)
Training pipeline. `GNNTrainer` handles the two phases advertised in CHANGES.md (Phase 1 pretrain, Phase 2 finetune; the third "head-only fine-tuning" phase is in the enum but not in `pretrain()`/`finetune()`). Uses AdamW + CosineAnnealingLR or ReduceLROnPlateau. Mixed-precision AMP on CUDA. Checkpoints go to `checkpoints/methylation_gnn/` (configurable). Returns `TrainingHistory` dataclass with per-epoch loss/AUC/AUPRC. `evaluate()` returns `{"auc", "loss_recon", "loss_anomaly"}` and **only checks key presence** — see §3.

### `gnn_inference.py` (366 lines)
Lightweight inference wrapper. `GNNInference.load(path)` reconstructs model + config from a checkpoint dict; `predict(graph)` returns scalar score; `predict_details(...)` returns score + node-anomaly vector + reconstruction-error heatmap. `MethylationGNNPredictor` bundles `RegulatoryGraphBuilder` + loaded model so callers can do `predictor.predict_sample(name, methylation_data, cfDNA_coverage)`.

### `integration.py` (408 lines)
Bridges the GNN branch to `src/multimodal_fusion/advanced_fusion.py`. `MethylationBranchAdapter(checkpoint, n_regions, device)` wraps builder + predictor; `predict_batch(samples)` returns `(n_samples,)` array of field-defect scores matching the interface expected by `CrossAttentionFusion.fit(...)` and `EarlyLateFusion.fit(...)`. `ModularArmsBuilder` is a composite extractor for all Stage 1 modalities (frag, cnv, sero, mfr, **gnn**) producing standardised scores. `extend_fusion_with_gnn(frag_scores, cnv_scores, sero_scores, mfr_scores, gnn_scores)` returns the 5-modality list. **All wiring is local — no call site uses these adapters in production code.**

### `config.py` (346 lines)
20-entry `NODE_FEATURE_SPEC` documenting every input feature (mean_methylation, methylation_entropy, methylation_variance, cpg_density, cpg_obs_exp, gc_content, coverage_depth, fragment_size_mean, fragment_short_frac, end_motif_diversity, then 5 chromatin marks + 5 region-type one-hot). `GNNConfig` dataclass with validation (`__post_init__`), auto device detection (CUDA → MPS → CPU), `estimated_params` property (calculator targeting ~2M params), and `to_dict` for serialization. Three preset configs: `DEFAULT_GNN_CONFIG` (50K nodes, full-size), `PROTOTYPE_GNN_CONFIG` (5K nodes, for tests), `PRODUCTION_GNN_CONFIG` (100K nodes, larger hidden dims).

### `CHANGES.md` (127 lines)
Polished v2.1 changelog with: rationale (fragmentomics has a 0.01% detection floor, methylation field defects precede any cancer cell), ASCII architecture diagram, key design decisions table, integration pseudocode showing `n_modalities=4 → 5` migration, dependencies added (`torch >= 2.0`, `torch_geometric >= 2.5`), test-coverage checklist, and an **explicit "Remaining Work (Data Acquisition)" section** listing the 6 dataset downloads still needed (UCSC CpG, GENCODE v44, FANTOM5, GM12878 Hi-C, ENCODE, TCGA). Rollback note: "simply don't import from `src.methylation_gnn`" — confirming the module has no side effects on existing code.

---

## 3. Test inventory — `src/methylation_gnn/test_integration.py`

**Total: 46 tests across 8 classes.** All run via stdlib `unittest` (no pytest fixtures); module is importable without torch/PyG (graceful skip).

| Class | # tests | PyG required? | What it covers |
|---|---|---|---|
| `TestConfig` | 5 | No | default config fields, prototype config, `to_dict`, invalid-input `ValueError`, `NODE_FEATURE_SPEC` length=20 |
| `TestGraphBuilder` | 8 | Some | synthetic-region fallback, `load_regions_from_list`, node features (basic / with chromatin / with no data), edge fallback when no Hi-C, `GenomicRegion.distance_to`, full PyG graph construction |
| `TestGNNModel` | 7 | Yes | model creation, forward pass, pretrain forward (masked), `predict()`, homogeneous fallback (`n_edge_types=1`), `FieldDefectLoss` total > 0, smallest-model param count |
| `TestGNNTrainer` | 5 | Yes | trainer creation, single pretrain epoch, full 2-phase cycle (2 epochs each), checkpoint save/load, predict with details |
| `TestGNNInference` | 4 | Yes | `load`, predict returns float in [0,1], `predict_details` returns node-anomaly scores, high-level predictor round-trip |
| `TestIntegration` | 5 | Some | `extend_fusion_with_gnn` shape check, `ModularArmsBuilder` without GNN, `extract_all`, `extend_then_fit_fusion`, `MethylationBranchAdapter` predict |
| `TestDataUtils` | 10 | No | reference catalog length, dataset metadata, download URLs sanity, synthetic data shape, `preprocess_methylation_betas` (with/without regions), `build_comethylation_matrix`, BED parsing, `generate_cpg_regions`, `print_catalog` |
| `TestImports` | 2 | No | module importability, all 17 `__all__` names resolvable |

### AUC / accuracy thresholds asserted

**There are no numeric AUC or accuracy thresholds asserted anywhere.** All AUC-related assertions are *key-presence* only:

- `test_full_training_cycle` (line 457): `self.assertIn("auc", metrics)` + `self.assertIn("loss_recon", metrics)`.
- `test_predict` (line 487): asserts `field_defect_score` key + `node_anomaly_scores` key exist in details dict.
- No test asserts `metrics["auc"] > 0.5` or any analogous threshold.

The pretraining / finetuning cycles use only **2 epochs each** with synthetic 20-sample datasets (`make_synthetic_methylation_data` at line 64). The "data" is `np.random.RandomState(42)` Beta-distributed with hypermethylation noise — *not* real cfDNA methylation.

### Smoke vs real-data

- **100% synthetic smoke tests.** `make_synthetic_methylation_data` is the only data source for the training/integration tests.
- **No fixtures use real .bedmethyl, .bismark, .cov, .bed, or TCGA files.** All references to external datasets are URL strings in `data.py:ReferenceDataCatalog` with no download-or-cache logic.
- Tests can be run as `python src/methylation_gnn/test_integration.py` (per CHANGES.md).

### CI integration status

`RESULTS.md` line 516 marks the directory as **"✅ in CI"** (46 tests). However, attempting `python -m pytest src/methylation_gnn/test_integration.py` against the current `.venv` fails at collection time because the venv does not have PyTorch installed — the file is *intended* to run, not actively running on this machine.

---

## 4. Other methylation references across both repos

### `deepcatch` — non-`methylation_gnn` methylation files

| Path | Type | Real or synthetic? | Wired into pipeline? | Notes |
|---|---|---|---|---|
| `results/gse185307_methylation_validation.json` | Real-data output | **Real** (GSE185307 ONT cfDNA, 13 samples) | **No** | LOOCV AUC = **0.357** (below 0.5 chance) — produced by `validation/py/two_stage_cet.py`-style tooling, **not** by the GNN branch. Sens@95% spec = 0.83 but overall AUC suggests failure. This is the only real methylation artefact on disk. |
| `src/foundation/config.py` & `data.py` | Modality dictionary entries | n/a (config only) | **No — placeholder only** | Declares `"gnn": 1` (1-dim modality for field-defect score) and a healthy/cancer prior range, but **nothing populates it**. No import of `src.methylation_gnn` anywhere under `src/foundation/`. |
| `src/multimodal_fusion/advanced_fusion.py:71` | Docstring comment | n/a | No | Mentions "methylation bins" in a heterogeneous-graph description; **no actual methylation code** in the fusion layer. |
| `src/tissue_deconv/*` | Separate module | Synthetic (with optional cfSort-style atlas hooks) | Sort of — independent methylation pipeline | `tissue_atlas.py` has synthetic tissue-specific methylation profiles; `tissue_features.py:273` and `integration.py` extract tissue fractions from methylation. **This is the only other methylation consumer, and it does not share code with `methylation_gnn`.** |
| `docs/USER_GUIDE.md`, `README.md`, `RESULTS.md`, `MODEL.md`, `USAGE.md`, `TEAM.md`, `AUDIT_REPORT*.md`, `PIPELINE_AUDIT.md`, `wiki/*`, `paper/PAPER.md`, `review/agent_review_2026-08-10.md` | Documentation | n/a | n/a | All mention methylation in prose; `README.md` and `USER_GUIDE.md` describe the GNN API at a high level, but no doc is dedicated to methylation (no `METHYLATION.md`). |

### `cfdna-fragmentomics-pipeline`

**No code files reference methylation.** Only docs:
- `USAGE.md` — Galleri PATHFINDER comparator note.
- `RESULTS.md` — mentions "GATv2 methylation GNN scaffolding" pointing back to `../../deepcatch/src/methylation_gnn/`.
- `AUDIT_REPORT_2.md` — flags head-to-head vs Galleri/Liu 2020 as missing work.
- `BENCHMARK.md` — mentions methylation as a "different feature class" worth testing.
- `TEAM.md` §2.2 — "Methylation / Galleri-style expert — **HIGH PRIORITY**"; explicitly says *"Audit the existing methylation GNN scaffolding"* is the next step; marks "Methylation channel" as **"❌ scaffolding only"** in the production-readiness summary table.

### No other artefacts found

- No `*.bedmethyl` files.
- No `*.bismark` / `*.bismark_cov` files.
- No `*methylation*.py` outside `src/methylation_gnn/` and `src/tissue_deconv/`.
- No `scripts/*methylation*` or `scripts/*gnn*`.
- No `test/test_*methylation*` or `test/test_*gnn*` (tests live exclusively inside the package).

---

## 5. Per-finding readiness assessment

| Finding | Documented? | Real-data tests? | Wired into headline pipeline? | Production-ready gap |
|---|---|---|---|---|
| **Methylation GNN architecture & training** (`gnn_model.py`, `gnn_trainer.py`) | ✅ `CHANGES.md` is thorough | ❌ All 7 model + 5 trainer tests use synthetic `np.random.RandomState(42)` data | ❌ Never invoked by `run_full_validation.py`, `fusion_ablation.py`, or any headline script | Needs: real β-value matrix → node-feature path; training on TCGA-LUAD/LIHC + matched healthy controls; AUC threshold tests |
| **Graph builder** (`graph_builder.py`) | ✅ Class docstrings; integration demo in module header | ⚠️ Tests use synthetic regions + random β-values; no real BED parsing test against UCSC CpG track | ❌ Not called by any pipeline | Needs: end-to-end test on a real *.bed.gz CpG island file; co-methylation edges from a real TCGA cohort |
| **Reference-data catalog** (`data.py`) | ✅ Inline dataset table | ⚠️ Tested for URL presence and shape of preprocessing, never for "download succeeds" | ❌ No code path downloads anything; download functions return `(url, path)` tuples only | Needs: actual download/caching layer (e.g. via pooch or explicit script); checksums; smoke-test on at least one downloaded file |
| **Inference** (`gnn_inference.py`) | ✅ Class docstrings + module example | ⚠️ Predictor round-trip is tested but on synthetic 20-region graphs | ❌ No checkpoint exists; inference cannot run end-to-end without a trained model | Needs: at least one trained checkpoint on a real cohort; integration test with `predictor.predict_sample` hitting a real methylation dict |
| **Fusion integration** (`integration.py`, `extend_fusion_with_gnn`) | ✅ CHANGES.md has integration pseudocode; `integration.py` has worked example | ⚠️ `TestIntegration` calls `extend_fusion_with_gnn` with synthetic arrays | ❌ **`src/multimodal_fusion/advanced_fusion.py` does not call it. `src/fragmentomics/fusion_ablation.py` does not import it. `src/foundation/` has a `"gnn"` placeholder key but no producer.** | Needs: a call site that loads `methylation_gnn` outputs into the actual fusion layer of `run_full_validation.py`; an ablation that includes/excludes the gnn branch and reports AUC delta |
| **GCLOUD/HEADLINE fusion_ablation result** | n/a | n/a | Uses synthetic mutation score; no GNN | Headline benchmark currently uses `(frag, cnv, sero, mfr)` + a synthetic mutation channel (marginal AUC ≈ 0.92). The GNN branch would be a 5th or 6th modality |
| **Real-data validation artefact** (`results/gse185307_methylation_validation.json`) | ✅ File exists on disk | ✅ Real ONT cfDNA methylation, 13 samples, 6 cancer / 7 healthy | ⚠️ Generated by `validation/py/` tooling, not by `methylation_gnn`; LOOCV AUC = **0.357** (failed) | A baseline methylation feature module already exists and **fails**. The GNN branch has never been benchmarked on this same data — this is the obvious next experiment |

---

## 6. What's needed to make it production-ready

Based on the explicit gap list in `CHANGES.md` plus the integration findings above:

### Critical path (no ordering assumed)

1. **Download + cache the 6 reference datasets** (UCSC CpG, GENCODE v44, FANTOM5, GM12878 Hi-C, ENCODE chromatin, TCGA 450K/850K β-values). `data.py` has URLs; there is no downloader script.
2. **Build a `methylation_features → graph` pipeline** that consumes real β-value matrices (e.g. TCGA Illumina 450K/EPIC) and emits a 50K-node graph with measured CpG density / GC / coverage / fragment-length features — currently `graph_builder.build_node_features` is exercised only on `np.random.rand()`.
3. **Train a checkpoint on real data** with held-out validation; save to `checkpoints/methylation_gnn/finetune_best.pt`. Without a checkpoint, `gnn_inference.py` is dead code.
4. **Add numeric-threshold tests**: assert `test_full_training_cycle`'s `metrics["auc"] > 0.7` (or whatever the real-data baseline supports) instead of only `assertIn("auc", metrics)`.
5. **Wire `MethylationBranchAdapter` (or `ModularArmsBuilder`) into `src/foundation/data.py`'s `"gnn"` modality slot**, so the foundation module actually populates the placeholder that already exists.
6. **Add a `gnn` ablation row to `src/fragmentomics/fusion_ablation.py`** (or a new `methylation_ablation.py`) that runs the 627-sample FinaleDB benchmark with vs without the GNN branch and reports ΔAUC / ΔSens@95.
7. **Re-run the failed `gse185307_methylation_validation.json` experiment using the GNN branch** as a follow-up to `validation/py/two_stage_cet.py` (which produced the 0.357 AUC).

### Defensive cleanup

8. **Fix `data.py` ChromHandlerMetadata gap** (no explicit error if a download fails — currently silent).
9. **Schema-validate `GNNConfig` against `NODE_FEATURE_SPEC`** — `n_node_features=20` is hard-coded but the spec has 20 entries; a regression that drops an entry would silently corrupt training.
10. **Resolve the test-suite collection error in this venv** (no `torch` installed → `python -m pytest` fails before any test runs). Either add `torch` to a CI requirements file or guard `__init__.py` imports so `import src.methylation_gnn` succeeds without PyTorch.

---

## 7. Bottom line

The methylation GNN scaffolding in `deepcatch/src/methylation_gnn/` is **well-engineered, well-documented, and well-tested at the smoke level** — but it is a **fully self-contained scaffolding layer with no production consumer**. The headline fusion pipeline (`fusion_ablation.py`, `multimodal_fusion/advanced_fusion.py`, `foundation/`) treats methylation as either an absent modality or a placeholder key. The one real-data methylation artefact in either repo (`gse185307_methylation_validation.json`, AUC = 0.36) was produced by a *different* code path that did not leverage this GNN.

To make this a usable capability rather than scaffolding, the next move is **not more code in `methylation_gnn/`** — it is **data acquisition + a fusion-ablation experiment that wires the branch into `run_full_validation.py`** (per `TEAM.md §2.2 HIGH PRIORITY` and `AUDIT_REPORT_2.md` open work item).
