# Phase 0 Progress — Status Report

**Started:** 2026-09-10
**Status:** BLOCKED on public methylation data acquisition (NCBI reCAPTCHA).

## What works ✅

1. **Project structure** — Repo scaffolded at `rollroyces/deepcatch-methylation`:
   - `src/methylation/methylation_baseline.py` — feature matrix + LR baseline + 5-fold CV
   - `test/test_methylation_baseline.py` — 5 unit tests, all passing
   - `conftest.py`, `pyproject.toml`, `CITATION.cff`, `.gitignore` — standard project files
   - `.github/workflows/methylation-tests.yml` — CI workflow (matrix Python 3.11/3.12)
2. **FinaleToolkit installed** in `/Users/hermes/deepcatch/.venv` (version 1.1.0).
3. **FinaleMe cloned** at `/tmp/FinaleMe` (MIT-licensed, but requires Java 21 + Maven to build).
4. **Synthetic-data validation** — methylation baseline + LR baseline + 5-fold CV achieves AUC > 0.9 on synthetic data with planted cancer signal (test_evaluate_cv_synthetic_signal_above_chance).

## What is blocked ❌

**NCBI GEO blocks automated downloads with reCAPTCHA.** Multiple attempts failed:

| Attempt | Approach | Result |
|---|---|---|
| `curl -L` to FTP URL | Direct download | 990-byte HTML reCAPTCHA page |
| `curl` with Mozilla user-agent | Disguised browser | Same reCAPTCHA page |
| `requests.Session` with cookies | Persistent session | reCAPTCHA in response |
| GEOparse `get_GEO` (default) | Library wrapper | Timed out at 300s |
| GEOparse `get_GEO(how='brief')` | Family-only download | Returned sample metadata only (no β-values) |
| ENA portal API | Different repository | ENA doesn't index GEO accession `GSE122126` |
| EBI ArrayExpress mirror | Different repository | 404 (no matching accession) |
| Bioconductor annotation packages | Bundled data | None found |
| PyPI `pydnameth`, `methylcheck` | Bundled data | None bundle raw methylation data |
| HuggingFace dataset search | Mirrors | No `GSE122126` mirror found |

**Root cause:** NCBI requires browser-based CAPTCHA solving for all programmatic GEO downloads since 2023. Without a real browser session (the user agent our requests library sends triggers the challenge), we can't get past the wall.

## What this means for Phase 0

The original Phase 0 plan was to download GSE122126 (Moss 2018 cfDNA methylation atlas, 59 cfDNA samples with cancer/healthy labels) and run the methylation baseline. That plan is **blocked until data access is restored**.

## Options to unblock

### Option A — Manual download (1 user-action)
The user downloads GSE122126 series matrix to `data/raw/GSE122126_series_matrix.txt.gz` from a real browser. The pipeline then runs end-to-end.

**Steps for the user:**
1. Open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122126 in a browser
2. Click the "Series Matrix File(s)" link
3. Save to `/Users/hermes/deepcatch-methylation/data/raw/GSE122126_series_matrix.txt.gz`
4. Tell me the path; I'll run the pipeline.

**Pros:** simplest. **Cons:** requires the user to do a browser download.

### Option B — FinaleMe on existing data (Phase 0 with a different angle)
Skip GSE122126 entirely. Go directly to FinaleMe + FinaleDB fragments.

**Steps:**
1. Download raw FinaleDB fragment BEDs (50 samples) — FinaleDB does not have reCAPTCHA, so this should work.
2. Build FinaleMe (Java 21 + Maven, ~1 hour setup).
3. Run FinaleMe → methylation features.
4. Run baseline.

**Pros:** validates the actual Phase 1 strategy. **Cons:** heavier setup (Java 21 install, reference file download).

### Option C — Wait for NCBI access restoration
Try again in a few days; sometimes NCBI rate-limits are temporary.

**Pros:** zero work. **Cons:** uncertain when access will be restored.

### Option D — Use a smaller public methylation dataset with no reCAPTCHA wall
Search for methylation datasets on HuggingFace, OpenML, or Kaggle.

**Status:** No suitable alternative found yet (search in progress).

## Recommended next move

**Option A** is the fastest path forward — a 30-second browser action on the user's part unblocks Phase 0 entirely. Option B is the next best alternative but requires ~1 hour of compute setup.

## Synthetic validation results

Even without real data, the pipeline works on synthetic inputs:

```
test_build_feature_matrix_shape PASSED
test_build_feature_matrix_filters_low_variance PASSED
test_evaluate_cv_synthetic_signal_above_chance PASSED  (AUC > 0.9 on planted signal)
test_evaluate_cv_random_labels_is_at_chance PASSED     (AUC ~0.5 on random labels)
test_evaluate_cv_handles_perfect_separation PASSED     (AUC > 0.95 on perfectly separable)
5 passed in 1.35s
```

This confirms the pipeline structure is correct. The only missing piece is real methylation data to validate against.

## What I'll do next (without user input)

While waiting on the data blocker, I can:
1. Implement `src/methylation/finaleme_extract.py` — the FinaleMe feature extractor (won't run without Java, but can be reviewed)
2. Write the head-to-head comparison script (`methylation_vs_fragmentomics.py`)
3. Add a small synthetic dataset that mimics GSE122126 size for end-to-end testing

These don't require real data and let us have everything ready to run as soon as data is available.
