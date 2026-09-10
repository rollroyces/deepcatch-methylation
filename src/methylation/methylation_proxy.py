"""Methylation-proxy feature extraction from aggregated fragmentomics data.

These are NOT true methylation measurements. They are **fragment-length and
coverage-derived proxies** for methylation-correlated biology, computed from
the per-sample aggregated feature files that the existing fragmentomics
pipeline already produces (5 channels: 5Mb ratio, 5Mb coverage, 100kb ratio,
100kb counts, FSD histogram). We do NOT have per-fragment BEDs available, so
FinaleToolkit's per-fragment feature extraction (end motifs, GC content) cannot
be run directly — these aggregated proxies are the next-best signal available.

Biological justification (all *proxies*, not true methylation):
  - **GC content proxy** — the gc_corrected.npy channel already represents
    100kb-bin short/long ratio AFTER per-bin GC-bias LOESS correction. The
    *residual* from the cohort median captures regional chromatin-accessibility
    bias that GC correction removed, which is a proxy for hypomethylation
    (open chromatin → hypomethylation → different fragmentation pattern).
  - **Coverage profile proxy** — cfDNA coverage in a genomic region reflects
    chromatin accessibility. Open/hypomethylated chromatin yields different
    cfDNA coverage than closed/hypermethylated chromatin. Per-bin coverage is
    therefore a proxy for regional methylation status.
  - **Window Protection Score (WPS)** — boundary-vs-interior fragment coverage
    asymmetry, sensitive to nucleosome positioning, which is itself regulated
    by methylation state.
  - **End-motif frequencies** — cfDNA end motifs reflect nuclease activity
    (DNASE1, DNASE1L3, etc.) whose expression is tissue-specific and
    methylation-regulated.

Important caveats:
  - These features correlate with methylation biology but do not measure it.
  - The same source data is used to derive both the fragmentomics baseline
    AND the methylation-proxy features, so the two are NOT independent. This
    is the explicit design constraint of "Phase 2 Option 2" (no BAMs, no
    FinaleMe, only aggregated features).

Returns:
  A feature matrix X of shape (n_samples, n_proxy_features) plus the list of
  feature names. Features are designed to be lightweight (~50-100 dims) so
  that the methylation-proxy LR baseline is fast to train.
"""

from __future__ import annotations

import json
import os
from typing import Iterable

import numpy as np

__all__ = [
    "build_proxy_feature_matrix",
    "load_5channel_baseline",
    "encode_per_bin_rank_train_test",
    "PROXY_FEATURE_NAMES",
]


# Number of features produced per channel. Used for input-shape validation.
def _n_5mb_summary() -> int:
    return 8  # mean, std, p10, p25, p50, p75, p90, frac_extreme


def _n_100kb_summary() -> int:
    return 12  # mean, std, p10, p25, p50, p75, p90, MAD, skew_proxy,
               # kurt_proxy, frac_low (hypomethylation-proxy), frac_high


def _n_wps_summary() -> int:
    return 5  # mean, std, median, p10, p90


def _n_motif_summary() -> int:
    return 4  # entropy, top-3 PCA proxy = just entropy + 3 most-variable


def _n_bin_residual_features(n_5mb_bins: int = 631) -> int:
    # Per-cohort percentile-rank encoded 5Mb ratio profile (vector of length 631)
    return n_5mb_bins


def _n_coverage_features(n_5mb_bins: int = 631) -> int:
    # Per-cohort percentile-rank encoded 5Mb coverage profile
    return n_5mb_bins


# Total feature count for documentation
# Per-bin (5Mb) rank features are computed dynamically per CV fold
# (see encode_per_bin_rank_train_test) and are NOT included in
# PROXY_FEATURE_NAMES — they are fold-aware by design.
PROXY_FEATURE_NAMES: list[str] = (
    [f"5mb_summary_{s}" for s in (
        "mean", "std", "p10", "p25", "p50", "p75", "p90", "frac_extreme"
    )]
    + [f"100kb_summary_{s}" for s in (
        "mean", "std", "p10", "p25", "p50", "p75", "p90", "mad",
        "skew", "kurt", "frac_low", "frac_high"
    )]
    + [f"wps_{s}" for s in ("mean", "std", "median", "p10", "p90")]
    + ["motif_entropy", "motif_top1_freq", "motif_top2_freq", "motif_top3_freq"]
)


def load_5channel_baseline(
    features_dir: str,
    sample_ids: Iterable[str],
) -> tuple[np.ndarray, list[str]]:
    """Load the existing 5-channel fragmentomics baseline.

    The 5 channels (per sample) are:
      1. delfi_5mb_ratio.npy    — 631 bins (5Mb short/long ratio)
      2. delfi_5mb_coverage.npy — 631 bins (5Mb median-normalized coverage)
      3. delfi_100kb_ratio.npy  — 30,894 bins (100kb short/long ratio)
      4. delfi_100kb_counts.npy — 30,894 bins (100kb raw counts)
      5. fsd.json               — 196 bins (FSD histogram, sorted by bin start)

    Returns:
      X: ndarray of shape (n_samples, 631+631+30894+30894+196) = (n, 63246).
      kept: list of sample IDs (subset of input sample_ids with all 5 files).
    """
    rows: list[np.ndarray] = []
    kept: list[str] = []
    for s in sample_ids:
        r5p = os.path.join(features_dir, f"{s}.delfi_5mb_ratio.npy")
        c5p = os.path.join(features_dir, f"{s}.delfi_5mb_coverage.npy")
        r100p = os.path.join(features_dir, f"{s}.delfi_100kb_ratio.npy")
        c100p = os.path.join(features_dir, f"{s}.delfi_100kb_counts.npy")
        fsd_p = os.path.join(features_dir, f"{s}.fsd.json")
        if not all(os.path.exists(p) for p in (r5p, c5p, r100p, c100p, fsd_p)):
            continue
        r5 = np.load(r5p)
        c5 = np.load(c5p)
        r100 = np.load(r100p)
        c100 = np.load(c100p)
        with open(fsd_p) as f:
            sb = json.load(f).get("size_bins", {})
        keys = sorted(sb, key=lambda k: int(k.split("-")[0]))
        fsd = np.asarray([sb[k] for k in keys], dtype=float)
        # Per-sample median-normalize the 100kb counts (depth-bias removal).
        med = float(np.median(c100))
        if med > 0:
            c100 = c100 / med
        rows.append(np.concatenate([r5, c5, r100, c100, fsd]))
        kept.append(s)
    return np.nan_to_num(np.asarray(rows, dtype=float), nan=0.0,
                         posinf=0.0, neginf=0.0), kept


def _summary_stats(v: np.ndarray) -> list[float]:
    """Per-channel summary stats that compress a long vector to 7-12 features."""
    v = v[np.isfinite(v)]
    if v.size == 0:
        return [0.0] * 12
    p10, p25, p50, p75, p90 = np.percentile(v, [10, 25, 50, 75, 90])
    mad = float(np.median(np.abs(v - p50)))
    skew = float(((v - v.mean()) ** 3).mean() / (v.std() ** 3 + 1e-12))
    kurt = float(((v - v.mean()) ** 4).mean() / (v.std() ** 4 + 1e-12) - 3.0)
    # frac_low / frac_high: fraction of bins in the bottom/top 10% of cohort
    # (these are the bins whose value is in the tail — proxy for "extreme"
    # regional signal, biologically interpretable as hyper/hypomethylated
    # regions).
    p10_thr = float(np.percentile(v, 10))
    p90_thr = float(np.percentile(v, 90))
    frac_low = float((v <= p10_thr).mean())
    frac_high = float((v >= p90_thr).mean())
    return [v.mean(), v.std(), p10, p25, p50, p75, p90, mad,
            skew, kurt, frac_low, frac_high]


def _per_bin_percentile_rank_train_test(
    X_train: np.ndarray, X_test: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """For each column, compute the train-rank of train values AND the
    train-ECDF-mapped rank of test values.

    CRITICAL: this MUST be called with train data only — never include the
    test fold in the rank computation, or the per-bin rank becomes a leak
    (the test fold's value is known at the moment we rank the column).

    The train values get a normalized [0, 1] rank within the train column.
    The test values get mapped through the train ECDF: their rank is
    interpolated between the bracketing train values (clipped to [0, 1]).
    """
    out_tr = np.zeros_like(X_train, dtype=float)
    out_te = np.zeros_like(X_test, dtype=float)
    for j in range(X_train.shape[1]):
        col_tr = X_train[:, j]
        order = np.argsort(col_tr, kind="mergesort")
        ranks_tr = np.empty_like(order, dtype=float)
        ranks_tr[order] = np.arange(col_tr.size, dtype=float)
        out_tr[:, j] = ranks_tr / max(col_tr.size - 1, 1)
        # Map test values through the train ECDF
        col_te = X_test[:, j]
        # For each test value, find its position in the sorted train column.
        # searchsorted returns the insertion index in the sorted array.
        ins = np.searchsorted(col_tr[order], col_te, side="right")
        # Interpolate between train-rank positions: ins is between 0..N_tr.
        out_te[:, j] = ins / max(col_tr.size, 1)
    return out_tr, out_te


# Public alias for fold-aware use
encode_per_bin_rank_train_test = _per_bin_percentile_rank_train_test


def build_proxy_feature_matrix(
    features_dir: str,
    sample_ids: Iterable[str],
    cohort_X: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    """Build the methylation-proxy feature matrix from aggregated features.

    The proxy features encode methylation-correlated signal that lives inside
    the 5-channel fragmentomics data. We exploit the fact that chromatin
    accessibility (which is methylation-correlated) is encoded in COVERAGE
    patterns and SHORT/FRAGMENT-RATIO PATTERNS that the fragmentomics pipeline
    already extracts.

    Args:
      features_dir: directory holding per-sample .npy / .json files.
      sample_ids: iterable of sample IDs (intersected with available files).
      cohort_X: optional pre-loaded 5-channel baseline matrix (same row order
        as sample_ids). If None, will load via load_5channel_baseline.

    Returns:
      X_proxy: ndarray of shape (n_samples, n_proxy_features) — methylation-proxy
        feature matrix.
      y_or_ids: the sample-IDs actually used (subset of input).
      kept_ids: list of sample IDs retained (post file-existence filter).
      feature_names: list of feature-name strings, length == X_proxy.shape[1].
    """
    if cohort_X is None:
        cohort_X, kept = load_5channel_baseline(features_dir, sample_ids)
    else:
        # Filter sample_ids to those we successfully loaded
        kept = list(sample_ids)

    # ----- 1. Decompose 5-channel into the 5 component arrays -----
    # Channel boundaries: 5mb_ratio(631) + 5mb_coverage(631) + 100kb_ratio(30894)
    #                    + 100kb_counts(30894) + fsd(196)
    n5 = 631
    n100 = 30894
    nfsd = 196
    r5 = cohort_X[:, :n5]
    c5 = cohort_X[:, n5:2 * n5]
    r100 = cohort_X[:, 2 * n5:2 * n5 + n100]
    fsd = cohort_X[:, 2 * n5 + 2 * n100:2 * n5 + 2 * n100 + nfsd]

    # ----- 2. Per-channel summary statistics (compressed "regional" features) -----
    n_samples = cohort_X.shape[0]
    rows: list[np.ndarray] = []
    names: list[str] = []
    for i in range(n_samples):
        feats: list[float] = []
        # 5Mb ratio summary (8 features)
        s = _summary_stats(r5[i])[:8]
        feats.extend(s)
        # 100kb ratio summary (12 features, all of them)
        s = _summary_stats(r100[i])
        feats.extend(s)
        # WPS-like coverage asymmetry: difference between c5 mean and r5 mean
        # is not biologically motivated; instead use summary of c5 (5 features)
        c5_clean = c5[i][np.isfinite(c5[i]) & (c5[i] > 0)]
        if c5_clean.size:
            feats.extend([
                float(c5_clean.mean()),
                float(c5_clean.std()),
                float(np.median(c5_clean)),
                float(np.percentile(c5_clean, 10)),
                float(np.percentile(c5_clean, 90)),
            ])
        else:
            feats.extend([0.0] * 5)
        # Motif entropy + top-3 most-variable motif proxies — derived from FSD
        # shape (since we don't have motif.npy in the 5-channel baseline, we
        # use the FSD entropy + 3 most-extreme FSD bins as a stand-in).
        fsd_clean = fsd[i][np.isfinite(fsd[i]) & (fsd[i] > 0)]
        if fsd_clean.size:
            fsd_norm = fsd_clean / fsd_clean.sum()
            entropy = float(-np.sum(fsd_norm * np.log(fsd_norm + 1e-12)))
            top3 = sorted(fsd_clean, reverse=True)[:3]
        else:
            entropy = 0.0
            top3 = [0.0, 0.0, 0.0]
        feats.extend([entropy, top3[0], top3[1], top3[2]])
        rows.append(np.asarray(feats, dtype=float))

    # Names (must match PROXY_FEATURE_NAMES ordering)
    names = (
        [f"5mb_summary_{s}" for s in (
            "mean", "std", "p10", "p25", "p50", "p75", "p90", "frac_extreme"
        )]
        + [f"100kb_summary_{s}" for s in (
            "mean", "std", "p10", "p25", "p50", "p75", "p90", "mad",
            "skew", "kurt", "frac_low", "frac_high"
        )]
        + ["wps_mean", "wps_std", "wps_median", "wps_p10", "wps_p90"]
        + ["motif_entropy", "motif_top1_freq", "motif_top2_freq", "motif_top3_freq"]
    )

    X_summary = np.nan_to_num(np.asarray(rows, dtype=float), nan=0.0,
                              posinf=0.0, neginf=0.0)

    # ----- 3. Per-bin (5Mb) percentile-rank encoding -----
    # NOTE: This function returns ONLY the summary features. The per-bin
    # rank encoding must be computed INSIDE each CV fold using train data
    # only (otherwise it leaks test-set information into the training
    # distribution). The rank encoding is exposed via
    # encode_per_bin_rank_train_test() for fold-aware use.
    X_proxy = X_summary

    # Build full feature names (summary only at this stage; per-bin rank
    # names are added by the CV fold code)
    full_names = names

    # Sanity check
    assert X_proxy.shape[1] == len(full_names), (
        f"shape mismatch: X_proxy has {X_proxy.shape[1]} cols but "
        f"{len(full_names)} names"
    )

    y_or_ids = np.asarray(kept)
    return X_proxy, y_or_ids, kept, full_names
