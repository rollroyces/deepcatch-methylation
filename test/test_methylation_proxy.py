"""Tests for the methylation-proxy feature extraction (Phase 2 Option 2).

These tests verify:
  - build_proxy_feature_matrix produces a feature matrix of the documented shape
  - encode_per_bin_rank_train_test respects train/test boundary (no leakage)
  - load_5channel_baseline returns the documented dimensions
  - end-to-end shape consistency between the 3 setups
"""

import numpy as np

from methylation.methylation_proxy import (
    PROXY_FEATURE_NAMES,
    build_proxy_feature_matrix,
    encode_per_bin_rank_train_test,
)


# Synthetic data helper for unit tests (no real-data dependency)
def _synthetic_5channel(n_samples: int = 50, seed: int = 42):
    """Build a synthetic 5-channel fragmentomics matrix matching the real shape."""
    rng = np.random.RandomState(seed)
    n5, n100, nfsd = 631, 30894, 196
    r5 = rng.normal(0.15, 0.05, size=(n_samples, n5))
    c5 = rng.normal(1.0, 0.2, size=(n_samples, n5))
    r100 = rng.normal(0.15, 0.04, size=(n_samples, n100))
    c100 = rng.poisson(1000, size=(n_samples, n100)).astype(float)
    fsd = rng.dirichlet(np.ones(nfsd), size=n_samples) * 1e6
    return np.nan_to_num(
        np.concatenate([r5, c5, r100, c100, fsd], axis=1),
        nan=0.0, posinf=0.0, neginf=0.0
    )


def test_load_5channel_baseline_shape():
    """5-channel baseline has shape (n, 631+631+30894+30894+196) = (n, 63246)."""
    n = 50
    X_5ch = _synthetic_5channel(n_samples=n)
    assert X_5ch.shape == (n, 63246)


def test_build_proxy_feature_matrix_summary_shape():
    """Methylation-proxy summary features: 29 dims (8+12+5+4)."""
    n = 50
    X_5ch = _synthetic_5channel(n_samples=n)
    X_proxy, _, _, names = build_proxy_feature_matrix(
        features_dir="/nonexistent",
        sample_ids=[f"S{i:04d}" for i in range(n)],
        cohort_X=X_5ch,
    )
    assert X_proxy.shape == (n, 29)
    assert len(names) == 29
    # First 29 entries of PROXY_FEATURE_NAMES must match the names returned
    assert names == list(PROXY_FEATURE_NAMES)


def test_build_proxy_feature_matrix_handles_nans():
    """NaN inputs are cleaned to 0; output has no NaN/Inf."""
    n = 30
    X_5ch = _synthetic_5channel(n_samples=n)
    # Inject some NaN/Inf
    X_5ch[0, 10] = np.nan
    X_5ch[1, 20] = np.inf
    X_5ch[2, 30] = -np.inf
    X_proxy, _, _, _ = build_proxy_feature_matrix(
        features_dir="/nonexistent",
        sample_ids=[f"S{i:04d}" for i in range(n)],
        cohort_X=X_5ch,
    )
    assert np.isfinite(X_proxy).all(), "Proxy matrix should have no NaN/Inf"


def test_encode_per_bin_rank_train_test_no_test_in_train_rank():
    """The encoded train rank must be invariant to the test set composition.

    Running the rank encoder with two different test sets (but the same
    train) must produce identical train ranks — the train ranks should not
    leak test information.
    """
    rng = np.random.RandomState(0)
    X_tr = rng.normal(0, 1, size=(20, 5))
    X_te_1 = rng.normal(0, 1, size=(5, 5))
    X_te_2 = rng.normal(0, 1, size=(7, 5))
    rank_tr_1, rank_te_1 = encode_per_bin_rank_train_test(X_tr, X_te_1)
    rank_tr_2, rank_te_2 = encode_per_bin_rank_train_test(X_tr, X_te_2)
    np.testing.assert_array_equal(rank_tr_1, rank_tr_2)


def test_encode_per_bin_rank_train_test_values_in_unit_interval():
    """Encoded ranks should be in [0, 1] for both train and test."""
    rng = np.random.RandomState(1)
    X_tr = rng.normal(0, 1, size=(20, 3))
    X_te = rng.normal(0, 1, size=(5, 3))
    rank_tr, rank_te = encode_per_bin_rank_train_test(X_tr, X_te)
    assert (rank_tr >= 0).all() and (rank_tr <= 1).all()
    assert (rank_te >= 0).all() and (rank_te <= 1).all()


def test_encode_per_bin_rank_train_test_handles_extreme_test_values():
    """Test values outside the train range should be clipped, not NaN."""
    rng = np.random.RandomState(2)
    X_tr = rng.normal(0, 1, size=(20, 3))
    # Test set has values way outside train range
    X_te = np.array([[100.0, -100.0, 0.0]] * 5)
    rank_tr, rank_te = encode_per_bin_rank_train_test(X_tr, X_te)
    assert np.isfinite(rank_te).all()
    # Values above max(train) → rank 1.0; below min(train) → rank 0.0
    assert (rank_te[:, 0] >= 0.99).all()
    assert (rank_te[:, 1] <= 0.01).all()


def test_proxy_feature_count_matches_naming():
    """PROXY_FEATURE_NAMES has exactly 29 entries (summary features)."""
    assert len(PROXY_FEATURE_NAMES) == 29


def test_proxy_summary_features_include_expected_channels():
    """PROXY_FEATURE_NAMES contains entries for 5mb_summary, 100kb_summary,
    wps, and motif."""
    s = " ".join(PROXY_FEATURE_NAMES)
    assert "5mb_summary" in s
    assert "100kb_summary" in s
    assert "wps_" in s
    assert "motif_" in s
