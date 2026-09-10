"""Tests for the methylation baseline pipeline (Phase 0)."""

import numpy as np
import pandas as pd

from methylation.methylation_baseline import build_feature_matrix, evaluate_cv


def _make_synthetic(n_samples=100, n_probes=500, n_cancer=50, seed=42):
    """Generate synthetic β-values with cancer/healthy separation."""
    rng = np.random.RandomState(seed)
    # Cancer samples have a subset of probes hypermethylated
    beta = rng.beta(2, 2, size=(n_samples, n_probes)).astype(np.float32)
    cancer_mask = np.zeros(n_samples, dtype=bool)
    cancer_mask[:n_cancer] = True
    # First 50 probes: shift cancer up
    beta[cancer_mask, :50] += rng.uniform(0.3, 0.5, size=(n_cancer, 50))
    beta = np.clip(beta, 0.0, 1.0)
    sample_ids = [f"S{i:04d}" for i in range(n_samples)]
    beta_df = pd.DataFrame(beta, index=sample_ids, columns=[f"cg{i:07d}" for i in range(n_probes)])
    labels_df = pd.DataFrame({"sample_id": sample_ids,
                              "label": cancer_mask.astype(int)})
    return beta_df, labels_df


def test_build_feature_matrix_shape():
    beta_df, labels_df = _make_synthetic()
    X, y, names = build_feature_matrix(beta_df, labels_df)
    assert X.shape[0] == 100
    assert y.shape[0] == 100
    assert len(names) == X.shape[1]
    assert np.all((y == 0) | (y == 1))
    assert 0 <= X.min() <= X.max() <= 1.0


def test_build_feature_matrix_filters_low_variance():
    beta_df, labels_df = _make_synthetic(n_probes=200)
    # Drop low-variance probes — synthetic data has nonzero variance everywhere
    X, y, names = build_feature_matrix(beta_df, labels_df, min_var=1e-6)
    assert X.shape[1] > 100  # Should keep most probes


def test_evaluate_cv_synthetic_signal_above_chance():
    """With strong cancer signal in synthetic data, AUC should be ~1.0."""
    beta_df, labels_df = _make_synthetic(n_samples=100, n_probes=500, n_cancer=50)
    X, y, _ = build_feature_matrix(beta_df, labels_df)
    result = evaluate_cv(X, y, n_seeds=2, n_splits=5, C=1.0)
    assert result["auc_mean"] > 0.9, f"Synthetic signal too weak: {result}"
    assert result["n_samples"] == 100
    assert result["n_features"] == X.shape[1]


def test_evaluate_cv_random_labels_is_at_chance():
    """With no signal (random labels), AUC should be ~0.5."""
    rng = np.random.RandomState(123)
    X = rng.beta(2, 2, size=(100, 50)).astype(np.float32)
    y = rng.randint(0, 2, size=100).astype(np.int32)
    result = evaluate_cv(X, y, n_seeds=2, n_splits=5, C=1.0)
    # Wide tolerance: random labels give AUC in [0.3, 0.7]
    assert 0.3 < result["auc_mean"] < 0.7, f"Expected chance-level AUC, got {result['auc_mean']}"


def test_evaluate_cv_handles_perfect_separation():
    """If features perfectly separate classes, AUC should be ~1.0."""
    X = np.zeros((60, 5), dtype=np.float32)
    y = np.zeros(60, dtype=np.int32)
    X[:30, 0] = 1.0
    y[:30] = 1
    result = evaluate_cv(X, y, n_seeds=1, n_splits=5, C=1.0)
    assert result["auc_mean"] > 0.95
