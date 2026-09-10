#!/usr/bin/env python3
"""
Methylation baseline classifier — Phase 0/1.

Validates that methylation features produce a meaningful cancer-vs-healthy
classification signal. Designed for arrays (450K/EPIC) in Phase 0 and for
FinaleMe-imputed features in Phase 1.

Usage:
    from methylation_baseline import build_feature_matrix, train_lr_baseline

    X, y, feature_names = build_feature_matrix(beta_df, labels_df)
    auc = train_lr_baseline(X, y, n_seeds=5)
"""

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

__all__ = ["build_feature_matrix", "train_lr_baseline", "evaluate_cv"]


def build_feature_matrix(beta_df, labels_df, min_var=0.01, min_beta=0.0, max_beta=1.0):
    """Build (X, y, feature_names) from a β-value DataFrame and labels.

    Args:
        beta_df: DataFrame of shape (n_samples, n_probes) with β-values in [0,1].
        labels_df: DataFrame with 'sample_id' and 'label' (1=cancer, 0=healthy).
        min_var: drop probes with variance below this threshold.
        min_beta, max_beta: clip outliers (default: no clipping).

    Returns:
        X: ndarray of shape (n_samples, n_filtered_probes).
        y: ndarray of shape (n_samples,) with labels.
        feature_names: list of probe IDs that passed filtering.
    """
    # Align samples on sample_id
    merged = labels_df.merge(
        beta_df.reset_index().rename(columns={"index": "sample_id"}),
        on="sample_id",
        how="inner",
    )

    feature_cols = [c for c in merged.columns if c not in {"sample_id", "label"}]
    X = merged[feature_cols].values.astype(np.float32)
    y = merged["label"].values.astype(np.int32)

    # Clip outliers
    if min_beta is not None or max_beta is not None:
        X = np.clip(X, min_beta if min_beta is not None else X.min(),
                       max_beta if max_beta is not None else X.max())

    # Drop low-variance probes
    variances = X.var(axis=0)
    keep = variances >= min_var
    X = X[:, keep]
    feature_names = [f for f, k in zip(feature_cols, keep) if k]

    # Median-impute any NaN (per-feature)
    for j in range(X.shape[1]):
        col = X[:, j]
        if np.isnan(col).any():
            med = np.nanmedian(col)
            col[np.isnan(col)] = med
            X[:, j] = col

    return X, y, feature_names


def evaluate_cv(X, y, n_seeds=5, n_splits=5, C=1.0):
    """Run pooled OOF CV across multiple seeds. Returns per-seed AUCs + std.

    Uses stratified k-fold with per-fold StandardScaler + L2-LR.
    Per-fold scaling prevents train/test leakage.
    """
    from sklearn.metrics import roc_auc_score

    aucs = []
    for seed in range(n_seeds):
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        oof_pred = np.zeros(len(y), dtype=np.float32)

        for train_idx, test_idx in skf.split(X, y):
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X[train_idx])
            X_test = scaler.transform(X[test_idx])

            clf = LogisticRegression(C=C, solver="lbfgs", max_iter=1000, random_state=seed)
            clf.fit(X_train, y[train_idx])
            oof_pred[test_idx] = clf.predict_proba(X_test)[:, 1]

        auc = roc_auc_score(y, oof_pred)
        aucs.append(auc)

    return {
        "aucs": aucs,
        "auc_mean": float(np.mean(aucs)),
        "auc_std": float(np.std(aucs)),
        "n_samples": int(len(y)),
        "n_features": int(X.shape[1]),
        "n_seeds": n_seeds,
        "C": C,
    }


def train_lr_baseline(X, y, n_seeds=5, C=1.0):
    """Shorthand: build features → train LR baseline → return summary."""
    return evaluate_cv(X, y, n_seeds=n_seeds, C=C)
