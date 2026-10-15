"""Tests for the v0.11 methylation-array read pipeline."""

from __future__ import annotations

import pytest

from deepcatch_methylation.methylation_reader import (
    MethylationArrayInput,
    TISSUE_OF_ORIGIN_MARKERS,
    build_cfrna_persistence_input,
    extract_cfrna_signature,
    make_synthetic_liver_beta_values,
    make_synthetic_non_liver_beta_values,
)


# -----------------------------------------------------------------------------
# Public surface
# -----------------------------------------------------------------------------


def test_public_imports_resolve():
    assert callable(build_cfrna_persistence_input)
    assert callable(extract_cfrna_signature)
    assert isinstance(TISSUE_OF_ORIGIN_MARKERS, dict)
    assert "liver" in TISSUE_OF_ORIGIN_MARKERS
    assert "lung" in TISSUE_OF_ORIGIN_MARKERS
    assert "cns" in TISSUE_OF_ORIGIN_MARKERS


# -----------------------------------------------------------------------------
# extract_cfrna_signature
# -----------------------------------------------------------------------------


def test_extract_signature_liver_match():
    betas = make_synthetic_liver_beta_values("liver")
    tissue, mean_beta = extract_cfrna_signature(betas, "liver")
    assert tissue == "liver"
    assert mean_beta == pytest.approx(0.85, abs=1e-9)


def test_extract_signature_lung_match():
    betas = make_synthetic_liver_beta_values("lung")
    tissue, mean_beta = extract_cfrna_signature(betas, "lung")
    assert tissue == "lung"
    assert mean_beta == pytest.approx(0.85, abs=1e-9)


def test_extract_signature_mismatch_returns_unknown():
    betas = make_synthetic_liver_beta_values("liver")
    tissue, mean_beta = extract_cfrna_signature(betas, "lung")
    assert tissue == "unknown"
    assert mean_beta == 0.0


def test_extract_signature_too_few_markers_returns_unknown():
    betas = {"generic_marker_001": 0.9}
    tissue, mean_beta = extract_cfrna_signature(betas, "default")
    assert tissue == "unknown"


def test_extract_signature_unknown_tissue_falls_back_to_default():
    betas = {"generic_marker_001": 0.9, "generic_marker_002": 0.9}
    tissue, mean_beta = extract_cfrna_signature(betas, "unknown_tissue_xyz")
    assert tissue == "unknown"


# -----------------------------------------------------------------------------
# build_cfrna_persistence_input
# -----------------------------------------------------------------------------


def test_build_input_liver_match_medium_persistence():
    inp = MethylationArrayInput(
        sample_id="liver_001",
        tissue="liver",
        beta_values=make_synthetic_liver_beta_values("liver"),
        cfrna_concentration_pg_ml=500.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfrna_at_t0_pg_ml=10000.0,
    )
    pred = build_cfrna_persistence_input(inp)
    assert pred["persistence_category"] == "medium"
    assert pred["target_tissue"] == "liver"
    assert pred["feature_breakdown"]["concentration_pg_ml"] == 500.0
    assert pred["predicted_half_life_hours"] == 8.0


def test_build_input_tissue_mismatch_low_persistence():
    inp = MethylationArrayInput(
        sample_id="lung_001",
        tissue="lung",
        beta_values=make_synthetic_liver_beta_values("liver"),
        cfrna_concentration_pg_ml=200.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfrna_at_t0_pg_ml=10000.0,
    )
    pred = build_cfrna_persistence_input(inp)
    assert pred["persistence_category"] == "low"
    assert pred["target_tissue"] == "lung"


def test_build_input_mode1_fallback_derives_concentration_from_beta():
    inp = MethylationArrayInput(
        sample_id="liver_mode1",
        tissue="liver",
        beta_values=make_synthetic_liver_beta_values("liver"),
        cfrna_concentration_pg_ml=0.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfrna_at_t0_pg_ml=10000.0,
    )
    pred = build_cfrna_persistence_input(inp)
    assert pred["feature_breakdown"]["concentration_pg_ml"] == pytest.approx(8500.0, abs=1.0)


def test_build_input_mode2_uses_pre_quantified_concentration():
    inp = MethylationArrayInput(
        sample_id="liver_mode2",
        tissue="liver",
        beta_values={},
        cfrna_concentration_pg_ml=123.4,
        sampling_time_post_dose_hours=24.0,
        dose_cfrna_at_t0_pg_ml=10000.0,
    )
    pred = build_cfrna_persistence_input(inp)
    assert pred["feature_breakdown"]["concentration_pg_ml"] == 123.4


def test_build_input_zero_concentration_no_match_returns_low():
    inp = MethylationArrayInput(
        sample_id="empty",
        tissue="liver",
        beta_values={},
        cfrna_concentration_pg_ml=0.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfrna_at_t0_pg_ml=10000.0,
    )
    pred = build_cfrna_persistence_input(inp)
    assert pred["persistence_category"] == "low"


# -----------------------------------------------------------------------------
# End-to-end: build_cfrna_persistence_input → couple_cfrna_to_cdm
# -----------------------------------------------------------------------------


def test_end_to_end_methylation_to_cdm_observation():
    from deepcatch_methylation.cfrna_coupling import couple_cfrna_to_cdm

    inp = MethylationArrayInput(
        sample_id="e2e_001",
        tissue="liver",
        beta_values=make_synthetic_liver_beta_values("liver"),
        cfrna_concentration_pg_ml=500.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfrna_at_t0_pg_ml=10000.0,
    )
    pred = build_cfrna_persistence_input(inp)
    obs = couple_cfrna_to_cdm(
        prediction=pred,
        observed_cfRNA_concentration_pg_ml=400.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=10000.0,
    )
    # expected_residual at t=24h with half_life=8h is 0.5^3 = 0.125.
    # observed_residual = 400 / 10000 = 0.04.
    # ratio = 0.04 / 0.125 = 0.32 — below ANOMALY_RATIO_LOW=0.5 → anomaly.
    assert obs.anomaly_flag is True
    assert obs.dosing_recommendation in ("increase", "escalate")


# -----------------------------------------------------------------------------
# Synthetic data generators
# -----------------------------------------------------------------------------


def test_make_synthetic_liver_beta_values():
    betas = make_synthetic_liver_beta_values("liver")
    assert all(v == 0.85 for v in betas.values())
    for m in TISSUE_OF_ORIGIN_MARKERS["liver"]:
        assert m in betas


def test_make_synthetic_non_liver_beta_values_lung():
    betas = make_synthetic_non_liver_beta_values("lung")
    for m in TISSUE_OF_ORIGIN_MARKERS["liver"]:
        assert betas[m] == 0.30
    for m in TISSUE_OF_ORIGIN_MARKERS["lung"]:
        assert betas[m] == 0.85
