"""Tests for the v0.10.0 deepcatch-methylation cfRNA → CDM coupling.

Coverage
--------

- Public surface (imports, public types, public functions).
- Happy path: predicted vs observed in range → ``"maintain"``, no anomaly.
- Anomaly detection (faster clearance → ``"increase"``; slower clearance
  → ``"decrease"``).
- Compound signal (anomaly + ``"low"`` category → ``"escalate"``).
- Hold rule (predicted half-life below ``MIN_HALFLIFE_HOURS`` → ``"hold"``).
- Edge cases: zero observation, very large observation, sampling at t=0,
  sampling at t=half_life.
- Negative paths: non-positive half-life, non-positive dose_at_t0,
  negative sampling time, unknown persistence_category.
- Output schema invariants: frozen dataclass, JSON-friendly primitives.
"""

from __future__ import annotations

import math

import pytest

from deepcatch_methylation.cfrna_coupling import (
    ANOMALY_RATIO_HIGH,
    ANOMALY_RATIO_LOW,
    CfRNAPersistenceInput,
    MIN_HALFLIFE_HOURS,
    RealTimeCDMObservation,
    couple_cfrna_to_cdm,
)


#: Standard test prediction (mirrors a v0.7.4 mRNA-3927 IV-liver prediction).
STANDARD_PREDICTION: CfRNAPersistenceInput = {
    "predicted_half_life_hours": 54.0,
    "persistence_category": "medium",
    "target_tissue": "liver",
    "codon_optimization": "balanced",
    "feature_breakdown": {"cai_factor": 0.85, "gc_factor": 0.92, "rare_factor": 0.78},
}


# -----------------------------------------------------------------------------
# Public surface
# -----------------------------------------------------------------------------


def test_public_imports_resolve():
    """All public symbols importable from the module."""
    assert callable(couple_cfrna_to_cdm)
    assert isinstance(ANOMALY_RATIO_LOW, float)
    assert isinstance(ANOMALY_RATIO_HIGH, float)
    assert MIN_HALFLIFE_HOURS > 0
    assert ANOMALY_RATIO_LOW < ANOMALY_RATIO_HIGH


def test_real_time_cdm_observation_is_frozen_dataclass():
    """RealTimeCDMObservation is a frozen dataclass — mutating raises."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=100.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    with pytest.raises(Exception):
        obs.dosing_recommendation = "increase"  # type: ignore[misc]


# -----------------------------------------------------------------------------
# Happy paths
# -----------------------------------------------------------------------------


def test_in_range_maintain_no_anomaly():
    """At sampling_time = 0, observed_residual == 1.0, expected_residual == 1.0."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=1000.0,
        sampling_time_post_dose_hours=0.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.expected_residual_fraction == pytest.approx(1.0)
    assert obs.observed_residual_fraction == pytest.approx(1.0)
    assert obs.observed_vs_expected_ratio == pytest.approx(1.0)
    assert obs.anomaly_flag is False
    assert obs.dosing_recommendation == "maintain"
    assert obs.notes == ()


def test_halflife_sampling_yields_half_residual():
    """At sampling_time = half_life, expected_residual == 0.5 (within rounding)."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=500.0,  # observed also = half
        sampling_time_post_dose_hours=54.0,  # exactly half_life
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.expected_residual_fraction == pytest.approx(0.5, abs=1e-9)
    assert obs.observed_residual_fraction == pytest.approx(0.5)
    assert obs.observed_vs_expected_ratio == pytest.approx(1.0)
    assert obs.anomaly_flag is False
    assert obs.dosing_recommendation == "maintain"


def test_passthrough_predicted_half_life_and_category():
    """Half-life and category are passed through unchanged."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=1000.0,
        sampling_time_post_dose_hours=0.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.predicted_half_life_hours == 54.0
    assert obs.predicted_persistence_category == "medium"


# -----------------------------------------------------------------------------
# Anomaly detection
# -----------------------------------------------------------------------------


def test_faster_clearance_triggers_increase():
    """observed_residual much lower than expected → 'increase' recommendation."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=100.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.observed_vs_expected_ratio < ANOMALY_RATIO_LOW
    assert obs.anomaly_flag is True
    assert obs.dosing_recommendation == "increase"
    assert any("faster clearance" in n for n in obs.notes)


def test_slower_clearance_triggers_decrease():
    """observed_residual much higher than expected → 'decrease' recommendation."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=2000.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.observed_vs_expected_ratio > ANOMALY_RATIO_HIGH
    assert obs.anomaly_flag is True
    assert obs.dosing_recommendation == "decrease"
    assert any("slower clearance" in n for n in obs.notes)


def test_anomaly_low_with_low_category_escalates():
    """Compound signal: faster-clearance anomaly + 'low' category → 'escalate'."""
    low_prediction: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    low_prediction["persistence_category"] = "low"
    obs = couple_cfrna_to_cdm(
        prediction=low_prediction,
        observed_cfRNA_concentration_pg_ml=100.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.dosing_recommendation == "escalate"
    assert any("compound signal" in n for n in obs.notes)


def test_anomaly_high_with_low_category_escalates():
    """Compound signal: slower-clearance anomaly + 'low' category → 'escalate'."""
    low_prediction: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    low_prediction["persistence_category"] = "low"
    obs = couple_cfrna_to_cdm(
        prediction=low_prediction,
        observed_cfRNA_concentration_pg_ml=2000.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.dosing_recommendation == "escalate"


def test_anomaly_at_exact_boundary_is_not_flagged():
    """ratio exactly at ANOMALY_RATIO_LOW or _HIGH is inclusive → not flagged."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=ANOMALY_RATIO_LOW * 1000.0,
        sampling_time_post_dose_hours=0.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.observed_vs_expected_ratio == pytest.approx(ANOMALY_RATIO_LOW)
    assert obs.anomaly_flag is False  # inclusive boundary
    assert obs.dosing_recommendation == "maintain"


# -----------------------------------------------------------------------------
# Hold rule
# -----------------------------------------------------------------------------


def test_hold_when_half_life_below_minimum():
    """Predicted half-life < MIN_HALFLIFE_HOURS → 'hold' regardless of ratio."""
    short_prediction: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    short_prediction["predicted_half_life_hours"] = 0.5  # below MIN = 1.0
    obs = couple_cfrna_to_cdm(
        prediction=short_prediction,
        observed_cfRNA_concentration_pg_ml=500.0,
        sampling_time_post_dose_hours=1.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.dosing_recommendation == "hold"
    assert any("below MIN_HALFLIFE_HOURS" in n for n in obs.notes)


def test_hold_fires_before_anomaly_rules():
    """Hold overrides faster-clearance / slower-clearance even if both would fire."""
    short_prediction: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    short_prediction["predicted_half_life_hours"] = 0.1
    short_prediction["persistence_category"] = "low"
    obs = couple_cfrna_to_cdm(
        prediction=short_prediction,
        observed_cfRNA_concentration_pg_ml=50.0,
        sampling_time_post_dose_hours=2.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.dosing_recommendation == "hold"


# -----------------------------------------------------------------------------
# Edge cases
# -----------------------------------------------------------------------------


def test_sampling_at_zero_yields_unit_residual():
    """At t=0, expected_residual == 1.0 regardless of half-life."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=1000.0,
        sampling_time_post_dose_hours=0.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.expected_residual_fraction == pytest.approx(1.0)


def test_very_long_sampling_yields_tiny_residual():
    """At t = 10 × half_life, expected_residual = 2^(-10) ≈ 9.77e-4."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=0.0,
        sampling_time_post_dose_hours=540.0,  # 10 × half_life
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    # 2^(-10) ≈ 9.77e-4; verify the formula is correct.
    assert obs.expected_residual_fraction == pytest.approx(2.0 ** -10, rel=1e-9)
    assert obs.expected_residual_fraction < 1e-3


def test_zero_observation_is_anomaly():
    """observed_cfRNA == 0 yields observed_residual = 0 → ratio = 0 → anomaly."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=0.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert obs.observed_residual_fraction == 0.0
    assert obs.anomaly_flag is True
    assert obs.dosing_recommendation in {"increase", "escalate"}


# -----------------------------------------------------------------------------
# Negative paths
# -----------------------------------------------------------------------------


def test_non_positive_half_life_raises():
    """Zero or negative predicted half-life raises ValueError."""
    bad: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    bad["predicted_half_life_hours"] = 0.0
    with pytest.raises(ValueError, match="half_life_hours"):
        couple_cfrna_to_cdm(
            prediction=bad,
            observed_cfRNA_concentration_pg_ml=100.0,
            sampling_time_post_dose_hours=24.0,
            dose_cfRNA_at_t0_pg_ml=1000.0,
        )


def test_negative_half_life_raises():
    bad: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    bad["predicted_half_life_hours"] = -5.0
    with pytest.raises(ValueError):
        couple_cfrna_to_cdm(
            prediction=bad,
            observed_cfRNA_concentration_pg_ml=100.0,
            sampling_time_post_dose_hours=24.0,
            dose_cfRNA_at_t0_pg_ml=1000.0,
        )


def test_non_finite_half_life_raises():
    bad: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    bad["predicted_half_life_hours"] = math.inf
    with pytest.raises(ValueError, match="finite"):
        couple_cfrna_to_cdm(
            prediction=bad,
            observed_cfRNA_concentration_pg_ml=100.0,
            sampling_time_post_dose_hours=24.0,
            dose_cfRNA_at_t0_pg_ml=1000.0,
        )


def test_unknown_category_raises():
    bad: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
    bad["persistence_category"] = "extreme"
    with pytest.raises(ValueError, match="Unknown persistence_category"):
        couple_cfrna_to_cdm(
            prediction=bad,
            observed_cfRNA_concentration_pg_ml=100.0,
            sampling_time_post_dose_hours=24.0,
            dose_cfRNA_at_t0_pg_ml=1000.0,
        )


def test_zero_dose_at_t0_raises():
    """dose_cfRNA_at_t0 must be positive — observed_residual is undefined."""
    with pytest.raises(ValueError, match="dose_cfRNA_at_t0_pg_ml"):
        couple_cfrna_to_cdm(
            prediction=STANDARD_PREDICTION,
            observed_cfRNA_concentration_pg_ml=100.0,
            sampling_time_post_dose_hours=24.0,
            dose_cfRNA_at_t0_pg_ml=0.0,
        )


def test_negative_sampling_time_raises():
    with pytest.raises(ValueError, match="sampling_time_post_dose_hours"):
        couple_cfrna_to_cdm(
            prediction=STANDARD_PREDICTION,
            observed_cfRNA_concentration_pg_ml=100.0,
            sampling_time_post_dose_hours=-1.0,
            dose_cfRNA_at_t0_pg_ml=1000.0,
        )


# -----------------------------------------------------------------------------
# Schema invariants
# -----------------------------------------------------------------------------


def test_output_is_json_friendly():
    """All fields are primitive numbers / bools / str / tuple of str."""
    obs = couple_cfrna_to_cdm(
        prediction=STANDARD_PREDICTION,
        observed_cfRNA_concentration_pg_ml=100.0,
        sampling_time_post_dose_hours=24.0,
        dose_cfRNA_at_t0_pg_ml=1000.0,
    )
    assert isinstance(obs.observed_cfRNA_concentration_pg_ml, float)
    assert isinstance(obs.sampling_time_post_dose_hours, float)
    assert isinstance(obs.dose_cfRNA_at_t0_pg_ml, float)
    assert isinstance(obs.predicted_half_life_hours, float)
    assert isinstance(obs.predicted_persistence_category, str)
    assert isinstance(obs.expected_residual_fraction, float)
    assert isinstance(obs.observed_residual_fraction, float)
    assert isinstance(obs.observed_vs_expected_ratio, float)
    assert isinstance(obs.anomaly_flag, bool)
    assert isinstance(obs.dosing_recommendation, str)
    assert isinstance(obs.notes, tuple)
    assert all(isinstance(n, str) for n in obs.notes)


def test_dosing_recommendation_is_in_canonical_set():
    """dosing_recommendation is one of the canonical literals."""
    valid_recs = {"maintain", "increase", "decrease", "hold", "escalate"}
    for category in ("high", "medium", "low"):
        for ratio_pick in ("fast", "match", "slow"):
            pred: CfRNAPersistenceInput = dict(STANDARD_PREDICTION)
            pred["persistence_category"] = category
            if ratio_pick == "fast":
                obs_pg, time_h = 100.0, 24.0
            elif ratio_pick == "slow":
                obs_pg, time_h = 2000.0, 24.0
            else:
                obs_pg, time_h = 641.0, 24.0
            obs = couple_cfrna_to_cdm(
                prediction=pred,
                observed_cfRNA_concentration_pg_ml=obs_pg,
                sampling_time_post_dose_hours=time_h,
                dose_cfRNA_at_t0_pg_ml=1000.0,
            )
            assert obs.dosing_recommendation in valid_recs
