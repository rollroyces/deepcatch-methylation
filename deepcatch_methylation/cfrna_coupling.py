"""cfRNA → real-time CDM feature coupling (v0.10.0).

This module is the **read-side** counterpart to
``mrnavax_rare.safety.cfrna.predict_cfrna_persistence``. It consumes the
prediction shape emitted by the design side and lifts it into a
real-time chronic-disease-management (CDM) feature schema that
downstream clinical-decision tooling can act on.

Why a separate project?
-----------------------

``mrnavax-rare`` is the mRNA-therapeutic design side — it predicts how
long a designed CDS should persist in plasma. ``deepcatch-methylation``
is the cfRNA read side — it observes how much cfRNA is actually
circulating in the patient's plasma sample and flags anomalies against
the prediction. The two projects share NO runtime code (no shared
imports); they couple over a strict :class:`CfRNAPersistenceInput`
shape. This keeps the read side deployable independently of the design
side (e.g. for sites that already have mrnavax-rare outputs in a
research database but don't need the full design pipeline).

Feature schema
--------------

For every observed plasma cfRNA measurement the coupling produces a
:class:`RealTimeCDMObservation` carrying:

- ``observed_cfRNA_concentration_pg_ml`` — measured cfRNA in plasma.
- ``sampling_time_post_dose_hours`` — wall-clock hours since the dose.
- ``predicted_half_life_hours`` — passthrough from the input prediction.
- ``predicted_persistence_category`` — passthrough (``"high"`` /
  ``"medium"`` / ``"low"``).
- ``expected_residual_fraction`` — ``exp(-ln(2) * sampling_time /
  half_life)``, the idealised first-order half-life decay residual
  (equals 1.0 at t=0 and 0.5 at t=half_life).
- ``observed_residual_fraction`` — ``observed / dose_cfRNA_at_t0``.
- ``observed_vs_expected_ratio`` — ``observed_residual / expected_residual``.
- ``anomaly_flag`` — True when ratio is outside the
  ``[ANOMALY_RATIO_LOW, ANOMALY_RATIO_HIGH]`` interval.
- ``dosing_recommendation`` — :class:`DosingRecommendation` action.

The anomaly thresholds and dosing-recommendation rules are
**research-use defaults** in v0.10.0; they will need clinical
calibration before any production use.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping, TypedDict


# -----------------------------------------------------------------------------
# Input shape — strict typed-dict contract with mrnavax-rare.
# -----------------------------------------------------------------------------


class CfRNAPersistenceInput(TypedDict):
    """The minimal input shape deepcatch-methylation expects from mrnavax-rare.

    Mirrors :class:`mrnavax_rare.safety.cfrna.CfRNAPersistencePrediction.as_dict()`.
    Defined as a :class:`TypedDict` (not a dataclass) so this project
    does not need to import mrnavax-rare at runtime — the contract is
    structural, not nominal.
    """

    predicted_half_life_hours: float
    persistence_category: str  # "high" | "medium" | "low"
    target_tissue: str
    codon_optimization: str
    feature_breakdown: Mapping[str, float]


# -----------------------------------------------------------------------------
# Constants — honest-disclosure research-use defaults.
# -----------------------------------------------------------------------------


#: Lower bound of the anomaly-flag interval. An ``observed_vs_expected_ratio``
#: below this threshold flags the observation as "faster clearance than
#: predicted" (potential under-dosing, accelerated clearance, or assay
#: artefact).
#:
#: v0.10.0 default: 0.5 (the patient retains <=50% of the expected residual).
#: Cite as research-use; clinical calibration required.
ANOMALY_RATIO_LOW: float = 0.5

#: Upper bound of the anomaly-flag interval. An ``observed_vs_expected_ratio``
#: above this threshold flags the observation as "slower clearance than
#: predicted" (potential accumulation, saturation, or measurement
#: artefact).
#:
#: v0.10.0 default: 2.0 (the patient retains >=200% of the expected residual).
#: Cite as research-use; clinical calibration required.
ANOMALY_RATIO_HIGH: float = 2.0

#: Minimum predicted half-life (hours) used to guard the residual-fraction
#: computation. Below this we treat the predicted half-life as effectively
#: instantaneous and the dosing recommendation fires as ``"hold"``.
MIN_HALFLIFE_HOURS: float = 1.0


# -----------------------------------------------------------------------------
# Output schema — real-time CDM observation.
# -----------------------------------------------------------------------------


#: Canonical dosing-recommendation literals. Kept as a plain string union
#: for JSON-serialisability.
DosingRecommendation = str  # "maintain" | "increase" | "decrease" | "hold" | "escalate"


@dataclass(frozen=True)
class RealTimeCDMObservation:
    """One cfRNA plasma observation lifted into the real-time CDM feature schema."""

    observed_cfRNA_concentration_pg_ml: float
    sampling_time_post_dose_hours: float
    dose_cfRNA_at_t0_pg_ml: float
    predicted_half_life_hours: float
    predicted_persistence_category: str
    expected_residual_fraction: float
    observed_residual_fraction: float
    observed_vs_expected_ratio: float
    anomaly_flag: bool
    dosing_recommendation: DosingRecommendation
    notes: tuple[str, ...] = field(default_factory=tuple)


__all__ = [
    "ANOMALY_RATIO_HIGH",
    "ANOMALY_RATIO_LOW",
    "CfRNAPersistenceInput",
    "DosingRecommendation",
    "MIN_HALFLIFE_HOURS",
    "RealTimeCDMObservation",
    "couple_cfrna_to_cdm",
]


# -----------------------------------------------------------------------------
# Coupling function
# -----------------------------------------------------------------------------


def _validate_input(prediction: CfRNAPersistenceInput) -> tuple[float, str]:
    """Validate the input shape; return ``(half_life, category)``."""
    half_life = float(prediction["predicted_half_life_hours"])
    category = str(prediction["persistence_category"]).strip().lower()
    if category not in {"high", "medium", "low"}:
        raise ValueError(
            f"Unknown persistence_category: {category!r}; expected one of "
            f"'high', 'medium', 'low'."
        )
    if half_life <= 0.0:
        raise ValueError(
            f"predicted_half_life_hours must be positive; got {half_life!r}."
        )
    if not math.isfinite(half_life):
        raise ValueError(
            f"predicted_half_life_hours must be finite; got {half_life!r}."
        )
    return half_life, category


def _compute_residuals(
    observed_conc: float,
    dose_at_t0: float,
    sampling_time_h: float,
    half_life_h: float,
) -> tuple[float, float, float]:
    """Compute ``(expected, observed, ratio)``.

    Uses the canonical first-order half-life decay: at t == half_life
    the residual is exactly 0.5. ``exp(-ln(2) * t / t_half_life)`` is
    equivalent to ``2 ** (-t / t_half_life)`` and matches the
    half-life definition.
    """
    if dose_at_t0 <= 0.0:
        raise ValueError(
            f"dose_cfRNA_at_t0_pg_ml must be positive; got {dose_at_t0!r}."
        )
    if sampling_time_h < 0.0:
        raise ValueError(
            f"sampling_time_post_dose_hours must be non-negative; "
            f"got {sampling_time_h!r}."
        )
    effective_half_life = max(half_life_h, MIN_HALFLIFE_HOURS)
    expected_residual = math.exp(-math.log(2) * sampling_time_h / effective_half_life)
    observed_residual = observed_conc / dose_at_t0
    if expected_residual <= 0.0:
        # Should not happen with finite sampling_time and effective_half_life > 0,
        # but guard anyway to avoid 0/0 in the ratio.
        ratio = float("inf")
    else:
        ratio = observed_residual / expected_residual
    return expected_residual, observed_residual, ratio


def _recommend(
    ratio: float,
    anomaly: bool,
    half_life_h: float,
    category: str,
) -> tuple[DosingRecommendation, tuple[str, ...]]:
    """Map ``(ratio, anomaly, half_life, category)`` to ``(rec, notes)``."""
    notes: list[str] = []

    # Hold recommendation fires before ratio-based rules — vanishingly
    # short half-life predictions are a compound signal of either upstream
    # numerical failure or a CDS that won't persist.
    if half_life_h < MIN_HALFLIFE_HOURS:
        notes.append(
            f"predicted half-life {half_life_h} h is below MIN_HALFLIFE_HOURS "
            f"= {MIN_HALFLIFE_HOURS} h — dosing held pending re-prediction."
        )
        return ("hold", tuple(notes))

    if anomaly:
        if ratio < ANOMALY_RATIO_LOW:
            notes.append(
                f"observed_vs_expected_ratio = {ratio:.3f} below "
                f"ANOMALY_RATIO_LOW = {ANOMALY_RATIO_LOW} — faster clearance "
                f"than predicted; potential under-dosing."
            )
            if category == "low":
                notes.append(
                    "compound signal: predicted persistence is also 'low' — "
                    "recommend escalating to clinical review."
                )
                return ("escalate", tuple(notes))
            return ("increase", tuple(notes))
        if ratio > ANOMALY_RATIO_HIGH:
            notes.append(
                f"observed_vs_expected_ratio = {ratio:.3f} above "
                f"ANOMALY_RATIO_HIGH = {ANOMALY_RATIO_HIGH} — slower clearance "
                f"than predicted; potential accumulation / saturation."
            )
            if category == "low":
                notes.append(
                    "compound signal: predicted persistence is 'low' but "
                    "actual is much higher — review for compliance / dosing "
                    "interval errors."
                )
                return ("escalate", tuple(notes))
            return ("decrease", tuple(notes))

    return ("maintain", tuple(notes))


def couple_cfrna_to_cdm(
    prediction: CfRNAPersistenceInput,
    observed_cfRNA_concentration_pg_ml: float,
    sampling_time_post_dose_hours: float,
    dose_cfRNA_at_t0_pg_ml: float,
) -> RealTimeCDMObservation:
    """Couple a mrnavax-rare cfRNA prediction to a real-time CDM observation.

    Parameters
    ----------
    prediction:
        The input shape (mirrors
        :class:`mrnavax_rare.safety.cfrna.CfRNAPersistencePrediction.as_dict()`).
        Required keys: ``predicted_half_life_hours``, ``persistence_category``.
    observed_cfRNA_concentration_pg_ml:
        The measured cfRNA concentration in the patient's plasma sample.
    sampling_time_post_dose_hours:
        Wall-clock hours between the most recent dose and the sample draw.
    dose_cfRNA_at_t0_pg_ml:
        Reference cfRNA concentration immediately post-dose — the
        patient's per-dose baseline. Used to compute the residual
        fraction.

    Returns
    -------
    RealTimeCDMObservation
        The CDM feature schema row.

    Raises
    ------
    ValueError
        For invalid prediction shape (non-positive or non-finite
        half-life, unknown persistence_category), non-positive
        ``dose_cfRNA_at_t0_pg_ml``, or negative
        ``sampling_time_post_dose_hours``.
    """
    half_life_h, category = _validate_input(prediction)
    expected_residual, observed_residual, ratio = _compute_residuals(
        observed_conc=observed_cfRNA_concentration_pg_ml,
        dose_at_t0=dose_cfRNA_at_t0_pg_ml,
        sampling_time_h=sampling_time_post_dose_hours,
        half_life_h=half_life_h,
    )
    anomaly = not (ANOMALY_RATIO_LOW <= ratio <= ANOMALY_RATIO_HIGH)
    rec, notes = _recommend(ratio, anomaly, half_life_h, category)
    return RealTimeCDMObservation(
        observed_cfRNA_concentration_pg_ml=float(observed_cfRNA_concentration_pg_ml),
        sampling_time_post_dose_hours=float(sampling_time_post_dose_hours),
        dose_cfRNA_at_t0_pg_ml=float(dose_cfRNA_at_t0_pg_ml),
        predicted_half_life_hours=half_life_h,
        predicted_persistence_category=category,
        expected_residual_fraction=expected_residual,
        observed_residual_fraction=observed_residual,
        observed_vs_expected_ratio=ratio,
        anomaly_flag=anomaly,
        dosing_recommendation=rec,
        notes=notes,
    )
