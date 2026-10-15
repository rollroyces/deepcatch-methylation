"""Methylation-array → cfRNA persistence input read pipeline (v0.11).

This module is the **read-side** complement to the v0.10.0 cfRNA coupling
(``deepcatch_methylation.cfrna_coupling.couple_cfrna_to_cdm``). It defines
a standardized read pipeline that ingests methylation array data
(Illumina EPIC 850K or processed beta values) and produces a
:class:`~deepcatch_methylation.cfrna_coupling.CfRNAPersistenceInput`
that the coupling can consume.

Architecture
------------

::

    methylation array (IDAT / processed beta values)
        → ``MethylationArrayInput`` (TypedDict)
            → ``extract_cfrna_signature()`` (tissue-of-origin identification)
                → ``build_cfrna_persistence_input()`` (normalizer)
                    → ``CfRNAPersistenceInput`` (TypedDict)
                        → ``couple_cfrna_to_cdm()`` (anomaly + dosing)
                            → ``RealTimeCDMObservation`` (CDM feature schema)

Honest disclosure
-----------------

The v0.11 read pipeline is a **structural scaffold** that defines the
input contract for cfRNA persistence coupling. It does NOT call a
real methylation-array analysis pipeline (e.g. minfi, methylKit, SeSAMe)
and does NOT extract tissue-of-origin markers from raw IDAT files.
Real production deployment would need:

1. **Real methylation-array processing** — IDAT file parsing via
   ``minfi`` (Illumina EPIC 850K), ``methylKit``, or ``SeSAMe``; raw
   signal → beta values / M values.
2. **Tissue-of-origin marker panel** — cfRNA carries methylation
   signatures from its tissue of origin (e.g. liver-specific
   methylation markers for an mRNA-3927 / mRNA-3139 patient with
   PCCA / OTC; lung-specific markers for a CF patient on 4D-710;
   brain-specific markers for a Dravet patient on STK-001).
3. **m5C detection** — 5-methylcytosine in cfRNA can be detected
   via bisulfite sequencing (BS-seq) or Nanopore direct RNA sequencing.
   The tissue-of-origin signal lives in the m5C pattern.
4. **cfRNA concentration quantification** — spike-in calibrated
   RT-qPCR or digital-droplet PCR (ddPCR) on plasma samples.

v0.11 ships the **input contract** (TypedDict shapes, normalizer
function, fake data generators for tests) so the coupling has a real
read-side to consume. The wet-lab read pipeline is a v0.12+ target.
"""

from __future__ import annotations

import math
from typing import Mapping, TypedDict


# -----------------------------------------------------------------------------
# Standardized input shapes
# -----------------------------------------------------------------------------


#: Canonical tissue-of-origin marker panel for the v0.11 cfRNA read
#: pipeline. The keys are tissue tokens (matching
#: ``mrnavax_rare.aav.AAV_CAPSID_TROPISM`` and
#: ``mrnavax_rare.safety.cfrna.TISSUE_FACTORS``); the values are
#: lists of methylation marker names (placeholder — real production
#: uses Illumina EPIC 850K probe IDs like ``cg00000029``).
#:
#: Honest v0.11 disclosure: the marker names below are **structural
#: placeholders**, not real probe IDs. A real production deployment
#: would use the validated tissue-of-origin marker panels from the
#: cfRNA / cell-free DNA literature (e.g. Guo 2017 *Nat Commun* 8:1822
#: for liver / lung / brain cfDNA methylation markers; Moss 2018 *Nat
#: Commun* 9:5068 for placenta cfRNA markers).
TISSUE_OF_ORIGIN_MARKERS: dict[str, tuple[str, ...]] = {
    "liver": ("liver_marker_001", "liver_marker_002", "liver_marker_003"),
    "lung": ("lung_marker_001", "lung_marker_002", "lung_marker_003"),
    "cns": ("cns_marker_001", "cns_marker_002", "cns_marker_003"),
    "muscle": ("muscle_marker_001", "muscle_marker_002"),
    "cardiac": ("cardiac_marker_001", "cardiac_marker_002"),
    "default": ("generic_marker_001",),
}


class MethylationArrayInput(TypedDict):
    """The standard input shape for the v0.11 methylation read pipeline."""

    sample_id: str
    tissue: str
    beta_values: Mapping[str, float]
    cfrna_concentration_pg_ml: float
    sampling_time_post_dose_hours: float
    dose_cfrna_at_t0_pg_ml: float


#: v0.11 — minimum number of tissue-specific markers required to call a
#: tissue-of-origin match with confidence.
_MIN_MARKERS_FOR_TISSUE_CALL: int = 2

#: v0.11 — beta-value threshold above which a marker is considered
#: "methylated" in the calling tissue. The standard threshold for
#: tissue-of-origin calls is 0.7 (cfDNA methylation papers); cfRNA
#: methylation is generally lower than cfDNA so 0.5 is the practical
#: threshold.
_TISSUE_BETA_THRESHOLD: float = 0.5


# -----------------------------------------------------------------------------
# Normalizer
# -----------------------------------------------------------------------------


def extract_cfrna_signature(
    beta_values: Mapping[str, float],
    tissue: str,
) -> tuple[str, float]:
    """Identify the tissue-of-origin from methylation beta values."""
    markers = TISSUE_OF_ORIGIN_MARKERS.get(tissue, TISSUE_OF_ORIGIN_MARKERS["default"])
    hit_count = 0
    beta_sum = 0.0
    for marker in markers:
        if marker in beta_values:
            beta = float(beta_values[marker])
            if beta >= _TISSUE_BETA_THRESHOLD:
                hit_count += 1
                beta_sum += beta
    if hit_count >= _MIN_MARKERS_FOR_TISSUE_CALL:
        return tissue, beta_sum / hit_count
    return "unknown", 0.0


def build_cfrna_persistence_input(
    methylation_input: MethylationArrayInput,
) -> "CfRNAPersistenceInput":
    """Build a :class:`CfRNAPersistenceInput` from methylation data."""
    # Defer the import to avoid a circular import at module load.
    from .cfrna_coupling import CfRNAPersistenceInput  # noqa: PLC0415

    tissue = methylation_input["tissue"].strip().lower()
    beta_values = methylation_input.get("beta_values", {}) or {}
    cfrna_conc = float(methylation_input.get("cfrna_concentration_pg_ml", 0.0))

    called_tissue, mean_beta = extract_cfrna_signature(beta_values, tissue)

    # Mode 1 fallback: if no pre-quantified concentration, derive a
    # structural stub from the mean beta.
    if cfrna_conc <= 0.0:
        cfrna_conc = mean_beta * float(methylation_input.get("dose_cfrna_at_t0_pg_ml", 0.0))

    # Persistence category (structural)
    if called_tissue == tissue and cfrna_conc > 0.0:
        persistence_category = "medium"
    elif cfrna_conc > 0.0:
        persistence_category = "low"  # detected but tissue mismatch
    else:
        persistence_category = "low"  # below detection

    # Half-life (structural default 8 h)
    half_life_h = 8.0

    return CfRNAPersistenceInput(
        predicted_half_life_hours=half_life_h,
        persistence_category=persistence_category,
        target_tissue=tissue,
        codon_optimization="balanced",
        feature_breakdown={
            "mean_beta": mean_beta,
            "called_tissue": 1.0 if called_tissue == tissue else 0.0,
            "concentration_pg_ml": cfrna_conc,
        },
    )


# -----------------------------------------------------------------------------
# Synthetic data generators (for tests)
# -----------------------------------------------------------------------------


def make_synthetic_liver_beta_values(tissue: str = "liver") -> dict[str, float]:
    """Generate synthetic beta values for a tissue-derived cfRNA sample."""
    markers = TISSUE_OF_ORIGIN_MARKERS.get(tissue, TISSUE_OF_ORIGIN_MARKERS["default"])
    return {m: 0.85 for m in markers}


def make_synthetic_non_liver_beta_values(tissue: str = "lung") -> dict[str, float]:
    """Generate synthetic beta values for a non-liver cfRNA sample.

    The liver markers are below the tissue-call threshold (≤ 0.5)
    but the tissue-specific markers are above (≥ 0.5) — exercises
    the "tissue mismatch" path in ``extract_cfrna_signature``.
    """
    liver_markers = TISSUE_OF_ORIGIN_MARKERS["liver"]
    tissue_markers = TISSUE_OF_ORIGIN_MARKERS.get(tissue, TISSUE_OF_ORIGIN_MARKERS["default"])
    return {
        **{m: 0.30 for m in liver_markers},  # below threshold
        **{m: 0.85 for m in tissue_markers},  # above threshold
    }


__all__ = [
    "MethylationArrayInput",
    "TISSUE_OF_ORIGIN_MARKERS",
    "build_cfrna_persistence_input",
    "extract_cfrna_signature",
    "make_synthetic_liver_beta_values",
    "make_synthetic_non_liver_beta_values",
]
