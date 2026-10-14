"""deepcatch-methylation — real-time cfRNA read-side companion to mrnavax-rare.

This is the v0.10.0 initial scaffold. The deepcatch-methylation project
consumes the cfRNA persistence predictions emitted by
``mrnavax_rare.safety.cfrna.predict_cfrna_persistence`` and lifts them
into a real-time chronic-disease-management (CDM) feature schema suitable
for downstream clinical-decision tooling.

Architecture
------------

``mrnavax-rare`` is the **design side** (writes the mRNA therapeutic and
predicts how long it should persist in plasma). ``deepcatch-methylation``
is the **read side** (observes how much cfRNA is actually circulating
in the patient's plasma sample and flags anomalies against the
prediction).

The coupling is one-directional:

    mrnavax_rare.safety.cfrna.predict_cfrna_persistence
        → CfRNAPersistencePrediction
            → deepcatch_methylation.cfrna_coupling.couple_cfrna_to_cdm
                → RealTimeCDMObservation

Public surface::

    from deepcatch_methylation.cfrna_coupling import (
        CfRNAPersistenceInput,
        RealTimeCDMObservation,
        couple_cfrna_to_cdm,
        DosingRecommendation,
        ANOMALY_RATIO_HIGH,
        ANOMALY_RATIO_LOW,
    )

Honest disclosure
-----------------

The v0.10.0 release is **structural scaffolding only**. The real-time CDM
feature schema is defined; the dosing-recommendation thresholds are
honest-disclosure research-use defaults that will need clinical
calibration before any production use. No patient data has been
processed; no anomaly has been flagged in a real cohort. All
quantitative anchors cite the upstream
``mrnavax_rare.safety.cfrna.data`` module — no fabricated numbers
are introduced in this project.
"""

from __future__ import annotations

__version__ = "0.10.0"
__all__ = ["__version__", "get_version"]


def get_version() -> str:
    """Return the package version string (helper for tooling)."""
    return __version__
