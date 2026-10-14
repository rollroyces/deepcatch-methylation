# deepcatch-methylation

> Real-time cfRNA read-side companion to [mrnavax-rare](https://github.com/rollroyces/mrnavax-rare).
> Observes plasma cfRNA persistence and lifts it into a real-time
> chronic-disease-management (CDM) feature schema.

## v0.10.0 status

This is the **v0.10.0 initial scaffold** of the deepcatch-methylation
project. It ships:

- One coupling module: `deepcatch_methylation.cfrna_coupling`
- The `RealTimeCDMObservation` CDM feature schema
- A `CfRNAPersistenceInput` typed-dict contract that mirrors the
  output shape of `mrnavax_rare.safety.cfrna.predict_cfrna_persistence`
- Research-use anomaly + dosing-recommendation defaults

It does NOT yet ship:

- Wet-lab-validated anomaly thresholds (`ANOMALY_RATIO_LOW` /
  `ANOMALY_RATIO_HIGH` are research-use defaults)
- Real-patient data processing
- A live integration with a methylation-array or RNA-seq read pipeline

## Architecture

`mrnavax-rare` is the **design side**: it predicts how long a designed
mRNA therapeutic should persist in plasma.

`deepcatch-methylation` is the **read side**: it observes how much
cfRNA is actually circulating in the patient's plasma sample and
flags anomalies against the prediction.

```
mrnavax_rare.safety.cfrna.predict_cfrna_persistence
    → CfRNAPersistencePrediction
        → deepcatch_methylation.cfrna_coupling.couple_cfrna_to_cdm
            → RealTimeCDMObservation
```

The two projects share **no runtime code**. The coupling is a strict
TypedDict contract — the read side is deployable independently of the
design side.

## Quick start

```python
from deepcatch_methylation.cfrna_coupling import (
    CfRNAPersistenceInput,
    couple_cfrna_to_cdm,
)

# 1. Run the design-side prediction in mrnavax-rare:
#    prediction = mrnavax_rare.safety.cfrna.predict_cfrna_persistence(...)
#    prediction_dict = prediction.as_dict()

# 2. Lift into the CDM feature schema:
prediction: CfRNAPersistenceInput = {
    "predicted_half_life_hours": 54.0,
    "persistence_category": "medium",
    "target_tissue": "liver",
    "codon_optimization": "balanced",
}

obs = couple_cfrna_to_cdm(
    prediction=prediction,
    observed_cfRNA_concentration_pg_ml=200.0,
    sampling_time_post_dose_hours=24.0,
    dose_cfRNA_at_t0_pg_ml=1000.0,
)

print(obs.dosing_recommendation)  # "increase" / "maintain" / "decrease" / "escalate" / "hold"
print(obs.anomaly_flag)            # True / False
```

## Honest disclosure

- The v0.10.0 anomaly thresholds (`ANOMALY_RATIO_LOW = 0.5`,
  `ANOMALY_RATIO_HIGH = 2.0`) are **research-use defaults**. They
  will need clinical calibration before any production use.
- The dosing-recommendation rules are heuristic. They are NOT a
  substitute for clinician judgment.
- No patient data has been processed by this code path.
- No numerical anchor in this project is fabricated — every
  quantitative input comes from the upstream mrnavax-rare
  `cfrna_persistence_classifier` or is a research-use structural
  default.

## Tests

```bash
cd deepcatch-methylation
python3 -m pytest tests/
```

## Sister repos

- [mrnavax-rare](https://github.com/rollroyces/mrnavax-rare) — mRNA therapeutics design toolkit (v0.10.0 ships the upstream cfRNA persistence classifier).
- [deepcatch](https://github.com/rollroyces/deepcatch) — fragmentomics + multimodal cfDNA cancer detection.
