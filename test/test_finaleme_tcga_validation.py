"""Tests for the FinaleMe vs TCGA-LIHC HM450 validation script.

These tests verify the validation pipeline ran end-to-end and produced
honest, well-structured results — they do NOT make claims about
biological accuracy (a different sample and a different tissue type
than the training data; this is a zero-shot validation, not a
paired plasma↔tissue study).

Skip everything if the upstream FinaleMe decoded output or the TCGA
HM450 β-value TSVs are not present.
"""

import json
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "finaleme_validation" / "run_finaleme_tcga_validation.py"
OUT_JSON = REPO / "results" / "finaleme_tcga_validation.json"
DECODE_DIR = REPO / "results" / "finaleme_decode"
TCGA_DIR = REPO / "data" / "raw" / "tcga_lihc_subset"


def _has_real_data():
    """Both FinaleMe decoded output and TCGA HM450 files must exist."""
    if not DECODE_DIR.exists():
        return False
    if not (DECODE_DIR / "BH01_chr22_healthy_decoded.bed.gz").exists():
        return False
    if not (DECODE_DIR / "BH01_chr22_cancer_decoded.bed.gz").exists():
        return False
    if not (TCGA_DIR / "normal_file_ids.txt").exists():
        return False
    return len(list(TCGA_DIR.glob("*.tsv"))) >= 1


def _venv_python():
    """Project interpreter; matches the user's runbook convention."""
    return "/Users/hermes/deepcatch/.venv/bin/python"


@pytest.mark.skipif(
    not _has_real_data(),
    reason="FinaleMe decoded output or TCGA HM450 files not present",
)
def test_validation_script_runs_and_produces_json():
    """Re-run the validation script end-to-end and confirm the JSON exists."""
    if OUT_JSON.exists():
        OUT_JSON.unlink()
    result = subprocess.run(
        [str(_venv_python()), str(SCRIPT)],
        capture_output=True,
        text=True,
        env={"PATH": "/usr/bin:/bin", "PYTHONPATH": ""},
        timeout=300,
    )
    assert result.returncode == 0, (
        f"validation script failed:\nSTDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )
    assert OUT_JSON.exists(), "JSON output not produced"


@pytest.mark.skipif(
    not OUT_JSON.exists(),
    reason="Validation JSON not produced yet",
)
def test_validation_json_structure():
    """The JSON must have the required keys and reasonable ranges."""
    with open(OUT_JSON) as fh:
        d = json.load(fh)

    # Top-level schema
    assert d["schema_version"] == "1.0"
    for k in [
        "input",
        "data_summary",
        "per_region_stats",
        "tcga_correlation",
        "tcga_metadata",
        "caveats",
    ]:
        assert k in d, f"missing key: {k}"

    # Data summary checks
    ds = d["data_summary"]
    assert ds["finaleme_n_cpgs_chr22"] == 489_370
    assert ds["finaleme_n_cpgs_chr22_cancer"] == 489_370
    # Both means are β in [0, 100]
    assert 0.0 <= ds["finaleme_healthy_mean_pct"] <= 100.0
    assert 0.0 <= ds["finaleme_cancer_mean_pct"] <= 100.0
    # Global diff: cancer < healthy (cancer hypomethylation pattern)
    assert ds["global_mean_diff_cancer_minus_healthy_pp"] < 0.0

    # Per-region stats — Island β should be markedly lower than non-island
    p = d["per_region_stats"]
    if "Island" in p["healthy_model"] and "GeneBody" in p["healthy_model"]:
        island_mean = p["healthy_model"]["Island"]["mean_pct"]
        body_mean = p["healthy_model"]["GeneBody"]["mean_pct"]
        # Healthy islands are typically < 30 % methylated
        assert island_mean < body_mean, (
            f"Island mean β ({island_mean:.1f}) should be < GeneBody mean β ({body_mean:.1f})"
        )

    # KS test global — should be highly significant for n>100K
    ks = p["ks_healthy_vs_cancer"]["_overall"]
    assert ks["ks_statistic"] > 0.0
    assert ks["pvalue"] < 1e-10
    assert ks["n_healthy"] == ks["n_cancer"] == 489_370

    # TCGA correlation — if it ran, must be in [-1, 1]
    tc = d["tcga_correlation"]
    if tc is not None:
        ov = tc["overall"]
        if ov.get("spearman_rho_healthy_vs_tcga") is not None:
            assert -1.0 <= ov["spearman_rho_healthy_vs_tcga"] <= 1.0
            assert -1.0 <= ov["spearman_rho_cancer_vs_tcga"] <= 1.0
            assert ov["n_cpgs"] >= 100  # must have a reasonable overlap


@pytest.mark.skipif(
    not OUT_JSON.exists(),
    reason="Validation JSON not produced yet",
)
def test_validation_json_honest_caveats():
    """The JSON must include honest caveats — not silent extrapolation."""
    with open(OUT_JSON) as fh:
        d = json.load(fh)
    assert len(d["caveats"]) >= 5
    caveat_text = " ".join(d["caveats"]).lower()
    # Required honest statements
    assert "plasma" in caveat_text or "cfdna" in caveat_text
    assert "tissue" in caveat_text
    assert "qualitative" in caveat_text
