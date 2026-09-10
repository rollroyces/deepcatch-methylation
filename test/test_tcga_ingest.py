"""Tests for the real-data ingestion path (TCGA GDC → baseline)."""

import sys
from pathlib import Path

import pytest

# Add scripts/ to path so we can import run_phase0_baseline
_scripts = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(_scripts))

from run_phase0_baseline import load_tcga_subset  # noqa: E402


def _has_real_data():
    """Check if the real TCGA data has been downloaded."""
    data_dir = Path(__file__).resolve().parent.parent / "data" / "raw" / "tcga_lihc_subset"
    return (data_dir / "tumor_file_ids.txt").exists() and len(list(data_dir.glob("*.tsv"))) > 5


@pytest.mark.skipif(
    not _has_real_data(),
    reason="TCGA data not downloaded (run download_tcga_methylation.py first)",
)
def test_load_tcga_subset_real_data():
    """If TCGA data is present, verify the loader returns the expected shape and labels."""
    data_dir = Path(__file__).resolve().parent.parent / "data" / "raw" / "tcga_lihc_subset"
    beta_df, labels_df = load_tcga_subset(data_dir)
    assert beta_df.shape[0] >= 10  # at least 10 samples
    assert beta_df.shape[1] >= 100_000  # at least 100K CpG probes (450K array has ~485K)
    assert labels_df.shape == (beta_df.shape[0], 2)
    assert set(labels_df.columns) == {"sample_id", "label"}
    assert set(labels_df["label"].unique()).issubset({0, 1})
    assert labels_df["label"].sum() > 0  # at least one tumor
    assert (labels_df["label"] == 0).sum() > 0  # at least one normal


def test_load_tcga_subset_no_data_graceful():
    """If data dir doesn't exist, load should fail with a clear error."""
    fake_dir = Path("/tmp/nonexistent_tcga_test_dir")
    if fake_dir.exists():
        return  # skip if by some miracle it exists
    with pytest.raises(FileNotFoundError):
        load_tcga_subset(fake_dir)
