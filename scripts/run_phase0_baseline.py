#!/usr/bin/env python3
"""
Phase 0: load TCGA methylation β-values + run LR baseline.

Reads TCGA-LIHC methylation files (already-downloaded TSVs) into a feature
matrix, joins with sample-type labels (Primary Tumor = 1, Solid Tissue
Normal = 0), and runs a 5-fold pooled-OOF cross-study evaluation.

Usage:
    python scripts/run_phase0_baseline.py \
        --data-dir data/raw/tcga_lihc_subset \
        --output results/phase0_baseline.json
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# Add src/ to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from methylation.methylation_baseline import build_feature_matrix, evaluate_cv


def load_tcga_subset(data_dir: Path):
    """Load TCGA methylation TSVs into (beta_df, labels_df) pair.

    The filename convention:
      {file_id}.tsv

    The sample type is recorded separately in {tumor,n}_file_ids.txt
    written by the GDC download script.
    """
    # Parse which files are tumor vs normal from the sidecar files
    tumor_ids = set(open(data_dir / "tumor_file_ids.txt").read().split())
    normal_ids = set(open(data_dir / "normal_file_ids.txt").read().split())

    betas = {}
    labels = []
    for tsv in sorted(data_dir.glob("*.tsv")):
        fid = tsv.stem
        if fid in tumor_ids:
            label = 1
        elif fid in normal_ids:
            label = 0
        else:
            print(f"  WARNING: {fid} not in tumor/normal lists — skipping", file=sys.stderr)
            continue
        df = pd.read_csv(tsv, sep="\t", header=None, names=["probe_id", "beta"])
        df = df.set_index("probe_id")["beta"]
        betas[fid] = df
        labels.append({"sample_id": fid, "label": label})

    beta_df = pd.DataFrame(betas).T  # (n_samples, n_probes)
    labels_df = pd.DataFrame(labels)
    return beta_df, labels_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/raw/tcga_lihc_subset",
                    help="Directory containing TCGA *.tsv files + tumor/normal id lists")
    ap.add_argument("--output", default="results/phase0_baseline.json")
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--C", type=float, default=1.0)
    ap.add_argument("--min-var", type=float, default=0.01)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    if not (data_dir / "tumor_file_ids.txt").exists():
        sys.exit(f"Missing {data_dir / 'tumor_file_ids.txt'} — run the GDC download first")

    print(f"Loading TCGA-LIHC methylation data from {data_dir} ...")
    beta_df, labels_df = load_tcga_subset(data_dir)
    print(f"  Loaded {len(beta_df)} samples × {len(beta_df.columns)} CpG probes")
    n_tumor = int((labels_df['label'] == 1).sum())
    n_normal = int((labels_df['label'] == 0).sum())
    print(f"  Labels: {n_tumor} tumor, {n_normal} normal")

    print("Building feature matrix (filter low-variance probes) ...")
    X, y, feature_names = build_feature_matrix(beta_df, labels_df, min_var=args.min_var)
    print(f"  After filtering: {X.shape[0]} samples × {X.shape[1]} probes")

    print(f"Running {args.n_seeds}-seed × {args.n_splits}-fold pooled-OOF CV (C={args.C}) ...")
    result = evaluate_cv(X, y, n_seeds=args.n_seeds, n_splits=args.n_splits, C=args.C)
    print()
    print("=" * 60)
    print(f"Phase 0 baseline (TCGA-LIHC methylation only, LR L2 C={args.C}):")
    print(f"  AUC: {result['auc_mean']:.4f} ± {result['auc_std']:.4f}")
    print(f"  Per-seed AUCs: {[round(a, 4) for a in result['aucs']]}")
    print(f"  Samples: {result['n_samples']}, Features: {result['n_features']}")
    print("=" * 60)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({
            **result,
            "cohort": "TCGA-LIHC (12 tumor + 6 normal)",
            "data_source": "GDC API methylation β-values",
            "n_probes_total": int(beta_df.shape[1]),
            "n_probes_after_filter": int(X.shape[1]),
            "min_var_threshold": args.min_var,
        }, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
