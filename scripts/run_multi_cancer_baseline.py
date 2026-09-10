#!/usr/bin/env python3
"""Multi-cancer TCGA methylation baseline — Phase 1 Option A.

Loads methylation β-values from data/raw/tcga_multi_cancer/{PROJECT}/ for
four cancer types (TCGA-LUAD, TCGA-BRCA, TCGA-COAD, TCGA-PAAD), plus the
existing TCGA-LIHC subset. Builds a unified feature matrix keyed on the
intersection of CpG probes present in all files, and runs three classifier
setups:

  a. Per-cancer-type baseline (5 cancer types × 5-seed CV)
  b. Pooled multi-cancer baseline (5-seed CV)
  c. Tumor-vs-tumor tissue-of-origin test (binary label per cancer type,
     trained on tumor samples only)

Usage:
    python scripts/run_multi_cancer_baseline.py \
        --multi-root data/raw/tcga_multi_cancer \
        --lihc-dir data/raw/tcga_lihc_subset \
        --output results/multi_cancer_baseline.json
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# Add src/ to path so we can import methylation_baseline unchanged
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from methylation.methylation_baseline import (  # noqa: E402
    build_feature_matrix,
    evaluate_cv,
)

CANCER_TYPES = ["TCGA-LUAD", "TCGA-BRCA", "TCGA-COAD", "TCGA-PAAD", "TCGA-LIHC"]


def load_one_dir(data_dir: Path, project: str):
    """Return (beta_df, labels_df) for one project directory.

    Mirrors scripts/run_phase0_baseline.py::load_tcga_subset so we don't
    duplicate logic.
    """
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
            continue
        df = pd.read_csv(tsv, sep="\t", header=None, names=["probe_id", "beta"])
        df = df.set_index("probe_id")["beta"]
        betas[fid] = df
        labels.append({"sample_id": fid, "label": label,
                       "project": project})
    beta_df = pd.DataFrame(betas).T
    labels_df = pd.DataFrame(labels)
    return beta_df, labels_df


def load_multi_cancer(multi_root: Path, lihc_dir: Path):
    """Load all five cancer types into a single (beta_df, labels_df) pair.

    Probes that don't appear in *every* sample are dropped via inner-join
    on the probe index, so all samples share the same feature space.
    """
    project_dirs = {
        "TCGA-LUAD": multi_root / "TCGA-LUAD",
        "TCGA-BRCA": multi_root / "TCGA-BRCA",
        "TCGA-COAD": multi_root / "TCGA-COAD",
        "TCGA-PAAD": multi_root / "TCGA-PAAD",
        "TCGA-LIHC": lihc_dir,
    }
    per_project = {}
    for project, d in project_dirs.items():
        if not (d / "tumor_file_ids.txt").exists():
            sys.exit(f"Missing {d}/tumor_file_ids.txt")
        beta_df, labels_df = load_one_dir(d, project)
        per_project[project] = (beta_df, labels_df)
        n_t = int((labels_df['label'] == 1).sum())
        n_n = int((labels_df['label'] == 0).sum())
        print(f"  {project}: {beta_df.shape[0]} samples × "
              f"{beta_df.shape[1]} probes ({n_t} tumor, {n_n} normal)")

    # Intersection of probe sets across all samples
    common_probes = None
    for beta_df, _ in per_project.values():
        s = set(beta_df.columns)
        common_probes = s if common_probes is None else (common_probes & s)
    assert common_probes is not None
    common_probes = sorted(common_probes)
    print(f"\nProbes common to every sample: {len(common_probes)}")

    # Stack
    beta_dfs, label_dfs = [], []
    for project, (beta_df, labels_df) in per_project.items():
        sub = beta_df[common_probes]
        beta_dfs.append(sub)
        label_dfs.append(labels_df)
    full_beta = pd.concat(beta_dfs, axis=0)
    full_labels = pd.concat(label_dfs, axis=0, ignore_index=True)
    return full_beta, full_labels, common_probes, per_project


def per_cancer_setup(per_project, common_probes, args):
    """5a. Per-cancer-type methylation baseline."""
    results = {}
    for project, (beta_df, labels_df) in per_project.items():
        sub = beta_df[common_probes]
        # build_feature_matrix expects only sample_id + label columns
        slim_labels = labels_df[["sample_id", "label"]].copy()
        X, y, feats = build_feature_matrix(sub, slim_labels, min_var=args.min_var)
        result = evaluate_cv(X, y, n_seeds=args.n_seeds,
                             n_splits=args.n_splits, C=args.C)
        n_t = int((y == 1).sum())
        n_n = int((y == 0).sum())
        results[project] = {
            **result,
            "n_tumor": n_t,
            "n_normal": n_n,
        }
        print(f"  {project}: AUC {result['auc_mean']:.4f} ± {result['auc_std']:.4f} "
              f"({n_t}t/{n_n}n, {result['n_features']} features)")
    return results


def pooled_setup(full_beta, full_labels, args):
    """5b. Pooled multi-cancer methylation baseline."""
    slim_labels = full_labels[["sample_id", "label"]].copy()
    X, y, feats = build_feature_matrix(full_beta, slim_labels,
                                        min_var=args.min_var)
    result = evaluate_cv(X, y, n_seeds=args.n_seeds,
                         n_splits=args.n_splits, C=args.C)
    print(f"  Pooled: AUC {result['auc_mean']:.4f} ± {result['auc_std']:.4f} "
          f"({result['n_samples']} samples, {result['n_features']} features)")
    return {
        **result,
        "n_tumor": int((y == 1).sum()),
        "n_normal": int((y == 0).sum()),
    }


def tissue_of_origin_setup(per_project, common_probes, args):
    """5c. Tumor-only tissue-of-origin classifier.

    Each pair of cancer types yields a separate binary task. We also run
    an OvR (one-vs-rest) pooled task where the label is "is cancer-type X".
    """
    results = {}
    # Concatenate tumor-only samples per project, tag with project_id
    tumor_dfs, tumor_labels = [], []
    proj_to_int = {p: i for i, p in enumerate(per_project.keys())}
    for project, (beta_df, labels_df) in per_project.items():
        tumor_rows = labels_df[labels_df["label"] == 1]
        sub = beta_df.loc[tumor_rows["sample_id"].values, common_probes]
        tumor_dfs.append(sub)
        lbl = tumor_rows.copy()
        lbl["project_id"] = proj_to_int[project]
        tumor_labels.append(lbl)
    tumor_beta = pd.concat(tumor_dfs, axis=0)
    tumor_lbl = pd.concat(tumor_labels, axis=0, ignore_index=True)

    # Pairwise binary classifiers
    pairwise = {}
    projects = list(per_project.keys())
    for i, pa in enumerate(projects):
        for pb in projects[i + 1:]:
            mask = tumor_lbl["project_id"].isin([proj_to_int[pa],
                                                  proj_to_int[pb]])
            sub_beta = tumor_beta.loc[mask.values]
            sub_lbl = tumor_lbl.loc[mask.values, ["sample_id", "label"]].copy()
            sub_lbl["label"] = (tumor_lbl.loc[mask.values, "project_id"]
                                == proj_to_int[pa]).astype(int).values
            # binary: 1 if pa, 0 if pb
            X, y, _ = build_feature_matrix(sub_beta, sub_lbl,
                                            min_var=args.min_var)
            res = evaluate_cv(X, y, n_seeds=args.n_seeds,
                              n_splits=args.n_splits, C=args.C)
            pairwise[f"{pa}_vs_{pb}"] = {
                **res,
                "n_class_pos": int((y == 1).sum()),
                "n_class_neg": int((y == 0).sum()),
            }
            print(f"  {pa} vs {pb}: AUC {res['auc_mean']:.4f} "
                  f"± {res['auc_std']:.4f} ({int((y==1).sum())}/{int((y==0).sum())})")
    results["pairwise"] = pairwise

    # Pooled OvR: pick the first project (TCGA-LUAD) as positive class vs
    # all others → a reasonable sanity-check single task.
    if "TCGA-LUAD" in proj_to_int:
        target = "TCGA-LUAD"
        sub_lbl = tumor_lbl[["sample_id"]].copy()
        sub_lbl["label"] = (tumor_lbl["project_id"]
                            == proj_to_int[target]).astype(int).values
        X, y, _ = build_feature_matrix(tumor_beta, sub_lbl,
                                        min_var=args.min_var)
        res = evaluate_cv(X, y, n_seeds=args.n_seeds,
                          n_splits=args.n_splits, C=args.C)
        results["ovr_TCGA-LUAD_vs_rest"] = {
            **res,
            "n_class_pos": int((y == 1).sum()),
            "n_class_neg": int((y == 0).sum()),
        }
        print(f"  OvR {target} vs rest: AUC {res['auc_mean']:.4f} "
              f"± {res['auc_std']:.4f} ({int((y==1).sum())}/{int((y==0).sum())})")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--multi-root", default="data/raw/tcga_multi_cancer")
    ap.add_argument("--lihc-dir", default="data/raw/tcga_lihc_subset")
    ap.add_argument("--output", default="results/multi_cancer_baseline.json")
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--C", type=float, default=1.0)
    ap.add_argument("--min-var", type=float, default=0.01)
    args = ap.parse_args()

    multi_root = Path(args.multi_root)
    lihc_dir = Path(args.lihc_dir)

    print("Loading multi-cancer TCGA methylation data ...")
    full_beta, full_labels, common_probes, per_project = load_multi_cancer(
        multi_root, lihc_dir)
    print(f"\nCombined: {full_beta.shape[0]} samples × {len(common_probes)} probes\n")

    print("Setup (a) per-cancer-type baseline ...")
    per_cancer = per_cancer_setup(per_project, common_probes, args)

    print("\nSetup (b) pooled multi-cancer baseline ...")
    pooled = pooled_setup(full_beta, full_labels, args)

    print("\nSetup (c) tumor-only tissue-of-origin test ...")
    tissue = tissue_of_origin_setup(per_project, common_probes, args)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "data_source": "GDC API methylation β-values",
        "cancer_types": CANCER_TYPES,
        "n_samples_per_project": {p: int(per_project[p][0].shape[0])
                                  for p in CANCER_TYPES},
        "n_tumor_per_project": {p: int((per_project[p][1]['label'] == 1).sum())
                                for p in CANCER_TYPES},
        "n_normal_per_project": {p: int((per_project[p][1]['label'] == 0).sum())
                                 for p in CANCER_TYPES},
        "n_probes_total_intersection": len(common_probes),
        "min_var_threshold": args.min_var,
        "C": args.C,
        "n_seeds": args.n_seeds,
        "n_splits": args.n_splits,
        "setup_a_per_cancer": per_cancer,
        "setup_b_pooled": pooled,
        "setup_c_tissue_of_origin": tissue,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
