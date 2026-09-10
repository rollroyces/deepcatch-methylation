#!/usr/bin/env python3
"""
Download TCGA-LIHC (or other project) methylation β-values from GDC API.

GDC has no reCAPTCHA unlike NCBI GEO. Uses the standard GDC REST API.

Usage:
    python scripts/download_tcga_methylation.py \
        --project TCGA-LIHC \
        --n-tumor 12 \
        --n-normal 6 \
        --output-dir data/raw/tcga_lihc_subset
"""

import argparse
import json
import sys
from pathlib import Path

import requests

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"


def query_gdc_files(project, data_type="Methylation Beta Value",
                     sample_types=None, size=100):
    """Query GDC for files matching project + data_type + sample_types."""
    if sample_types is None:
        sample_types = ["Primary Tumor", "Solid Tissue Normal"]
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "=", "content": {"field": "data_type", "value": data_type}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": sample_types}},
        ]
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.case_id,cases.samples.sample_type,file_size",
        "size": str(size),
        "format": "json",
    }
    r = requests.get(GDC_FILES_URL, params=params, timeout=30)
    r.raise_for_status()
    return r.json()["data"]["hits"]


def download_file(file_id, out_path):
    """Download one GDC file to out_path."""
    if out_path.exists():
        return out_path.stat().st_size
    dl_url = f"{GDC_DATA_URL}/{file_id}"
    r = requests.get(dl_url, timeout=120, stream=True)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        for chunk in r.iter_content(8192):
            f.write(chunk)
    return out_path.stat().st_size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default="TCGA-LIHC",
                    help="TCGA project ID (e.g. TCGA-LIHC, TCGA-LUAD, TCGA-BRCA)")
    ap.add_argument("--n-tumor", type=int, default=12)
    ap.add_argument("--n-normal", type=int, default=6)
    ap.add_argument("--output-dir", default="data/raw/tcga_subset")
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Querying GDC for {args.project} methylation files ...")
    hits = query_gdc_files(args.project)
    tumor_hits = [h for h in hits if any(s.get("sample_type") == "Primary Tumor"
                                          for s in h.get("cases", [{}])[0].get("samples", []))]
    normal_hits = [h for h in hits if any(s.get("sample_type") == "Solid Tissue Normal"
                                            for s in h.get("cases", [{}])[0].get("samples", []))]
    print(f"  Available: {len(tumor_hits)} tumor + {len(normal_hits)} normal")
    print(f"  Requested: {args.n_tumor} tumor + {args.n_normal} normal")

    selected_tumor = tumor_hits[:args.n_tumor]
    selected_normal = normal_hits[:args.n_normal]
    if len(selected_tumor) < args.n_tumor:
        sys.exit(f"Only {len(selected_tumor)} tumor files available; need {args.n_tumor}")
    if len(selected_normal) < args.n_normal:
        sys.exit(f"Only {len(selected_normal)} normal files available; need {args.n_normal}")

    tumor_ids = [h["file_id"] for h in selected_tumor]
    normal_ids = [h["file_id"] for h in selected_normal]
    with open(out_dir / "tumor_file_ids.txt", "w") as f:
        f.write("\n".join(tumor_ids) + "\n")
    with open(out_dir / "normal_file_ids.txt", "w") as f:
        f.write("\n".join(normal_ids) + "\n")

    print(f"\nDownloading {len(tumor_ids) + len(normal_ids)} files to {out_dir} ...")
    for label, ids in [("tumor", tumor_ids), ("normal", normal_ids)]:
        for fid in ids:
            out_path = out_dir / f"{fid}.tsv"
            size_kb = download_file(fid, out_path) // 1024
            print(f"  [{label:6s}] {fid[:18]} {size_kb:6d} KB")
    print(f"\nDone. Files in {out_dir}: {len(list(out_dir.glob('*.tsv')))}")


if __name__ == "__main__":
    main()
