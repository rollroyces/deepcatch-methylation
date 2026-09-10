#!/usr/bin/env python3
"""Download missing TCGA files for an existing project directory.

Used as a recovery tool when a previous download was interrupted.
Reads {tumor,normal}_file_ids.txt and downloads only the files not
already present (and complete) on disk.
"""
import argparse
import sys
from pathlib import Path

import requests

GDC_DATA_URL = "https://api.gdc.cancer.gov/data"


def download_file(file_id, out_path):
    # Consider a file complete if it has at least 20K probes (covers both 27K-filtered
    # and full 486K-probe GDC methylation files). Anything smaller is partial.
    if out_path.exists():
        with open(out_path, "rb") as f:
            nlines = sum(1 for _ in f)
        if nlines >= 20_000:
            return out_path.stat().st_size
        out_path.unlink()
    dl_url = f"{GDC_DATA_URL}/{file_id}"
    for attempt in range(3):
        try:
            r = requests.get(dl_url, timeout=180, stream=True)
            r.raise_for_status()
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(16384):
                    f.write(chunk)
            return out_path.stat().st_size
        except Exception as e:
            print(f"  retry {attempt+1} for {file_id[:12]}: {e}", file=sys.stderr)
            if out_path.exists():
                out_path.unlink()
    raise RuntimeError(f"failed to download {file_id}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--label", required=True, choices=["tumor", "normal"])
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    ids_path = data_dir / f"{args.label}_file_ids.txt"
    ids = [line.strip() for line in ids_path.read_text().splitlines() if line.strip()]
    missing = []
    for fid in ids:
        out = data_dir / f"{fid}.tsv"
        if out.exists():
            with open(out, "rb") as f:
                nlines = sum(1 for _ in f)
            if nlines >= 20_000:
                continue
            out.unlink()
        missing.append(fid)
    print(f"[{args.label}] {len(ids)} expected, {len(missing)} missing/incomplete")
    for fid in missing:
        out = data_dir / f"{fid}.tsv"
        size = download_file(fid, out)
        print(f"  {fid[:18]} {size//1024} KB")


if __name__ == "__main__":
    main()
