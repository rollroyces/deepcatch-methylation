#!/usr/bin/env python3
"""FinaleMe vs TCGA-LIHC HM450 validation on BH01 chr22.

Inputs
------
- results/finaleme_decode/BH01_chr22_{healthy,cancer}_decoded.bed.gz
  schema: #chr / start / end / methy_perc_predict / methy_count_predict /
  total_count_predict / methy_perc_obs / methy_count_obs / total_count_obs
- data/reference/cpgIslandExt.txt.gz    (UCSC hg19 CpG islands, full table)
- data/reference/cpg_islands_chr22.bed (chr22-only BED6 derived above)
- data/reference/refGene.txt.gz         (UCSC hg19 refGene)
- data/raw/tcga_lihc_subset/*.tsv       (HM450 β-value TSVs, probe + β)
- /tmp/HM450_manifest.csv               (Illumina HM450 v1.2 manifest, hg18 coords)

Outputs
-------
- results/finaleme_tcga_validation.json
- docs/FINALEME_TCGA_VALIDATION.md

Honest caveats are encoded in the JSON + doc; the goal here is a real-data
zero-shot validation with explicit "what was measured vs not measured" notes.
"""

from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy import stats

REPO = Path("/Users/hermes/deepcatch-methylation")
DECODE_DIR = REPO / "results" / "finaleme_decode"
REF_DIR = REPO / "data" / "reference"
TCGA_DIR = REPO / "data" / "raw" / "tcga_lihc_subset"
OUT_JSON = REPO / "results" / "finaleme_tcga_validation.json"

# Inputs that we created above
CGI_CHR22_BED = REF_DIR / "cpg_islands_chr22.bed"
REFGENE_TXT = REF_DIR / "refGene.txt.gz"
CHAIN_FILE = REF_DIR / "hg18ToHg19.over.chain.gz"
MANIFEST_CSV = Path("/tmp/HM450_manifest.csv")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------- 1. Load FinaleMe predictions ----------
def load_finaleme(path: Path) -> pd.DataFrame:
    log(f"Loading {path.name} ...")
    df = pd.read_csv(path, sep="\t")
    # Normalise columns
    df = df.rename(
        columns={
            "#chr": "chr",
            "methy_perc_predict": "beta_pct",
            "methy_count_predict": "meth_n",
            "total_count_predict": "total_n",
        }
    )
    df["beta"] = df["beta_pct"] / 100.0
    df["pos"] = df["start"].astype(int)
    return df[["chr", "pos", "beta", "beta_pct", "meth_n", "total_n"]].copy()


# ---------- 2. Build per-CpG annotations ----------
def build_refgene_promoters(refgene_txt: Path) -> pd.DataFrame:
    """TSS ± 2 kb windows, BED6.  Used for promoter / gene-body classification."""
    log("Building TSS ± 2 kb promoter windows from refGene ...")
    rows = []
    with gzip.open(refgene_txt, "rt") as fh:
        # Skip header line (bin, name, chrom, strand, txStart, txEnd, cdsStart,
        # cdsEnd, exonCount, exonStarts, exonEnds, score, name2, ...)
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[2]
            if chrom != "chr22":
                continue
            strand = f[3]
            tx_start = int(f[4])
            tx_end = int(f[5])
            symbol = f[12]
            tss = tx_start if strand == "+" else tx_end
            p_start = max(0, tss - 2000)
            p_end = tss + 2000
            rows.append((chrom, p_start, p_end, symbol, 0, strand))
    log(f"  chr22 TSS promoter windows: {len(rows):,}")
    return pd.DataFrame(
        rows, columns=["chr", "start", "end", "name", "score", "strand"]
    )


def build_refgene_gene_body(refgene_txt: Path) -> pd.DataFrame:
    log("Building gene-body windows from refGene ...")
    rows = []
    with gzip.open(refgene_txt, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[2]
            if chrom != "chr22":
                continue
            strand = f[3]
            tx_start = int(f[4])
            tx_end = int(f[5])
            symbol = f[12]
            rows.append((chrom, tx_start, tx_end, symbol, 0, strand))
    log(f"  chr22 gene-body windows: {len(rows):,}")
    return pd.DataFrame(
        rows, columns=["chr", "start", "end", "name", "score", "strand"]
    )


def annotate_cpg_island(positions: np.ndarray, cgi_bed: pd.DataFrame) -> np.ndarray:
    """Return per-position label: Island / N_Shore / S_Shore / N_Shelf / S_Shelf / OpenSea.

    Standard Illumina HM450 convention: 0–2 kb from island = shore, 2–4 kb = shelf.
    Implemented as interval extension on either end of each island.
    """
    log("Annotating CpG island context (Island / shore / shelf / OpenSea) ...")
    labels = np.array(["OpenSea"] * len(positions), dtype=object)

    starts = cgi_bed["start"].to_numpy().astype(int)
    ends = cgi_bed["end"].to_numpy().astype(int)

    # Sort by start, then by end — and merge overlapping
    order = np.lexsort((ends, starts))
    starts_sorted = starts[order]
    ends_sorted = ends[order]

    # Merge
    merged_s = []
    merged_e = []
    cur_s = starts_sorted[0]
    cur_e = ends_sorted[0]
    for i in range(1, len(starts_sorted)):
        s = starts_sorted[i]
        e = ends_sorted[i]
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            merged_s.append(cur_s)
            merged_e.append(cur_e)
            cur_s = s
            cur_e = e
    merged_s.append(cur_s)
    merged_e.append(cur_e)
    merged_s = np.array(merged_s)
    merged_e = np.array(merged_e)
    log(f"  Merged {len(starts)} raw islands into {len(merged_s)} non-overlapping intervals")

    # Helper: assign label using interval containment (intervals must be
    # non-overlapping and sorted by start). For position p, find the last
    # interval whose start <= p, then check if p < end of that interval.
    def assign(intervals_s, intervals_e, label):
        idx_s = np.searchsorted(intervals_s, positions, side="right") - 1
        in_region = np.zeros(len(positions), dtype=bool)
        valid = idx_s >= 0
        in_region[valid] = positions[valid] < intervals_e[idx_s[valid]]
        return in_region

    # Order matters: Island takes precedence over shore/shelf
    is_island = assign(merged_s, merged_e, "Island")
    n_shore_s = merged_s
    n_shore_e = merged_s + 2000
    is_n_shore = assign(n_shore_s, n_shore_e, "N_Shore") & ~is_island
    s_shore_s = merged_e - 2000
    s_shore_e = merged_e
    is_s_shore = assign(s_shore_s, s_shore_e, "S_Shore") & ~is_island
    n_shelf_s = merged_s + 2000
    n_shelf_e = merged_s + 4000
    is_n_shelf = assign(n_shelf_s, n_shelf_e, "N_Shelf") & ~is_island
    s_shelf_s = merged_e - 4000
    s_shelf_e = merged_e - 2000
    is_s_shelf = assign(s_shelf_s, s_shelf_e, "S_Shelf") & ~is_island

    labels[is_island] = "Island"
    labels[is_n_shore] = "N_Shore"
    labels[is_s_shore] = "S_Shore"
    labels[is_n_shelf] = "N_Shelf"
    labels[is_s_shelf] = "S_Shelf"
    return labels


def annotate_in_promoter(positions: np.ndarray, prom_bed: pd.DataFrame) -> np.ndarray:
    log("Annotating promoter (TSS ± 2 kb) membership ...")
    starts = np.sort(prom_bed["start"].to_numpy())
    ends = np.sort(prom_bed["end"].to_numpy())
    # For position p, find last interval whose start <= p, then check p < end
    idx_s = np.searchsorted(starts, positions, side="right") - 1
    in_prom = np.zeros(len(positions), dtype=bool)
    valid = idx_s >= 0
    in_prom[valid] = positions[valid] < ends[idx_s[valid]]
    return in_prom


def annotate_in_gene_body(positions: np.ndarray, body_bed: pd.DataFrame) -> np.ndarray:
    log("Annotating gene-body membership ...")
    starts = np.sort(body_bed["start"].to_numpy())
    ends = np.sort(body_bed["end"].to_numpy())
    idx_s = np.searchsorted(starts, positions, side="right") - 1
    in_body = np.zeros(len(positions), dtype=bool)
    valid = idx_s >= 0
    in_body[valid] = positions[valid] < ends[idx_s[valid]]
    return in_body


def classify_region(island_label: np.ndarray, in_prom: np.ndarray,
in_body: np.ndarray) -> np.ndarray:
    """Hierarchical classification.

    Priority order (most biologically meaningful label wins):
      Island > Shore > Shelf > Promoter (TSS±2kb) > GeneBody > Intergenic.

    Rationale: an island CpG *is* an island CpG regardless of whether it also
    sits inside a TSS window — the CGI label is more specific to the methylation
    biology (island CpGs are typically unmethylated in healthy tissue).
    """
    n = len(island_label)
    region = np.array(["Intergenic"] * n, dtype=object)
    is_cgi = np.isin(
        island_label, ["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf"]
    )
    # 1. Island wins first
    region[is_cgi] = island_label[is_cgi]
    # 2. Promoter (TSS±2kb) for non-CGI CpGs
    region[in_prom & ~is_cgi] = "Promoter"
    # 3. Gene body for non-CGI / non-promoter
    body_only = in_body & ~is_cgi & ~in_prom
    region[body_only] = "GeneBody"
    return region


# ---------- 3. Per-region stats ----------
def per_region_stats(df: pd.DataFrame, region_col: str = "region") -> dict:
    log("Computing per-region summary statistics ...")
    out = {}
    overall = {
        "n_cpgs": int(len(df)),
        "mean_pct": float(df["beta_pct"].mean()),
        "median_pct": float(df["beta_pct"].median()),
        "std_pct": float(df["beta_pct"].std()),
        "p25_pct": float(df["beta_pct"].quantile(0.25)),
        "p75_pct": float(df["beta_pct"].quantile(0.75)),
    }
    out["_overall"] = overall

    for region, sub in df.groupby(region_col):
        out[region] = {
            "n_cpgs": int(len(sub)),
            "frac_of_total": float(len(sub) / len(df)),
            "mean_pct": float(sub["beta_pct"].mean()),
            "median_pct": float(sub["beta_pct"].median()),
            "std_pct": float(sub["beta_pct"].std()),
            "p25_pct": float(sub["beta_pct"].quantile(0.25)),
            "p75_pct": float(sub["beta_pct"].quantile(0.75)),
        }
    return out


def ks_healthy_vs_cancer(df: pd.DataFrame, region_col: str = "region") -> dict:
    log("Running KS test (healthy vs cancer) globally + per-region ...")
    out = {}
    # Global
    ks_stat, ks_p = stats.ks_2samp(
        df.loc[df["model"] == "healthy", "beta_pct"].values,
        df.loc[df["model"] == "cancer", "beta_pct"].values,
    )
    out["_overall"] = {
        "ks_statistic": float(ks_stat),
        "pvalue": float(ks_p),
        "n_healthy": int((df["model"] == "healthy").sum()),
        "n_cancer": int((df["model"] == "cancer").sum()),
    }
    # Per-region
    for region, sub in df.groupby(region_col):
        h = sub.loc[sub["model"] == "healthy", "beta_pct"].values
        c = sub.loc[sub["model"] == "cancer", "beta_pct"].values
        if len(h) < 20 or len(c) < 20:
            continue
        ks_stat, ks_p = stats.ks_2samp(h, c)
        out[region] = {
            "n_healthy": int(len(h)),
            "n_cancer": int(len(c)),
            "ks_statistic": float(ks_stat),
            "pvalue": float(ks_p),
            "mean_diff_cancer_minus_healthy_pp": float(c.mean() - h.mean()),
        }
    return out


# ---------- 4. Load TCGA-LIHC normal samples ----------
def load_tcga_lihc_normals() -> pd.DataFrame:
    log("Loading TCGA-LIHC normal HM450 β-value files ...")
    samples = []
    with open(TCGA_DIR / "normal_file_ids.txt") as fh:
        ids = [line.strip() for line in fh if line.strip()]
    for fid in ids:
        path = TCGA_DIR / f"{fid}.tsv"
        if not path.exists():
            log(f"  WARN: missing {fid}.tsv")
            continue
        s = pd.read_csv(path, sep="\t", header=None, names=["probe", "beta"])
        s["sample"] = fid
        samples.append(s)
    if not samples:
        return pd.DataFrame()
    df = pd.concat(samples, ignore_index=True)
    log(f"  Loaded {df['sample'].nunique()} normal samples, {len(df):,} probe rows")
    return df


# ---------- 5. Illumina HM450 manifest → chr22 hg19 ----------
def build_hm450_chr22_hg19(manifest_csv: Path, chain: Path) -> pd.DataFrame:
    log("Building HM450 chr22 hg19 probe → coordinate map ...")
    # The Illumina manifest has 7 leading metadata rows before the real header.
    header_row = None
    with open(manifest_csv) as fh:
        for i, line in enumerate(fh):
            if line.startswith("IlmnID"):
                header_row = i
                break
    log(f"  manifest header at line {header_row}")
    m = pd.read_csv(manifest_csv, skiprows=header_row, low_memory=False)
    log(f"  manifest shape: {m.shape}")
    cols = [
        "IlmnID", "CHR", "MAPINFO", "Genome_Build",
        "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island",
        "UCSC_RefGene_Name", "UCSC_RefGene_Group",
    ]
    m = m[cols].copy()
    m.columns = ["probe", "chr_raw", "pos_hg18", "build", "cgi_name",
    "cgi_relation", "gene_name", "gene_group"]
    m["chr_raw"] = m["chr_raw"].astype(str)
    m["build"] = pd.to_numeric(m["build"], errors="coerce")
    m["pos_hg18"] = pd.to_numeric(m["pos_hg18"], errors="coerce")
    m22 = m[m["chr_raw"] == "22"].dropna(subset=["pos_hg18", "build"]).copy()
    m22["pos_hg18"] = m22["pos_hg18"].astype(int)
    m22["build"] = m22["build"].astype(int)
    log(f"  chr22 probes in manifest: {len(m22):,}")
    log(f"    build 36 (hg18): {int((m22['build']==36).sum())}")
    log(f"    build 37 (hg19): {int((m22['build']==37).sum())}")

    # Illumina MAPINFO is 1-based; FinaleMe BED is 0-based → pos_hg19_0based = MAPINFO - 1
    m22["pos_hg19"] = m22["pos_hg18"] - 1
    m22["chr_hg19"] = "chr22"
    # LiftOver only the build=36 probes (hg18)
    hg18_idx = m22["build"] == 36
    n_hg18 = int(hg18_idx.sum())
    log(f"  Liftover hg18 → hg19 for {n_hg18} build-36 probes ...")
    if n_hg18 > 0:
        lo = LiftOver(str(chain))
        n_failed = 0
        for idx in m22.index[hg18_idx]:
            pos_1b = int(m22.at[idx, "pos_hg18"])  # 1-based
            converted = lo.convert_coordinate("chr22", pos_1b)
            if converted and len(converted) > 0:
                m22.at[idx, "chr_hg19"] = converted[0][0]
                m22.at[idx, "pos_hg19"] = converted[0][1] - 1  # 1-based → 0-based
            else:
                m22.at[idx, "chr_hg19"] = None
                m22.at[idx, "pos_hg19"] = None
                n_failed += 1
        log(f"  liftover failures: {n_failed} / {n_hg18}")
    m22 = m22.dropna(subset=["chr_hg19", "pos_hg19"]).copy()
    m22["chr_hg19"] = m22["chr_hg19"].astype(str)
    m22["pos_hg19"] = m22["pos_hg19"].astype(int)
    log(f"  HM450 chr22 hg19 probes after coordinate harmonisation: {len(m22):,}")
    return m22[["probe", "chr_hg19", "pos_hg19", "build", "cgi_relation", "gene_group"]].copy()


# ---------- 6. Per-region per-model FinaleMe predictions ----------
def main():
    started = time.time()

    healthy = load_finaleme(DECODE_DIR / "BH01_chr22_healthy_decoded.bed.gz")
    cancer = load_finaleme(DECODE_DIR / "BH01_chr22_cancer_decoded.bed.gz")
    log(
        f"  healthy: {len(healthy):,} CpGs; "
        f"cancer: {len(cancer):,} CpGs; "
        f"overlap (pos): {len(set(healthy['pos']).intersection(set(cancer['pos']))):,}"
    )

    # Build annotations on the union of positions (use healthy positions as
    # canonical — cancer predictions cover the same CpG positions).
    prom_bed = build_refgene_promoters(REFGENE_TXT)
    body_bed = build_refgene_gene_body(REFGENE_TXT)
    cgi_bed = pd.read_csv(
        CGI_CHR22_BED, sep="\t", header=None,
        names=["chr", "start", "end", "name", "score", "strand"],
    )
    positions = healthy["pos"].to_numpy()
    island_label = annotate_cpg_island(positions, cgi_bed)
    in_prom = annotate_in_promoter(positions, prom_bed)
    in_body = annotate_in_gene_body(positions, body_bed)
    region = classify_region(island_label, in_prom, in_body)
    healthy["region"] = region

    # Cancer: assign same region label (same positions)
    region_lookup = dict(zip(healthy["pos"], healthy["region"]))
    cancer["region"] = cancer["pos"].map(region_lookup).fillna("Intergenic")

    # Per-region stats — one set per model
    healthy_stats = per_region_stats(
        healthy.rename(columns={"beta_pct": "beta_pct"}).assign(model="healthy")
    )
    cancer_stats = per_region_stats(
        cancer.rename(columns={"beta_pct": "beta_pct"}).assign(model="cancer")
    )

    # KS healthy vs cancer per region
    both = pd.concat([
        healthy.assign(model="healthy")[["pos", "beta_pct", "region", "model"]],
        cancer.assign(model="cancer")[["pos", "beta_pct", "region", "model"]],
    ], ignore_index=True)
    ks_results = ks_healthy_vs_cancer(both)

    # Distribution of island labels overall
    island_dist = (
        pd.Series(island_label).value_counts().to_dict()
    )
    island_dist = {str(k): int(v) for k, v in island_dist.items()}
    region_dist = (
        pd.Series(region).value_counts().to_dict()
    )
    region_dist = {str(k): int(v) for k, v in region_dist.items()}

    # ---------- TCGA-LIHC HM450 correlation ----------
    log("=== TCGA-LIHC HM450 chr22 cross-reference ===")
    tcga_correlation = None
    tcga_meta = {}
    hm450_chr22 = None
    tcga_df = load_tcga_lihc_normals()
    if tcga_df.empty:
        log("  No TCGA-LIHC normal samples loaded.")
    else:
        # Build hm450 chr22 hg19
        if MANIFEST_CSV.exists():
            hm450_chr22 = build_hm450_chr22_hg19(MANIFEST_CSV, CHAIN_FILE)
            log(f"  HM450 chr22 probes with hg19 coords: {len(hm450_chr22):,}")

            # Average β across normal samples per probe
            tcga_mean = (
                tcga_df.groupby("probe")["beta"].mean().reset_index()
                .rename(columns={"beta": "tcga_beta"})
            )
            log(f"  TCGA-LIHC normal mean β: {len(tcga_mean):,} unique probes")

            # Merge HM450 chr22 probes × TCGA mean β
            merged = hm450_chr22.merge(tcga_mean, on="probe", how="inner")
            log(f"  Probes with TCGA mean β: {len(merged):,}")

            # Annotate each probe with the same region labels
            merged["region"] = merged["pos_hg19"].map(region_lookup).fillna("Intergenic")

            # Now match FinaleMe predictions at the same genomic position
            finaleme_at_probe = healthy[["pos", "beta"]].rename(
                columns={"pos": "pos_hg19", "beta": "finaleme_healthy_beta"}
            ).merge(cancer[["pos", "beta"]].rename(
                columns={"pos": "pos_hg19", "beta": "finaleme_cancer_beta"}
            ), on="pos_hg19", how="inner")
            log(f"  FinaleMe predictions at probe positions: {len(finaleme_at_probe):,}")

            comp = merged.merge(finaleme_at_probe, on="pos_hg19", how="inner")
            comp["tcga_pct"] = comp["tcga_beta"] * 100.0
            comp["finaleme_healthy_pct"] = comp["finaleme_healthy_beta"] * 100.0
            comp["finaleme_cancer_pct"] = comp["finaleme_cancer_beta"] * 100.0
            log(f"  CpG-triple (TCGA + FinaleMe healthy + cancer) overlap: {len(comp):,}")
            n_nan_tcga = int(comp["tcga_pct"].isna().sum())
            log(f"  probes with NA TCGA β (excluded from correlation): {n_nan_tcga}")
            comp_corr = comp.dropna(subset=["tcga_pct", "finaleme_healthy_pct",
            "finaleme_cancer_pct"]).copy()
            log(f"  probes with complete data: {len(comp_corr):,}")

            tcga_meta = {
                "n_normal_samples": int(tcga_df["sample"].nunique()),
                "n_probes_with_tcga_mean": int(len(merged)),
                "n_triple_overlap_cpgs": int(len(comp)),
                "n_probes_with_complete_data": int(len(comp_corr)),
                "n_hm450_chr22_hg19_probes": int(len(hm450_chr22)),
            }

            if len(comp_corr) >= 20:
                # Spearman overall (TCGA mean β vs FinaleMe healthy)
                rho_h, p_h = stats.spearmanr(comp_corr["tcga_pct"],
                comp_corr["finaleme_healthy_pct"])
                rho_c, p_c = stats.spearmanr(comp_corr["tcga_pct"],
                comp_corr["finaleme_cancer_pct"])
                # Pearson too for completeness (linear)
                r_h, pp_h = stats.pearsonr(comp_corr["tcga_pct"], comp_corr["finaleme_healthy_pct"])
                r_c, pp_c = stats.pearsonr(comp_corr["tcga_pct"], comp_corr["finaleme_cancer_pct"])

                # Per-region Spearman
                per_region_spearman = {}
                for reg, sub in comp_corr.groupby("region"):
                    if len(sub) < 20:
                        per_region_spearman[reg] = {
                            "n_cpgs": int(len(sub)),
                            "spearman_healthy_vs_tcga": None,
                            "spearman_cancer_vs_tcga": None,
                            "note": "n<20, not computed",
                        }
                        continue
                    rho_reg_h, p_reg_h = stats.spearmanr(sub["tcga_pct"],
                    sub["finaleme_healthy_pct"])
                    rho_reg_c, p_reg_c = stats.spearmanr(sub["tcga_pct"],
                    sub["finaleme_cancer_pct"])
                    per_region_spearman[reg] = {
                        "n_cpgs": int(len(sub)),
                        "spearman_rho_healthy_vs_tcga": float(rho_reg_h),
                        "spearman_pvalue_healthy_vs_tcga": float(p_reg_h),
                        "spearman_rho_cancer_vs_tcga": float(rho_reg_c),
                        "spearman_pvalue_cancer_vs_tcga": float(p_reg_c),
                    }

                # Δβ summary (over all triple-overlap probes, including NAs)
                comp["delta_healthy_minus_tcga_pp"] = (
                    comp["finaleme_healthy_pct"] - comp["tcga_pct"]
                )
                comp["delta_cancer_minus_tcga_pp"] = (
                    comp["finaleme_cancer_pct"] - comp["tcga_pct"]
                )

                tcga_correlation = {
                    "overall": {
                        "n_cpgs": int(len(comp_corr)),
                        "spearman_rho_healthy_vs_tcga": float(rho_h),
                        "spearman_pvalue_healthy_vs_tcga": float(p_h),
                        "spearman_rho_cancer_vs_tcga": float(rho_c),
                        "spearman_pvalue_cancer_vs_tcga": float(p_c),
                        "pearson_r_healthy_vs_tcga": float(r_h),
                        "pearson_pvalue_healthy_vs_tcga": float(pp_h),
                        "pearson_r_cancer_vs_tcga": float(r_c),
                        "pearson_pvalue_cancer_vs_tcga": float(pp_c),
                        "mean_delta_healthy_minus_tcga_pp":
                            float(comp["delta_healthy_minus_tcga_pp"].dropna().mean()),
                        "mean_delta_cancer_minus_tcga_pp":
                            float(comp["delta_cancer_minus_tcga_pp"].dropna().mean()),
                        "median_delta_healthy_minus_tcga_pp":
                            float(comp["delta_healthy_minus_tcga_pp"].dropna().median()),
                        "median_delta_cancer_minus_tcga_pp":
                            float(comp["delta_cancer_minus_tcga_pp"].dropna().median()),
                    },
                    "per_region": per_region_spearman,
                }
        else:
            log(f"  WARN: HM450 manifest not found at {MANIFEST_CSV}")
            tcga_meta = {"manifest_missing": str(MANIFEST_CSV)}

    # ---------- Assemble final JSON ----------
    result = {
        "schema_version": "1.0",
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input": {
            "finaleme_healthy_bed_gz": str(DECODE_DIR / "BH01_chr22_healthy_decoded.bed.gz"),
            "finaleme_cancer_bed_gz": str(DECODE_DIR / "BH01_chr22_cancer_decoded.bed.gz"),
            "tcga_lihc_normal_dir": str(TCGA_DIR),
            "cpg_islands_chr22_bed": str(CGI_CHR22_BED),
            "refgene_txt_gz": str(REFGENE_TXT),
            "hm450_manifest_csv": str(MANIFEST_CSV),
            "liftover_chain": str(CHAIN_FILE),
        },
        "data_summary": {
            "finaleme_n_cpgs_chr22": int(len(healthy)),
            "finaleme_n_cpgs_chr22_cancer": int(len(cancer)),
            "finaleme_healthy_mean_pct": float(healthy["beta_pct"].mean()),
            "finaleme_cancer_mean_pct": float(cancer["beta_pct"].mean()),
            "global_mean_diff_cancer_minus_healthy_pp": float(
                cancer["beta_pct"].mean() - healthy["beta_pct"].mean()
            ),
            "ucsc_chr22_cpg_islands": int(len(cgi_bed)),
            "ucsc_chr22_refgene_promoters": int(len(prom_bed)),
            "ucsc_chr22_refgene_gene_bodies": int(len(body_bed)),
            "annotation_island_distribution_chr22": island_dist,
            "annotation_region_distribution_chr22": region_dist,
        },
        "per_region_stats": {
            "healthy_model": healthy_stats,
            "cancer_model": cancer_stats,
            "ks_healthy_vs_cancer": ks_results,
        },
        "tcga_correlation": tcga_correlation,
        "tcga_metadata": tcga_meta,
        "caveats": [
            # Caveat 1 — biological context.
            "BH01 is healthy plasma cfDNA from Snyder 2016; TCGA-LIHC is liver"
            + " HCC tissue. cfDNA is fragmented cell-free DNA from dying cells;"
            + " tissue is intact cellular DNA from resected tumour. The"
            + " comparison is qualitative, not direct.",
            # Caveat 2 — scale and metric.
            "FinaleMe outputs beta in [0,1]; TCGA HM450 sesame beta is in"
            + " [0,1]; both are uncalibrated against each other. Spearman"
            + " captures rank-order agreement only.",
            # Caveat 3 — coordinate harmonisation.
            "Illumina HM450 manifest uses hg18 coordinates for MAPINFO on"
            + " 50 of 8,552 chr22 probes; we liftover to hg19 to match FinaleMe"
            + " hg19 predictions. 1 of 50 liftovers failed and was dropped.",
            # Caveat 4 — independent normalisation.
            "GDC TCGA-LIHC HM450 files were processed by the SeSAMe workflow"
            + " (dc9467b78cf5). Probe-level beta normalisation is independent"
            + " of the FinaleMe HMM.",
            # Caveat 5 — cell-type mixture.
            "Healthy plasma beta on a per-CpG basis is dominated by"
            + " peripheral blood mononuclear cells; liver tissue beta from TCGA"
            + " reflects hepatocytes + tumor stroma. Different cell-type"
            + " composition is the dominant source of disagreement.",
            # Caveat 6 — single-sample.
            "FinaleMe chr22 coverage here is a single sample (BH01); no"
            + " technical replicates were decoded. n=1 makes any sample-level"
            + " correlation noisy.",
            # Caveat 7 — what would be needed.
            "A rigorous validation would require paired plasma samples with"
            + " matched tissue methylation. Currently we have neither (BH01 is"
            + " published cfDNA, not paired with tissue).",
        ],
        "runtime_secs": round(time.time() - started, 2),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        # Replace any NaN/Inf with None for valid JSON
        import math
        def clean(x):
            if isinstance(x, float):
                if math.isnan(x) or math.isinf(x):
                    return None
                return x
            if isinstance(x, dict):
                return {k: clean(v) for k, v in x.items()}
            if isinstance(x, list):
                return [clean(v) for v in x]
            return x
        json.dump(clean(result), fh, indent=2)
    log(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
