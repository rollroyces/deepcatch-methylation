"""Parse and validate the FinaleMe output."""

import json

import pandas as pd


def main():
    print("Loading FinaleMe output...")
    df = pd.read_csv(
        '/tmp/FinaleMe/results/BH01.cpg_features.hg19.bed.gz', sep='\t'
    )
    print(f"Loaded {len(df):,} CpG-fragment data points")
    print(f"Columns: {list(df.columns)}")
    print()

    print("=== Per-CpG summary ===")
    per_cpg = df.groupby(['chr', 'start']).size()
    print(f"Unique CpGs: {len(per_cpg):,}")
    print(f"Fragments per CpG: mean={per_cpg.mean():.2f}, "
          f"median={per_cpg.median()}, max={per_cpg.max()}, min={per_cpg.min()}")

    print()
    print("=== Methylation state distribution ===")
    meth = df['methy_stat'].value_counts()
    print(meth)
    pct_meth = float(100 * meth.get('m', 0) / len(df))
    print(f"% methylated: {pct_meth:.2f}%")

    print()
    print("=== Fragment length distribution ===")
    print(f"Mean: {df['FragLen'].mean():.1f} bp")
    print(f"Median: {df['FragLen'].median():.1f} bp")
    print(f"Std: {df['FragLen'].std():.1f} bp")
    print(f"Min: {df['FragLen'].min()}, Max: {df['FragLen'].max()}")
    sub = int(((df['FragLen'] < 100)).sum())
    mono = int(((df['FragLen'] >= 100) & (df['FragLen'] <= 250)).sum())
    di = int(((df['FragLen'] >= 250) & (df['FragLen'] <= 400)).sum())
    print(f"Sub-nucleosomal (<100bp): {sub:,} ({100*sub/len(df):.1f}%)")
    print(f"Mono-nucleosomal (100-250bp): {mono:,} ({100*mono/len(df):.1f}%)")
    print(f"Di-nucleosomal (250-400bp): {di:,} ({100*di/len(df):.1f}%)")

    print()
    print("=== Coverage distribution (per CpG) ===")
    cov_per_cpg = df.groupby(['chr', 'start'])['Norm_Frag_cov'].first()
    print(f"Mean normalized coverage per CpG: {cov_per_cpg.mean():.4f}")
    print(f"Median: {cov_per_cpg.median():.4f}")
    print(f"Distribution percentiles: 10%={cov_per_cpg.quantile(0.1):.4f}, "
          f"90%={cov_per_cpg.quantile(0.9):.4f}, 99%={cov_per_cpg.quantile(0.99):.4f}")
    print(f"CpGs with very low coverage (Norm_Frag_cov < 0.001): {(cov_per_cpg < 0.001).sum():,}")
    print(f"CpGs with coverage >= 0.1: {(cov_per_cpg >= 0.1).sum():,}")

    summary = {
        "data_source": "FinaleMe v0.61 (Liu et al. Nat Commun 15:2790, 2024)",
        "input": "Snyder 2016 BH01 chr22 BAM (Snyder et al. Cell 2016, ~27M reads)",
        "output_format": "CpG-fragment table with methylation features",
        "n_data_points": int(len(df)),
        "n_unique_cpgs": int(len(per_cpg)),
        "mean_fragments_per_cpg": float(per_cpg.mean()),
        "pct_methylated": pct_meth,
        "fragment_length_mean_bp": float(df['FragLen'].mean()),
        "fragment_length_median_bp": float(df['FragLen'].median()),
        "fragment_length_std_bp": float(df['FragLen'].std()),
        "n_subnucleosomal_lt100bp": sub,
        "n_mononucleosomal_100_250bp": mono,
        "n_dinucleosomal_250_400bp": di,
        "pct_subnucleosomal_lt100bp": float(100 * sub / len(df)),
        "pct_mononucleosomal_100_250bp": float(100 * mono / len(df)),
        "pct_dinucleosomal_250_400bp": float(100 * di / len(df)),
        "mean_normalized_coverage_per_cpg": float(cov_per_cpg.mean()),
        "median_normalized_coverage_per_cpg": float(cov_per_cpg.median()),
        "n_cpgs_coverage_above_0.1": int((cov_per_cpg >= 0.1).sum()),
        "runtime_minutes": 14.31,
        "runtime_secs": 858.44,
        "jvm": "Java 21.0.12.1 aarch64, -Xmx20G, -XX:+UseParallelGC, -t 10",
        "platform": "Apple M4 (10 physical cores)",
        "chr_target": "chr22 only (filtered from full hg19 CpG index)",
    }
    out_path = '/Users/hermes/deepcatch-methylation/results/finaleme_smoke_validation.json'
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print()
    print(f"Summary saved to {out_path}")


if __name__ == "__main__":
    main()
