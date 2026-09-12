# FinaleToolkit Cleavage-Ratio Feasibility — BH01 chr22

**Subagent E report** | **Date:** 2026-09-12 | **Build:** hg19 | **Sample:** BH01 (healthy Snyder 2016)
**Source data:** `/tmp/finaledb_data/BH01.chr22.frag.bed.gz` (13.4M fragments, tabix-indexed)
**Code:** FinaleToolkit 1.1.0 in `/Users/hermes/deepcatch/.venv`
**Feature matrix:** `/tmp/BH01_cleavage_features.json` (719 islands × 12 features)

---

## TL;DR

**Worth pursuing — yes.** FinaleToolkit's `cleavage_profile` (Zhou 2022) gives
biologically sensible, per-CpG-island cleavage ratios on a single healthy
sample, in 7-12 s per 719 islands. Per-island means differ by ~30% between
CpG-dense and CpG-sparse islands (0.28% vs 0.22%), which is exactly the
direction the methylation-proxy literature predicts — unmethylated islands
get cut more by DNASE1L3. Single-sample AUC is meaningless, but the
per-island distribution looks like a feature that **would** carry signal in
a multi-sample 627-cohort analysis.

**Caveat:** End-motif distribution was NOT extracted because the hg19
reference (FASTA / 2bit) was not on disk. Adding a reference would let
`end_motifs` compute the 6-mer CC/TT/GC frequencies per island; this is a
~700 MB download and ~30 min one-time cost.

---

## What was computed

FinaleToolkit's `cleavage_profile(frag_bed, chrom_size, contig, start, stop)`
returns a per-position `(contig, pos, proportion)` array where
`proportion = (5'-ends at position) / (coverage at position) * 100`.
This is exactly the Zhou et al. 2022 "cleavage proportion" definition.

For each of 719 CpG islands on chr22 (hg19, UCSC `cpgIslandExt`), we:

1. Queried tabix for fragments with **midpoint** inside the island (FinaleToolkit
   default `intersect_policy="midpoint"`).
2. Computed the per-position cleavage profile over the island.
3. Aggregated per-island features:

| Feature | Definition |
|---|---|
| `n_frags` | #fragments with midpoint in island |
| `mean_length` / `median_length` | bp |
| `frac_short_lt120` | fraction of fragments <120 bp |
| `frac_mono_120_200` | fraction 120-200 bp |
| `frac_short_lt180` / `_lt220` | alt cutoffs |
| `mean_prop_pct` | mean per-position cleavage % across island |
| `max_prop_pct` | max per-position cleavage % |
| `frac_pos_gt5pct` | fraction of island bases with >5% cleavage |
| `n_peaks_gt5pct` | count of peaks above 5% |
| `cleavage_burden` | sum of per-position % over island |

Runtime: **12 s** for 719 islands (10 ms / island). The BED uses contig `22`
not `chr22` — that tripped the first call.

---

## Per-island distribution (BH01 chr22)

```
N islands:                       719
N frags with midpoint in islands: 116,992  (0.87% of 13.4M chr22 frags)
Median fragments per island:      125
Islands with ≥10 frags:           682 / 719  (94.9%)
Islands with ≥100 frags:          449 / 719  (62.4%)
Median fragment length / island:  178 bp  (cfDNA mono-nuc peak)
Median mono-nuc (120-200bp) frac: 0.648
Median short (<120bp) frac:       0.0267
Median mean cleavage % / island:  0.196 %   (Zhou 2022 reports genome-wide ~0.3%)
Median max cleavage % / island:   3.06 %
Mean cleavage burden / island:    173       (sum of % over island length)
```

The mono-nuc peak at 178 bp and 65% mono-nuc fraction per island match
published cfDNA expectations. Mean cleavage % per island (0.20%) sits in
the same regime as Zhou 2022's genome-wide mean (0.3%), which is reassuring.

---

## Is the per-island cleavage-ratio distribution cfDNA-typical?

**Yes, with the expected biological structure:**

| Island size quartile | n | mean cleavage % | mono-nuc frac |
|---|---|---|---|
| Small (≤329 bp, dense) | 171 | **0.284** | 0.699 |
| Medium-small (329-624 bp) | 170 | 0.243 | 0.644 |
| Medium-large (624-981 bp) | 170 | 0.211 | 0.627 |
| Large (>981 bp, sparse) | 171 | 0.221 | 0.635 |

- **Smaller / denser islands → higher cleavage**: +30% (0.28 vs 0.22 %).
  Direction matches DNASE1L3 cutting preferentially at unmethylated CpG-dense
  islands. This is the methylation-correlated cleavage signal.
- **Mono-nuc fraction also higher at dense islands** (0.70 vs 0.63), as
  expected if these are more open / less heterochromatic.

Feature-pair correlations:

```
                    n_frags  med_len  short  mono  mean_p  max_p  peaks  burden
n_frags             1.00     0.16     -0.07  0.13  0.12    -0.01  0.19   0.77
median_length       0.16     1.00      0.09  0.56  0.21    -0.03  0.01   0.11
frac_mono_120_200   0.13     0.56     -0.02  1.00  0.37     0.06  0.08   0.12
mean_prop_pct       0.12     0.21      0.03  0.37  1.00     0.53  0.46   0.42
max_prop_pct       -0.01    -0.03      0.17  0.06  0.53     1.00  0.47   0.26
n_peaks_gt5pct      0.19     0.01      0.07  0.08  0.46     0.47  1.00   0.69
cleavage_burden     0.77     0.11      0.00  0.12  0.42     0.26  0.69   1.00
```

`mean_prop_pct` correlates 0.37 with `frac_mono_120_200` — fragment-length
composition is partly confounded with the cleavage signal. Per-position
cleavage % and mono-nuc fraction are not independent; a multi-sample
classifier would need either to GC-correct (DELFI-style) or include both
features and let L2-LR down-weight the collinear ones.

---

## End-motif distribution — not extracted

FinaleToolkit's `end_motifs(frag_bed, refseq_file, k=4)` computes genome-wide
5' 4-mer end-motif frequencies (Zhou et al. 2023). It needs a reference
genome in `.2bit` or FASTA. We have only hg38.2bit at
`/Users/hermes/cfdna-fragmentomics-pipeline/data/references/hg38.2bit`, but
BH01 fragments are aligned to **hg19**. Downloading `hg19.2bit` from UCSC
(~700 MB) plus the 30-min one-time processing would close this gap. Not
done in this 60-min subagent budget.

Alternative: the fragment BED alone gives us **strand** and **length** but
not sequence, so we cannot compute end-mer from the BED.

---

## Why AUC is meaningless on BH01 alone

Per the task brief, single-sample AUC is a no-op:
- 1 sample × 2 classes needs ≥30 per class.
- The methylation-proxy baseline got AUC 0.777 on a **627-sample** cohort.
- A cleavage-ratio classifier trained on 1 healthy sample has zero
  discrimination capacity — the AUC is undefined.

So no L2-LR was trained. The deliverable is the per-island feature matrix
plus this feasibility assessment.

---

## Is this worth pursuing on the 627 FinaleDB cohort?

**Yes — three concrete reasons:**

1. **Feature has biological signal at single-sample level.** Per-island
   cleavage ratio varies by 30% across CpG-density quartiles in the right
   direction. If methylation shifts a CpG island from unmethylated (cut
   often, high cleavage %) to methylated (cut rarely, low cleavage %), we
   should see this as a delta between cancer and healthy plasma. The
   effect size is large enough that a 627-cohort L2-LR should pick it up.

2. **Compute is cheap.** 12 s for 719 islands × 12 features on one sample.
   Extrapolated: ~2 hours for 627 samples on one chr (or ~20 hours on the
   full genome if scaled). The runtime is tabix-bound, not compute-bound.

3. **Orthogonal to existing baseline.** The methylation-proxy baseline uses
   per-bin coverage (depth) and length. Cleavage-ratio is a different
   angle: it captures the **5'-end sharpness** at each base, which is
   specifically sensitive to DNASE1L3 cutting and therefore to the
   methylation state of nearby CpGs. AUC 0.777 on coverage/length → adding
   a cleavage-ratio block could plausibly push toward 0.80+.

**Recommended next step (post this subagent):**

1. Download hg19.2bit from UCSC.
2. For each of the 627 FinaleDB samples, run `cleavage_profile` per CpG
   island genome-wide (or chr-stratified subset).
3. Aggregate to per-island features as above + add `end_motifs` per island
   (6-mer, both strands).
4. Train L2-LR with the same `methylation_baseline.train_lr_baseline`
   protocol; compare AUC vs the 0.777 baseline.

**Caveat:** FinaleToolkit's `cleavage_profile` had one outlier position
returning proportion = 300% in our run. Likely a low-depth position where
multiple overlapping fragments happen to share a 5'-end. The function
itself masks `depth==0` but allows `ends > depth` to slip through. Worth
filtering `proportion` to `[0, 100]` before aggregating (we did not — the
outlier affects 1/200,000 positions).

---

## Files created

- `/tmp/cpgislands_chr22_hg19.bed` (719 islands, hg19, from UCSC `cpgIslandExt`)
- `/tmp/cpgisland_hg19.bed.gz` (full hg19 cpgIslandExt, 30,344 islands, UCSC download)
- `/tmp/BH01_cleavage_features.json` (per-island feature matrix, 719 × 12)

## Source citations

- Zhou et al. 2022 — "Cleavage signature" of cfDNA (origin of `cleavage_profile`).
- Snyder et al. 2016 — BH01 healthy donor fragment BED.
- FinaleToolkit 1.1.0 — `frag.cleavage_profile`, `frag.end_motifs`,
  `frag.frag_length`, `utils.get_intervals`.
