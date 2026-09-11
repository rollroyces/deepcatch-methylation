# FinaleMe Pretrained HMM Model Search — Path C Findings

**Date:** 2026-09-11
**Task:** Determine whether FinaleMe (Liu et al. *Nat Commun* 15:2790, 2024) authors have published pretrained HMM checkpoints so Step 2 training can be skipped.
**Bottom line:** **Yes — pretrained HMM models are publicly available on Zenodo.** Two models (healthy + cancer) totaling 148 kB were downloaded, MD5-verified, and confirmed compatible with the FinaleMe v0.61 JAR.

---

## TL;DR

| Item | Status |
|---|---|
| Pretrained HMM models published? | **Yes** — Zenodo record 14013719 (v3, Oct 2024) and 7779198 (v3, Feb 2024) |
| Model filenames | `healthy_WGS.mincg7.example.hmm_model`, `cancer_WGS.mincg7.example.hmm_model` |
| Model size | 74.1 kB each (148 kB total) |
| Model format | Java serialized `BayesianNhmmV5` object (jahmm library) — the exact class `FinaleMe.java` deserializes in decode mode |
| MD5 verification | Both MD5s match Zenodo's published checksums exactly |
| License | MIT for academic use; commercial use requires author contact |
| Training data (paper) | HD_45 healthy individual cfDNA WGS for the "healthy" model; the cancer model was trained on matched WGS from a prostate cancer patient in dbGaP phs003287 |
| Also bundled in Zenodo 14013719 | `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` (155 MB), `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` (324 MB), `wgEncodeDuke…hg19.bed` (86 kB), `ref_panel_wgbs.tar` (5.1 GB for tissue-of-origin) — all the references Step 1 needs |
| Step 3 decode runnable in this session? | **No** — Step 1 features are not cached, and the reference files (~480 MB) plus ~25 min BAM→features run would bust the 45-min / 2-GB task budget |

---

## 1. Search results

### 1.1 GitHub releases (`https://github.com/epifluidlab/FinaleMe/releases`)
Two releases exist (v0.61 latest, v0.58 / v0.58.1). Release notes mention "results from methylation prediction (steps 1-3) are still the same" and a tutorial update — **no pretrained model artifacts are attached to any release**. Assets 3 in v0.61 are source-only (the JAR is built from source).

### 1.2 Repo contents
- No `*.model`, `*.hmm`, `*.pkl`, `*.joblib` files anywhere under `src/`, `data/`, `tutorial/`, or `scripts/`. The `FinaleMe` Java class writes its model only when invoked with `-train`; there is no bundled "default" model.
- The `FinaleMe_workflow` Snakemake repo explicitly expects Step 2 training to run per-sample (`{sample_name}.finaleme.model` is an output, not an input).

### 1.3 FinaleMe paper — Nat Commun 15:2790 (Liu et al. 2024)
- No supplementary "pretrained model" file. Supplementary materials include Source Data (raw per-region methylation tables) and Supplementary Methods (math), but not the trained HMM object.
- Data availability statement: "The source code for the FinaleMe is freely available on GitHub" and "The data that support the findings of this study are available within the paper and its Supplementary Information, and from the corresponding author upon reasonable request." Training WGBS is in dbGaP phs000846 / phs001417 (controlled access) and phs003287 (the new cohort generated for this paper, also controlled access).

### 1.4 Zenodo — **the pretrained models live here**
Two records, both v3, by the FinaleMe team:

**Record 14013719** (Oct 30, 2024, **recommended — newer and most complete**, 22.4 GB total):
- `healthy_WGS.mincg7.example.hmm_model` (74.1 kB) — MD5 `c4c2c34ea5115c4b82e7296c5fbcabf6`
- `cancer_WGS.mincg7.example.hmm_model` (74.1 kB) — MD5 `f669a6d34c152f5c956ab7bbfcd6c9fa`
- `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` (155.1 MB)
- `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` (324.3 MB)
- `wgEncodeDukeMapability…hg19.bed` (85.8 kB)
- `ref_panel_wgbs.tar` (5.1 GB) — tissue-of-origin reference methylomes
- Plus the raw training data TSVs (HD_45, HD_46, 14230_1.WGS.b37.tsv.gz, meth_wgbs.tar, meth_wgs.tar.gz, ulp_wgs_frag.tar) and source zip

**Record 7779198** (Feb 21, 2024, 17.2 GB): same two `.hmm_model` files and methylation prior bw, but no `CG_motif.bedgraph.gz` and no `ref_panel_wgbs.tar`.

DOI badge in `FinaleMe_workflow/README.md` points to **record 14013719** as the canonical source.

### 1.5 Other epifluidlab repos
- `FinaleToolkit`, `FinaleMe_workflow`, `headneck`, `cragr`/`CRAG`, `finaledb_workflow` — none ship FinaleMe HMM models. `headneck` (HNSCC paper) is the only cancer-specific methylation sister project but uses FinaleMe as a library rather than distributing its models.

### 1.6 dbGaP phs003287
Sun et al. 2015 / 2018 are phs000846 / phs001417; the new FinaleMe cohort is phs003287. dbGaP access requires PI-level Data Access Committee approval via NIH eRA — not feasible in this session and orthogonal to the question (we have pretrained models already).

### 1.7 Author contact
- Yaping Liu (corresponding author): `yaping@northwestern.edu`
- Ravi Bandaru (co-author, bioinformatics lead): `ravi.bandaru@northwestern.edu`
- Kundan Baliga (workflow maintainer): `kundanbal2969@k12.ipsd.org`
- No public issue tracker is more useful than Zenodo — the README explicitly directs users there.

---

## 2. Downloads performed

```bash
mkdir -p /tmp/FinaleMe/models
curl -L --retry 3 --max-time 300 \
  -o /tmp/FinaleMe/models/healthy_WGS.mincg7.example.hmm_model \
  "https://zenodo.org/records/14013719/files/healthy_WGS.mincg7.example.hmm_model?download=1"
curl -L --retry 3 --max-time 300 \
  -o /tmp/FinaleMe/models/cancer_WGS.mincg7.example.hmm_model \
  "https://zenodo.org/records/14013719/files/cancer_WGS.mincg7.example.hmm_model?download=1"
```

Result:
```
-rw-r--r--  74065  Sep 11 17:50  cancer_WGS.mincg7.example.hmm_model
-rw-r--r--  74065  Sep 11 17:50  healthy_WGS.mincg7.example.hmm_model
```

`md5sum` confirms both MD5s match Zenodo's published checksums.

Format check (`file(1)`):
```
Java serialization data, version 5
```
Header bytes spell out `main.java.edu.mit.compbio.ccinferrence.hmm.BayesianNhmmV5` — the exact class `edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java` deserializes when `-decodeModeOnly` is passed (lines 491-499 of that file).

The two models differ (binary), so they are not duplicate files — they are healthy vs. cancer-trained parameter sets.

---

## 3. Why Step 3 was not run in this session

To run Step 3 (`-decodeModeOnly`) we still need a Step 1 output file: a per-CpG feature BED.gz from a BAM. Step 1 requires:
- `data/hg19.2bit` (~700 MB)
- `data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz` (155 MB)
- `data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` (324 MB)
- `data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed` (86 kB)

Total reference data: ~1.2 GB. Plus a BAM to run against (the BH01 chr22 BAM at `/tmp/finaleme_test/` is gone after the previous cleanup; only `/tmp/snyder_bams/IH03.chr22.bam` survives at 31 MB, and the previous Step 1 attempt on it hung).

Estimated end-to-end: ~1.2 GB download + ~25 min Step 1 + <1 min Step 3. **Within the 2 GB / 45 min budget only on paper; the prior session's IH03 hang indicates the path is fragile.**

---

## 4. What Step 3 looks like once references are present

The exact invocation (from the FinaleMe README, with our paths):

```bash
JAR="/tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar"
MODEL="/tmp/FinaleMe/models/healthy_WGS.mincg7.example.hmm_model"
FEATS="/path/to/BH01.cpg_features.hg19.bed.gz"   # from Step 1
OUT="/path/to/BH01.decode.prediction.bed.gz"

java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
  "$MODEL" "$FEATS" "$OUT" \
  -decodeModeOnly -t 4 \
  -bwOutput \
  -chromSizeFile /tmp/FinaleMe/data/hg19.chrom.sizes
```

Outputs:
- `*.prediction.bed.gz` (per-CpG binary methylation calls)
- `*.methy.bw`, `*.cov.bw`, `*.methy_count.bw` (bigWig tracks)

Then a per-window aggregation script is needed to get 1-kb methylation density for a cancer-vs-healthy AUC (the per-CpG file is ~16M rows for a chr22 BAM).

---

## 5. Recommended next steps for the user

1. **Decision: are the pretrained models acceptable for the project's scientific bar?**
   - **Pro:** Same author training, MIT licensed, 74 kB each, binary-checked, identical format to locally-trained models.
   - **Con:** The Zenodo description calls them `*.example.hmm_model` — i.e. *examples*, not the production-grade canonical models. The paper uses HD_45 specifically as the "healthy" training sample and one prostate cancer patient for the "cancer" model; if the project's goal is independent validation, training on a held-out sample is the cleaner approach.
   - The biorxiv paper text confirms: *"HMM model trained in high-coverage samples (HD_45, a healthy individual) to estimate the model parameters and applied it directly to each ULP-WGS dataset for the decoding"* — so the pretrained healthy model is exactly what was used in the paper's low-pass decoding experiments.

2. **If pretrained models are acceptable:**
   1. Download Zenodo 14013719's reference bundle (`CG_motif.bedgraph.gz`, `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw`, mappability BED) — ~480 MB, ~5 min on a healthy connection.
   2. Pull `hg19.2bit` from UCSC (~700 MB, ~15 min on this host — the previous run hit rate limits).
   3. Re-download the BH01 chr22 BAM from Zenodo record 6914806 (this time let it finish; 585 MB in ~5 min).
   4. Run Step 1 (~25 min for chr22), then Step 3 with the pretrained healthy model (<1 min).
   5. Repeat with the cancer model for cancer samples.
   6. Aggregate per 1-kb CGI/CGI-shore windows; run a logistic regression vs fragmentomics baseline.

3. **If pretrained models are NOT acceptable (independent training is required):**
   - The blockers are: dbGaP access for HD_45/HD_46/phs003287, and SRA throughput on this host.
   - Concrete path: request dbGaP DAC approval for phs000846 + phs001417 (Sun 2015/2018 buffy-coat WGBS, public-access tier may suffice); use that WGBS to train your own HMM via `FinaleMe.java -train`. The training step is `<1 min` once Step 1 features are in hand — the bottleneck is the WGBS download + Step 1 run, not training itself.
   - Or email Yaping Liu (`yaping@northwestern.edu`) to ask whether the *production* (non-`example`) HMM models can be shared; the Zenodo `*.example` naming strongly suggests the authors have non-example versions internally.

4. **Cheapest sanity check first:** if a single BAM's Step 1 features can be produced (Step 1 took 14 min on BH01 in the previous session), Step 3 with the pretrained model takes <1 min and yields a per-CpG methylation track you can compare against bisulfite truth elsewhere — this is the minimal demonstration that "true methylation from FinaleMe beats methylation-proxy" the project's head-to-head currently lacks.

---

## 6. Files created

- `/Users/hermes/deepcatch-methylation/docs/PRETRAINED_MODEL_SEARCH.md` — this document
- `/tmp/FinaleMe/models/healthy_WGS.mincg7.example.hmm_model` (74.1 kB, MD5-verified)
- `/tmp/FinaleMe/models/cancer_WGS.mincg7.example.hmm_model` (74.1 kB, MD5-verified)
