# Alternative Methylation Tools — Survey

**Date:** 2026-09-11
**Author:** Survey by subagent (Hermes Agent), commissioned by Yu Ching Lam
**Scope:** Survey only — no code, no implementation. Honest assessment of whether any
public tool can replace FinaleMe for our use case (cfDNA WGS → methylation features
on a 627-sample cross-study cohort on a single M4 Mac mini).

---

## Bottom line up front

**For methylation calling from plain cfDNA WGS (no bisulfite), the realistic
landscape in 2026 is extremely narrow.** Only three tools even attempt the task:
**FinaleMe** (the one we tried), **WGS2meth** (a brand-new 2025 CpG-island-level
classifier with no published cfDNA validation), and a third path — using
**FinaleToolkit's existing cleavage-ratio and end-motif features as a methylation
proxy** (already partially implemented in this repo, AUC 0.777 on 627 samples).

There is **no drop-in alternative to FinaleMe** that works on low-pass cfDNA WGS
with the maturity, validation, and TOO-deconvolution support we need. The honest
path forward is:

1. **Stop pursuing new methylation imputation tools** for this hardware/data regime.
2. **Document the methylation-proxy result honestly** (AUC 0.777, proxy already
   incorporated; combined +0.0006 AUC vs fragmentomics-only — published null).
3. **Pivot the project's headline methylation contribution to a different signal
   source** (see "Recommended path forward" below) rather than waste the 1–2 week
   time budget on tools that can't materially improve on the proxy result.

A focused investigation of WGS2meth and a deeper FinaleMe single-sample smoke run
remain scientifically interesting (≈2–3 days combined), but neither changes the
recommendation above. Detail and reasoning follow.

---

## 1. The constraint that filters the candidate list

Our inputs are **plain cfDNA WGS** at low-to-mid coverage (≈5–10× from FinaleDB;
some samples closer to 0.1× from ultra-low-pass studies). They have **no
bisulfite conversion**, **no enzymatic conversion (EM-seq)**, **no long reads
(ONT/PacBio)**, **no reduced representation (RRBS)**, **no methylation array**.

That single constraint eliminates the vast majority of methylation tools:

| Excluded because | Tools |
|---|---|
| **Require bisulfite sequencing** | Bismark, BWA-meth, MethylDackel, methylKit (with BS data), Bismap, cfMethyl-Seq, MethAtlas, cfNOMe toolkit, CelFiE, CelFEER, UXM, cfTools |
| **Require long-read / ONT / PacBio** | nanomethphase, DeepSignal, DeepMod, Nanopolish, Tombo, METEORE |
| **Require single-cell methylome coverage** | DeepCpG, MethylNet, Melissa, CaMelia, CpG-Transformer |
| **Require targeted methylation assay / enrichment** | cfMeDIP-seq tools, EM-seq callers, Roche SeqCap, TruSeq Methyl Capture |
| **Require methylation array (450K/850K/EPIC)** | minfi, ChAMP, EpiDISH, Houseman method, RefFreeEWAS |

What remains is a small set of tools that work on **plain WGS** — almost all of
which are **fragmentation-based imputation methods** (the FinaleMe family of
ideas). That family has exactly **two published, named tools**: FinaleMe and
WGS2meth. Everything else is either (a) a research paper without a released tool,
or (b) a re-implementation of the same biological signal.

---

## 2. Comprehensive survey table

| Tool / paper | Input | Output | Maturity / Last release | WGS-no-bisulfite compatible? | Viability for our task |
|---|---|---|---|---|---|
| **FinaleMe** (Liu 2024, *Nat Commun*) | cfDNA WGS BAM/frag.gz, hg19/hg38 | Per-CpG β, bigWig, TOO deconv | Active, MIT, v0.61, ~10 commits in 2024 | ✅ **Yes — designed for this.** Validated 16–39× deep + decode-only at 0.1×. | **Already tried — failed at multi-sample training (JVM hang, disk). Single-sample decode remains plausible** but requires solving the disk + Java issues first. 1–2 weeks to a cohort feature matrix is optimistic. |
| **WGS2meth** (Abdullaev 2025, *BMC Genomics*) | Standard WGS BAM (any depth) | Per-CpG-island binary methylation status | **Just published.** Paper has a pre-trained XGBoost model, code in supplements but **no released GitHub repo found**. | ✅ **Yes — designed for this.** Trained on cell-line WGS; AUC 0.83–0.96 reported. | **Promising in principle** (newest tool, designed for the same problem). **In practice blocked**: no public repo to clone, no docs, no CI, no license, no per-fragment WGS-vs-WGBS validation. 3–7 days to reproduce + validate. High risk of orphan-code trap. |
| **FinaleToolkit (already installed)** | frag.gz / BAM from FinaleDB-style workflow | 10+ features including **cleavage proportion** (CpG cleavage ratio — Zhou 2022), end-motif, MDS, WPS, DELFI | Active, MIT, v0.10.7+, well-maintained | ✅ **Yes — designed for this.** Authored by the same lab as FinaleMe. | **Already implemented.** The "methylation-proxy" features we ran are exactly what FinaleToolkit extracts at the aggregated level. **Cleavage proportion over CpG-rich windows** is the closest thing to a methylation feature FinaleToolkit exposes. Re-extracting cleavage-ratio specifically at annotated CpG islands on raw fragments (not the existing aggregated .npy) could meaningfully improve the proxy — see Section 4. |
| **methylKit** (R/Bioconductor, Akalin 2012) | BS-seq, targeted-BS, or array-derived methylation calls | Per-CpG methylation stats / DMRs | Mature, MIT, actively maintained | ❌ **Requires pre-computed methylation calls** (or BS-seq input). Cannot ingest plain WGS. | **Not viable.** Useless for WGS without an upstream caller, and the upstream caller problem is exactly what we're trying to solve. |
| **Bismark / BWA-meth / MethylDackel** | Bisulfite-converted reads | Per-CpG methylation from C→T conversion | Mature, GPL | ❌ **Requires bisulfite sequencing.** No bisulfite → no C→T signature → no output. | **Not viable.** |
| **Bismap** | Bisulfite reads | Aligned reads for downstream methylation calling | Older | ❌ **Requires bisulfite sequencing.** | **Not viable.** |
| **nanomethphase** (ONT) | ONT BAM with modified-base tags | Phased methylation haplotypes | Active, MIT, 2024 | ❌ **Requires ONT reads.** Our WGS is short-read Illumina. | **Not viable.** |
| **DeepSignal / DeepMod / Nanopolish** | ONT/PacBio signal or BAM | Per-read methylation | Various | ❌ **Requires long-read data.** | **Not viable.** |
| **cfNOMe toolkit** | BS-seq (NOMe-style accessibility + methylation) | Cell-type proportions | 2020, MIT | ❌ **Requires nucleosome occupancy methylation (NOMe) assay.** | **Not viable** — wrong assay. |
| **cfMeDIP-seq tools** | cfMeDIP enriched reads | Methylation-enriched regions | Various | ❌ **Requires MeDIP enrichment.** | **Not viable** — wrong assay. |
| **cfTools** (Hu 2025, *Bioinformatics Advances*) | Bismark-processed WGBS / cfMethyl-Seq | Cell-type proportions, tumor fraction | Active, R/Bioconductor, MIT | ❌ **Requires pre-computed methylation calls.** CancerDetector + cfDeconvolve + cfSort. | **Not viable** — consumes methylation calls, doesn't produce them. Could become useful **after** we have methylation features from another source. |
| **CelFiE / CelFiE-ISH / MethAtlas / UXM / cfNOMe toolkit** | WGBS-derived methylation calls at TIMs (tissue-informative markers) | Cell-type deconvolution | All active, GPL/MIT | ❌ **Require methylation calls** from WGBS/cfMethyl-seq. | **Not viable** as primary caller. Same downstream-only status as cfTools. |
| **DISMIR** (Li 2021, *Brief Bioinform*) | WGBS | Per-read cancer/tissue origin classification | 2021, last commit 2021 | ❌ **Requires WGBS.** | **Not viable.** |
| **DeepCpG** (Angermueller 2017, *Genome Biol*) | Single-cell BS-seq (scBS-seq) | Imputed single-cell methylation | 2017, last commit 2018 | ❌ **Single-cell BS-seq.** | **Not viable** — wrong scale, wrong input. |
| **MethylNet** (Levy 2020, *BMC Bioinform*) | Methylation array or BS-seq β-value matrix | Deep models for downstream tasks (deconv, classification, age) | 2020, last commit 2020, unmaintained | ❌ **Consumes a methylation matrix.** | **Not viable** as primary caller. Same downstream-only status. |
| **CelFiE-ISH ONT extension** (Unterman 2026, bioRxiv) | ULP-WGS ONT | Cell-type deconvolution including unknown types | 2026 preprint | ❌ **Requires ONT.** | **Not viable.** |
| **cfMethyl-Seq** (Stackpole 2022, *Nat Commun*) | cfMethyl-Seq library prep (cost-optimized BS) | WGBS-equivalent methylation calls | Wet-lab protocol + bioinformatics | ❌ **Requires new library prep.** | **Not viable** for our existing WGS data. |
| **In-house XGBoost on fragment coordinates** (the WGS2meth approach) | Any WGS BAM | Per-CpG-island binary methylation | Not a tool — paper-only | ✅ In principle. | Could be reproduced from scratch in 1–2 weeks (the model is 11-bp windows around CpGs, XGBoost). But it's exactly what FinaleMe does at finer granularity. |

**Summary count:** Out of ~20 surveyed tools, **2 work on plain cfDNA WGS**
(FinaleMe + WGS2meth), **1 already provides methylation-correlated features from
the same data pipeline (FinaleToolkit)**, and the rest are either downstream-only
or require bisulfite / long-read / single-cell / targeted-assay data that we don't have.

---

## 3. The two WGS-compatible tools — deep-dive

### 3.1 FinaleMe (already attempted — known failure mode)

| Property | Value |
|---|---|
| Repository | https://github.com/epifluidlab/FinaleMe |
| License | MIT (academic); commercial license required |
| Input | BAM (or FinaleDB-style frag.gz, hg19/hg38); reference fasta/2bit + CpG annotation + WGBS-prior bigWig |
| Output | Per-CpG β-value BED; methylation bigWig; coverage bigWig; TOO deconvolution (per-tissue fraction) |
| Hardware | Java 21 + Maven. Per-sample Step 1 takes ~25 min for chr22 BAM (≈100×); the full-genome equivalent is multi-hour, multi-tens-of-GB-RAM. Recommended `-Xmx20G` per JVM. |
| Last activity | v0.61; active 2024 |
| Coverage regime validated | Deep 16–39× for training; **decode-only at 0.1× using a model trained on HD_45 healthy sample** |
| Our prior experience | Per-sample feature extraction on Snyder 2016 chr22 BAM hung at <1% CPU. Multi-sample training failed because disk + JVM combination kept OOMing on the M4 / 16 GB / 50 GB setup. |

**Honest re-assessment.** The JVM-hang failure we observed is probably a
**Step-1 (feature matrix) bug, not a fundamental limit**. Step 1 needs to
read each BAM once, emit CpG records; it should not hang. But diagnosing the
hang without re-attempting on real BAMs costs days we don't have.

**Coverage concern.** FinaleMe's headline auROC of **0.91 was measured at 16–39×
deep WGS**. At decode-only 0.1×, performance drops substantially — the paper
explicitly notes they trained on HD_45 (a healthy high-coverage sample) and
applied that model directly to ULP-WGS without retraining. Our cohort's 5–10×
is in the middle of these regimes: **the model probably works better than at
0.1× but worse than at 16×**. Whether the predicted β-values carry meaningful
cancer-vs-healthy signal at our coverage is **unknown** and would need a paired
WGBS ground truth (which we don't have access to) to validate. This is the
"validated auROC 0.91" being taken on faith in the project plan — that
number is for fragments with ≥5 CpGs in CpG-rich regions at deep coverage,
not for our actual data.

**Time to a usable cohort matrix (honest estimate).** Assuming the JVM hang is
fixable:
- **3–5 days** to get Step 1 working on a single Snyder 2016 BAM end-to-end
  (debug, JDK 21 install path, reference file resolution).
- **1–2 days per sample** for Step 1 once debugged (parallelism limited by RAM).
- **~627 samples × 1 day = 2 years single-threaded**, or **2–4 weeks** with
  parallelization constrained by 50 GB disk and 16 GB RAM.
- Plus Step 2–4 (training, decoding, deconvolution), another 2–3 days total.

This **violates the 1–2-week budget** even in the optimistic case. And the
payoff is uncertain because the validated auROC may not transfer to 5–10×
coverage.

### 3.2 WGS2meth (Abdullaev 2025, *BMC Genomics* 26:973)

| Property | Value |
|---|---|
| Paper | https://link.springer.com/article/10.1186/s12864-025-12257-7 |
| Repository | **Not found.** Author GitHub `Abdullaev6` only contains unrelated repos. Authors affiliated with Max Planck Institute for Molecular Genetics (Berlin) and BIH Charité. |
| License | Not stated (CC-BY 4.0 paper) |
| Input | WGS BAM aligned to hg38 (no bisulfite) |
| Output | Binary CpG-island methylated / unmethylated labels (~30k CpG islands) |
| Method | XGBoost classifier on dinucleotide fragmentation rates (ρ^XY) per CpG island, computed from 5'-end read mapping. Pre-trained on cell-line WGS + matched WGBS methylation. |
| Hardware | Pure Python + XGBoost. Genome-wide per-sample compute likely under an hour on a modern CPU. |
| Last activity | Published Sep 2025 — brand new. |
| Validation in paper | AUC 0.83–0.96 on cell-line WGS with matched WGBS. **No validation on cfDNA reported** (only cell-line / normal-tissue WGS). |
| Coverage regime | Trained at "standard WGS coverage"; not explicitly tested at low-pass cfDNA depths. |

**What's interesting.** WGS2meth's biological premise is the **same** as
FinaleMe's — fragmentation biases at CpG sites encode methylation — but the
implementation is **completely different**: an XGBoost classifier on
fragmentation rates per CpG island vs. a non-homogeneous HMM on per-fragment
feature vectors. The signal is real and replicated (FinaleMe auROC 0.91 at
deep coverage; WGS2meth AUC 0.83–0.96 on cell lines). WGS2meth is also
**dramatically simpler** to run — no Java, no Maven, no per-CpG output, just a
pre-trained XGBoost model applied per CpG island.

**Why it's blocked for us right now.**
1. **No public code repository.** The paper says the model is "distributed with
   the code" but I could not find it on GitHub, Bitbucket, Zenodo, or the
   authors' personal pages as of 2026-09-11. (Possible locations: paper
   supplementary, Max Planck GitLab — would need institutional access or
   emailing the corresponding author.)
2. **No validation on cfDNA.** Cell-line WGS is the test material. cfDNA has
   a completely different fragmentation biology (endogenous nucleases vs.
   sonication). The model may not transfer at all — this is exactly the issue
   FinaleMe addressed by training on healthy cfDNA samples.
3. **Output is binary CpG-island methylation** (~30k features per sample),
   not continuous β-values. Limited utility for downstream TOO deconvolution
   and incompatible with methylation-probability-style features.

**Time to a usable cohort matrix (honest estimate).**
- **1–3 days** to obtain the model (email author, wait, or recreate from paper
  description — the XGBoost features are well-specified in the methods).
- **2–5 days** to implement a minimal WGS BAM → CpG-island fragmentation-rate
  pipeline from scratch if no code is available.
- **< 1 hour per sample** to apply the model once implemented.
- **1–2 days** to train an LR baseline and benchmark against the existing
  methylation-proxy result (AUC 0.777) and fragmentomics (AUC 0.9746).

**Verdict.** **Scientifically the most interesting lead in the survey** but
practically a coin-flip on whether the missing code can be obtained and whether
cfDNA transfer will work. Worth **2–3 days of investigation** (email authors +
try to obtain the model), with a hard pivot back if the code doesn't materialize.

---

## 4. The third path: deeper use of FinaleToolkit (already in our hands)

This survey surfaces a finding the project has not yet exploited: **FinaleToolkit
explicitly implements a "Cleavage Proportion" feature designed by Zhou et al.
2022 that the authors document as correlated with DNA methylation at CpG sites**.

From the FinaleToolkit paper:
> "A recent study suggested a tight correlation between the cleavage ratios
> near CpG sites and the DNA methylation status at CpGs. Our FinaleToolkit
> enabled the efficient extraction of the genome-wide cleavage ratios feature
> from deep cfDNA WGS and observed the anti-correlated phasing patterns
> between nucleosomes (WPS) and cleavage ratio near CTCF."

What we already have from the methylation-proxy head-to-head:
- Fragmentomics-only: AUC 0.9746
- **Methylation-proxy-only: AUC 0.777** (29 summary features)
- Combined: AUC 0.9752 (Δ +0.0006, statistically significant but practically
  negligible)

What we **do not yet have**, but could extract with FinaleToolkit:
- **Per-CpG-island cleavage proportion** (Zhou 2022) — directly methylation-correlated,
  not yet in the proxy feature set.
- **Per-CpG-island end-motif distribution** (Jiang 2020, also noted by FinaleToolkit
  authors as methylation-related).
- **Cleavage ratios near CTCF / TSS / TFBS** (anti-correlated with methylation,
  same paper).

These are **distinct biological signals** from the existing DELFI short/long
ratios and FSD. The published evidence (Zhou 2022, Jiang 2020, FinaleToolkit
paper) suggests they carry orthogonal methylation information. **Time to
implement + run + LR baseline on the 627-sample cohort: 3–5 days**, on
hardware we already have, with software already installed.

**This is the highest-expected-value remaining methylation path.** It doesn't
require new tools, new raw data, new funding, or new collaborators. It does
require reaching back to raw FinaleDB `.frag.tsv.bgz` files (which we haven't
done yet — the existing pipeline only consumes the pre-aggregated 5-channel
features).

The honest null risk: these features may also fail to lift AUC beyond the
existing +0.0006, just like the 29-feature proxy did. **But** the test is
cheap (days, not weeks) and the signal has independent published support.

---

## 5. Recommended path forward

In priority order, with honest time estimates and exit criteria:

### Tier 1 — Run within the existing budget (3–5 days)

**Action:** Extract per-CpG-island cleavage proportion, per-CpG-island end-motif
distribution, and cleavage-ratio-near-CTCF features from the existing
FinaleDB `.frag.tsv.bgz` files using FinaleToolkit (already installed).
Train an LR baseline on the 627-sample cohort and compare against the
methylation-proxy AUC 0.777 and fragmentomics AUC 0.9746.

**Why first:** It's the cheapest test with the highest prior probability of
adding real signal (documented independent biology, orthogonal to existing
features, no new tool acquisition risk). It also doesn't compete with the
ongoing fragmentomics baseline work — same data, same labels, same CV protocol.

**Exit criteria:**
- If new cleavage-ratio / per-island end-motif features lift combined AUC by
  ≥ +0.005 over the current 0.9752 → publishable methylation contribution.
- If they don't (Δ < +0.002, like the existing proxy) → publish the null
  result and move to Tier 2.

### Tier 2 — Quick WGS2meth investigation (2–3 days)

**Action:** Email the corresponding author (Eldar T. Abdullaev, Max Planck
Institute) to request the model + reference data. Also check the BMC Genomics
supplementary materials for code archives. If obtained, run on 10 samples
(5 cancer + 5 healthy) and benchmark.

**Exit criteria:**
- Code obtained + 10-sample AUC > 0.6 on cancer-vs-healthy → Tier 3 (cohort
  extraction).
- Code not obtainable, or 10-sample AUC ≤ chance → drop WGS2meth permanently
  and publish the survey finding as "no public tool beyond FinaleMe can
  deliver methylation features on plain cfDNA WGS in 2026".

### Tier 3 — Conditional FinaleMe retry (only if Tiers 1 and 2 both null, 1–2 weeks)

**Action:** Retry FinaleMe's Step 1 on a single Snyder 2016 BAM (already
downloaded chr22 BAM exists). Diagnose the JVM hang. If a single sample runs
end-to-end in < 24 hours and produces above-chance signal on a known cancer
sample, scale to a 20-sample subset to measure time-to-cohort realistically.

**Exit criterion:** Even a successful 20-sample FinaleMe run with
documented signal would only validate that the tool works on our data — the
full 627-sample feature matrix is still 2–4 weeks of compute. The cost-benefit
of Tier 3 is poor and it's listed last on purpose.

### Stop condition

If Tier 1 (FinaleToolkit deep dive) and Tier 2 (WGS2meth) both return null
results, **declare the methylation imputation path exhausted** and reframe
the project's methylation contribution as:
- "Methylation-proxy (cleavage-ratio / end-motif derived from FinaleToolkit
  aggregated features) reaches AUC 0.777 on a 627-sample cross-study cohort."
- "Fragmentomics saturates at AUC 0.975; methylation-proxy adds +0.0006,
  below the practical-significance threshold."
- "Documenting the failure modes of FinaleMe + WGS2meth + FinaleToolkit deep
  proxy is itself a useful methodological contribution for the solo-developer
  cfDNA methylation community."

---

## 6. Honest assessment: does any tool replace FinaleMe for our use case?

**No, not in 2026.** Three independent findings support this:

1. **The only mature tool is FinaleMe itself.** Every other methylation tool
   requires bisulfite / long-read / single-cell / targeted-assay / array data
   that we don't have. The remaining alternative (WGS2meth) is brand new with
   no public code and no cfDNA validation.

2. **Even FinaleMe's published performance may not transfer.** The 0.91 auROC
   headline number was measured at 16–39× deep coverage; our cohort is 5–10×.
   The paper's own ULP-WGS (0.1×) decode results required a model trained on
   a high-coverage healthy donor — not a fresh multi-sample training — which
   is exactly what failed for us on JVM+disk constraints.

3. **The methylation-proxy result we already have (AUC 0.777) is in the same
   range as what any fragment-based methylation imputation would deliver.**
   If FinaleMe worked perfectly on our data, the realistic best-case output
   is AUC ≈ 0.78–0.85 (validated deep-coverage auROC for per-CpG calls is
   0.91, but that's per-CpG accuracy on fragments with ≥5 CpGs in CpG-rich
   regions — the aggregated per-sample cancer-vs-healthy AUC is necessarily
   lower because most of the genome is CpG-poor). And as the methylation-proxy
   head-to-head already showed, **even a 0.85 methylation AUC fused with a
   0.975 fragmentomics AUC adds ≤ +0.005 to the combined AUC** because
   fragmentomics is already near saturation.

**The honest framing for the user:** "Methylation imputation on plain cfDNA
WGS at our coverage is fundamentally a low-signal task that adds marginal
value to an already-saturated fragmentomics baseline. No public tool in 2026
overcomes this fundamental signal ceiling within a 1–2 week budget on solo
hardware. The most honest path forward is to (a) try the cheap Tier-1
FinaleToolkit deep-dive in 3–5 days, (b) investigate WGS2meth in 2–3 days,
and (c) accept the methylation-proxy result as the project's methylation
contribution if both return null."

---

## Appendix A — Quick reference for the surveyed tools

| Tool | Why surveyed | Why excluded (or note) |
|---|---|---|
| FinaleMe | Original target; failed | JVM hang + multi-sample training blocked; coverage regime concern |
| WGS2meth | Newest (2025); same biological premise as FinaleMe; simpler to run | No public code repo; no cfDNA validation |
| FinaleToolkit (cleavage ratio / end motif per CpG island) | Already installed; explicit methylation-correlated features | Not yet extracted at the per-CpG-island level — Tier-1 opportunity |
| methylKit | Common R methylation toolkit | Requires methylation calls or BS-seq, not WGS |
| Bismark / BWA-meth / MethylDackel | Standard bisulfite callers | No bisulfite → no signal |
| Bismap | Bisulfite aligner | Same |
| nanomethphase / DeepSignal / DeepMod / Nanopolish | Long-read methylation | No long reads |
| DISMIR | Per-read cancer/tissue from cfDNA | Requires WGBS |
| DeepCpG / MethylNet / Melissa / CaMelia / CpG-Transformer | Single-cell imputation | No single-cell data |
| cfTools / CelFiE / CelFiE-ISH / MethAtlas / UXM / cfNOMe toolkit | Cell-type deconvolution | Downstream only; require methylation calls from WGBS |
| cfMeDIP-seq tools / cfMethyl-Seq | Assay-specific | Wrong assay |
| EpiDISH / minfi / ChAMP / RefFreeEWAS | Methylation array | No array data |
| CelFiE-ISH ONT extension | ULP-WGS ONT deconv | Requires ONT |

---

## Appendix B — Key sources

- Liu Y et al. (2024). FinaleMe: Predicting DNA methylation by the fragmentation
  patterns of plasma cell-free DNA. *Nature Communications* 15:2790.
  https://www.nature.com/articles/s41467-024-47196-6
- Abdullaev ET, Haridoss DA, Arndt PF (2025). Predicting the methylation status
  of CpG islands from read distribution biases. *BMC Genomics* 26:973.
  https://link.springer.com/article/10.1186/s12864-025-12257-7
- Wang Y et al. (2024). FinaleToolkit: Accelerating Cell-Free DNA Fragmentation
  Analysis with a High-Speed Computational Toolkit. *Bioinformatics Advances* 5(1):vbaf236.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11160763/
- Zhou Z et al. (2022). Cleavage ratios near CpG sites — referenced by FinaleToolkit paper
  as methylation-correlated.
- Hu R et al. (2025). cfTools: an R/Bioconductor package for deconvolving
  cell-free DNA via methylation analysis. *Bioinformatics Advances* vbaf108.
  https://github.com/jasminezhoulab/cfTools
- Caggiano C et al. (2021). CelFiE. *Nat Commun* 12:2717.
- Sun K et al. (2024). Systematic evaluation of methylation-based cell type
  deconvolution methods for plasma cell-free DNA. *Genome Biology* 25:318.

---

## Appendix C — What this survey explicitly does NOT recommend

- **Do NOT spend more time on FinaleMe Java/Maven debugging** on solo hardware.
  The disk + JVM failure mode is not unique to our setup; even a working
  FinaleMe would take 2–4 weeks of compute to deliver a 627-sample matrix.
- **Do NOT acquire bisulfite data.** It's months of DAC approval and outside
  the project's scope.
- **Do NOT train a new XGBoost model from scratch on cfDNA WGS** to recreate
  WGS2meth's approach. The paper's own cell-line-only validation means the
  result may not transfer to cfDNA at all; even if it did, the per-CpG-island
  binary output is less useful than FinaleMe's continuous β-values.
- **Do NOT assume "methylation" can ever replace fragmentomics** on this data.
  Even a perfect methylation call would add ≤ +0.005 AUC to the existing
  fragmentomics baseline (cf. existing proxy result). The fragmentomics
  signal at AUC 0.975 has near-zero headroom for orthogonal channels.
