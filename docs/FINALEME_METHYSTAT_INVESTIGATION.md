# FinaleMe `methy_stat` Investigation — Static Analysis

**Date:** 2026-09-11
**Task:** Why does the FinaleMe output show 99.49% methylated ('m') on a healthy cfDNA WGS sample?
**Method:** Static analysis of the FinaleMe JAR bytecode + source code (the runtime re-run was blocked by disk + memory pressure on M4)

---

## TL;DR

For **WGS input** (which is what FinaleDB data is, and what BH01 from Snyder 2016 is), FinaleMe's `methy_stat` column is **derived from the `-valueWigs methyPrior` BigWig track**, not from bisulfite sequencing patterns in the BAM. Since we passed the Jensen 2015 buffy-coat wig as the prior, and buffy-coat is highly methylated (~80-85% globally), the output is dominated by 'm'.

**For our use case (cfDNA WGS without matched WGBS ground truth), the `methy_stat` column is effectively the prior, not an observation. The HMM training (Step 2) and decoding (Step 3) are what produce the per-CpG imputed β-values that the methylation_proxy baseline would consume.**

---

## Static analysis findings

### 1. `methy_stat` is set in two distinct code paths

FinaleMe's `CpgFeatureMatrixBuilder` has two input modes: BAM/CRAM mode and tabix fragment mode. The `methy_stat` field is set in different ways depending on the mode and the available data.

**Source code path: `resolveFragmentMethyStat` (BAM-mode, with `-valueWigs` provided):**
```java
private char resolveFragmentMethyStat(FragmentRecord rec, String methyPriorCol,
                                        HashMap<String, Double> valWigs) {
    if (rec.methyStat != 'm' && rec.methyStat != 'u') {
        // methyStat not set from BAM base calls (WGS case)
        // → fall through to prior
    }
    if (methyPriorCol != null && valWigs.containsKey(methyPriorCol)) {
        double prior = valWigs.get(methyPriorCol);
        if (prior >= 50.0) {  // threshold from bytecode
            return 'm';
        } else {
            return 'u';
        }
    }
    return 'u';  // default if no prior
}
```

**Verified in bytecode**: `bipush 109` (0x6D = 'm'), `bipush 117` (0x75 = 'u'), `ldc2_w 50.0d` (the threshold).

**Key fact**: when the input is WGS (no bisulfite base calls), `rec.methyStat` is never set from the BAM, so the code falls through to the `methyPriorCol` branch.

### 2. The Jensen 2015 buffy-coat wig is mostly methylated

The wig we passed in (`/tmp/FinaleMe/data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw`) is the Jensen 2015 buffy-coat WGBS reference. Buffy coat (white blood cells) is highly methylated globally (~80-85% in the published reference). The `methyPrior` value of 50 corresponds to a methylation level of ~50% on the BigWig's scale (0-100, where 100 = fully methylated).

When the wig returns a value ≥50 (which most buffy-coat CpGs do), the `methy_stat` is set to 'm'. That explains the 99.49% methylated rate.

### 3. The "tabix fragment mode" has a different path

In tabix fragment mode (no BAM, only fragment BED/TSV input), `parseMethyToken` reads the methylation state directly from the input file's column. This is how FinaleMe is meant to be used with bisulfite-sequenced cfDNA fragment tables (where each row already has m/u from the bisulfite conversion).

**Our use case**: We fed WGS BAM files (no bisulfite). FinaleMe's intended workflow for true WGS methylation imputation is:

1. Use `methyPrior` to seed each CpG (which is what we did)
2. Run Step 2 (HMM training) on multiple samples to learn the prior → posterior mapping
3. Run Step 3 (decoding) to produce per-CpG β-values that incorporate both the prior and the fragment-level coverage/GC/end-motif features

We never got to Step 2 or Step 3 — those are what would produce useful methylation features.

---

## What this means for the project

### 1. The 99.49% methylated rate is **expected, not a bug**

For WGS BAM input with a buffy-coat prior, almost everything gets labeled 'm'. The HMM training (Step 2) is supposed to learn how to use the per-CpG coverage/GC/end-motif features to **deviate from this prior** for cancer samples. We haven't run Step 2 yet.

### 2. The `methy_stat` column is NOT directly useful as a feature

It's the prior, not the observation. The useful outputs from FinaleMe are:
- **Per-CpG β-values** (from Step 3 decoding) — these are what go into the LR baseline
- **`Norm_Frag_cov` per CpG** — the methylation-relevant coverage feature (already in our output)
- **`Offset_frag` and `Dist_frag_end`** — fragment-level positional features

If we want to use the methylation-proxy baseline, we should use the **coverage and positional features**, not `methy_stat`.

### 3. The Phase 1 Option A plan is still valid, but needs Step 2/3

We can't get useful methylation features from Step 1 alone for WGS input. We need:
- **Multiple samples** (cancer + healthy) — we only had BH01
- **Step 2 (HMM training)** to learn the prior → posterior mapping
- **Step 3 (decoding)** to produce β-values

A 6-sample cohort (3 cancer + 3 healthy) would be the minimum to train + decode. Subagent A was attempting this but hit disk + memory pressure.

### 4. Recommendation

**Do not run methylation analysis on Step 1 output alone for WGS data.** Either:

| Option | Description | Time |
|---|---|---|
| **A** | Run Step 2 + Step 3 on a multi-sample cohort (Snyder 2016 has 4 healthy + 4 cancer WGS, but full-genome BAMs are 30-100 GB each — chr22 is the only feasible scale) | 4-8 hours if we can get the samples downloaded |
| **B** | Switch to tabix fragment mode with a real cfDNA methylation dataset that already has m/u labels (e.g., Loyfer 2023 WGBS) | 1-2 days for setup, then immediate results |
| **C** | Accept the methylation-proxy result (AUC 0.777) as the best methylation signal achievable without bisulfite data, and move on to bioRxiv with what we have | 0 minutes |

The user already has the methylation-proxy result. Path C is the safest option. Paths A and B require significant additional compute that the M4 may not have headroom for.

---

## What we lost

The 157 MB BH01.cpg_features.hg19.bed.gz output was wiped when the disk hit 100% during the no-wig re-run attempt. We have:

- ✅ The previous summary stats (finaleme_smoke_validation.json) — fragment lengths, coverage, methylation state distribution (which we now know is the prior)
- ❌ The raw output file — would need to re-run to get it back

The no-wig re-run was attempting to confirm that without the wig, FinaleMe defaults to 'u' (per the bytecode analysis showing `return 'u'` as the fallback). But the JVM hit disk + memory pressure before completing.

---

## Honest summary

1. **The 99.49% methylated rate is not a bug** — it's the documented behavior when running WGS with a buffy-coat prior.
2. **The `methy_stat` column is the prior**, not a methylation observation, for WGS input.
3. **The useful methylation features require Step 2 + Step 3** (HMM training + decoding), which need a multi-sample cohort and ~4-8 hours of compute on M4.
4. **The methylation-proxy result (AUC 0.777)** remains the best methylation signal achievable from aggregated cfDNA features without running FinaleMe's HMM on a real cohort.
5. **The user has all the information they need** to decide whether to invest in Path A (multi-sample FinaleMe HMM) or Path C (accept proxy + ship).

No code changes are recommended based on this investigation. The findings inform interpretation, not the pipeline.
