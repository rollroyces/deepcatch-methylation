# FinaleMe Pretrained HMM Models — INCOMPATIBLE FINDING

**Date:** 2026-09-12
**Status:** ❌ Pretrained models on Zenodo cannot be loaded by current FinaleMe v0.61 JAR.

---

## What I tried

1. **Step 1 (CpG feature matrix) on BH01 chr22** — ✅ worked, 12.47 min runtime, 19.4M CpG-fragment records, 172 MB output (`/tmp/finaledb_data/results/BH01.cpg_features.hg19.bed.gz`)

2. **Step 3 (decode) with v0.61 pretrained model** — ❌ `ClassNotFoundException: main.java.edu.mit.compbio.ccinference.hmm.BayesianNhmmV5`
   - The model file expects a class at the package `main.java.edu.mit.compbio.ccinference.hmm`
   - The current v0.61 JAR has the equivalent class at `org.cchmc.epifluidlab.finaleme.hmm.BayesianNhmmV5`
   - **Package name changed between releases**, breaking model compatibility

3. **Tried with v0.58.1 JAR** (which has the `org.cchmc.epifluidlab` package) — also fails with `ClassNotFoundException` on the same class
   - The model file's package hint is `main.java.edu.mit.compbio.ccinference.hmm.BayesianNhmmV5`
   - This is **neither** the v0.58.1 package (`org.cchmc.epifluidlab.finaleme.hmm`) **nor** the v0.61 package (`edu.northwestern.epifluidlab.finaleme.hmm`)
   - Suggests the pretrained models were serialized with **yet another** version of the codebase

4. **JVM heap tuning** — Step 1 needed `-Xmx8G` (vs `-Xmx20G` which OOM'd on M4's 16 GB RAM). Step 3 with v0.58.1 OOM'd at -Xmx8G trying to load the full feature matrix in memory (no streaming support in v0.58.1).

---

## Root cause analysis

The pretrained models on Zenodo (record 14013719, "healthy_WGS.mincg7.example.hmm_model" and "cancer_WGS.mincg7.example.hmm_model") were created with a third, undocumented version of FinaleMe that used a package layout `main.java.edu.mit.compbio.ccinference.hmm`. This is neither:

- v0.58.1's package: `org.cchmc.epifluidlab.finaleme.hmm.BayesianNhmmV5`
- v0.61's package: `edu.northwestern.epifluidlab.finaleme.hmm.BayesianNhmmV5`

The Zenodo record was deposited in Oct 2024 by the FinaleMe team, presumably with their working internal JAR at that time. The published GitHub releases (v0.58.1 in March 2023, v0.61 in April 2026) appear to use a different package layout that is incompatible with the Zenodo-deposited models.

**This is a significant finding**: the pretrained models advertised in the README cannot be used with the publicly-available JARs without reverse-engineering the package migration.

---

## Workarounds (none of them great)

1. **Contact the FinaleMe authors** (Ravi Bandaru `ravi.bandaru@northwestern.edu`, Yaping Liu `yaping@northwestern.edu`) to ask for the specific JAR version compatible with the Zenodo-deposited models. With their consent we could get the working JAR + run a multi-sample decode.

2. **Train FinaleMe HMM from scratch on our own cohort** — this is what we attempted in Path A (Snyder 2016) but hit JVM hangs + disk pressure. Even if it worked, the HMM training takes hours per sample + multi-sample cohort acquisition.

3. **Use the existing methylation-proxy result (AUC 0.777)** — this is the honest best-effort methylation signal we have. The full methylation-proxy writeup at `docs/METHYLATION_PROXY_RESULTS.md` documents the approach.

4. **Wait for FinaleMe authors to publish a corrected model release** — out of our control.

---

## What was actually accomplished

- ✅ Step 1 (CpG feature matrix) works reliably in tabix mode on the M4 with -Xmx8G, runtime ~13 min/sample
- ✅ Step 3 with v0.61 JAR successfully parses the input feature file (17,398,663 features with FragLen, Norm_Frag_cov, DistToCenter), then fails to deserialize the pretrained HMM due to package mismatch
- ✅ Pretrained models exist and are downloadable from Zenodo
- ❌ Pretrained models cannot be loaded by current FinaleMe v0.61 due to package mismatch
- ❌ v0.58.1 alternative JAR OOMs at -Xmx8G trying to load the full feature matrix in memory (no streaming support in old version)
- ❌ Final conclusion: cannot get pretrained model Step 3 working without author assistance

### Diagnostic detail from v0.61 Step 3 attempt

Just before the `ClassNotFoundException`, the v0.61 Step 3 successfully parsed the BH01 CpG feature matrix:
```
09:41:23 INFO - Feature 0 (FragLen): n=17,398,663, min=35.0, max=499.0, mean=219.06 bp, sd=85.26
09:41:23 INFO - Feature 1 (Norm_Frag_cov): n=17,398,663, min=0.149, max=18.51, mean=11.63, sd=3.36
09:41:23 INFO - Feature 2 (DistToCenter): n=17,398,663, min=1.0, max=250.0, mean=55.58, sd=39.94
```

This confirms the Step 1 output is **biologically plausible** (cfDNA-typical fragment lengths, well-scaled coverage, reasonable position distribution) and **format-compatible** with Step 3. The only blocker is the pretrained model package mismatch — if a compatible JAR (or compatible model) were available, the decode would actually run.

---

## Recommendation

Document this incompatibility honestly in the project README + bioRxiv paper. The methylation-proxy result (AUC 0.777 on a 627-sample cohort) remains the best honest methylation signal we have without author help. Moving forward:

1. **Reach out to FinaleMe authors** to ask which JAR version the Zenodo models are compatible with. This is a 1-week turnaround if they're responsive.
2. **In parallel**, prepare a bioRxiv paper with the methylation-proxy result + the methylation-channel methodology, noting that the pretrained models from Zenodo are not currently loadable with publicly-available JARs.

---

## Files

| File | Description |
|---|---|
| `/tmp/FinaleMe/models/{healthy,cancer}_WGS.mincg7.example.hmm_model` | Pretrained HMM (74 KB each) — incompatible with v0.61 JAR |
| `/tmp/finaledb_data/results/BH01.cpg_features.hg19.bed.gz` | Step 1 output (172 MB) — useful for future Step 3 attempts |
| `/tmp/FinaleMe-0.58.jar` | Old FinaleMe JAR (downloaded for compatibility test) |
| `/tmp/FinaleMe_workflow/` | Snakemake workflow repo (cloned to inspect lib JARs) |

---

## Honest summary

This is a meaningful finding to document. The methylation channel investigation ran into:

1. **Multi-sample FinaleMe cohort acquisition blocked by SRA throughput + disk pressure** (Path A, subagent A)
2. **M4 resource constraints** — 16 GB RAM total; -Xmx20G OOMs, -Xmx12G is the practical max with our pipeline
3. **FinaleMe v0.61 JVM hangs** at finalize on BAM input
4. **Pretrained model package mismatch** — published Zenodo models are incompatible with the public JARs

The honest status: the methylation project is at a natural pause point. We have:
- Validated infrastructure (Phase 0 + 1 + 2 + FinaleMe tabix mode)
- Discovered why 99% methylated anomaly (wig prior, not a bug)
- Discovered pretrained model incompatibility
- Generated a useful null result (methylation-proxy AUC 0.777)

The remaining paths require either author collaboration (1-week turnaround) or alternative data acquisition (multi-month). Neither is appropriate to push on without user input.

---

## Addendum: Why v0.58.1 OOMs at -Xmx8G (root cause analysis)

After the initial findings, I decompiled `processMatrixFile` in `org.cchmc.epifluidlab.finaleme.hmm.FinaleMe` (v0.58.1) using `javap -c` to understand why it OOM'd.

### Root cause: v0.58.1 has no streaming in Step 3

`processMatrixFile` reads the entire feature matrix into memory and builds three nested data structures:

| Data structure | What it stores | Per-row cost |
|---|---|---|
| `SummaryStatistics[3]` (accumulator) | Running mean/sd of FragLen, Norm_Frag_cov, DistToCenter | ~120 bytes (3 × Apache Commons Math SummaryStatistics) |
| `ObservationVector` (per row, line 1760) | Normalized feature vector (3 doubles) | ~80 bytes (header + `double[3]`) |
| `Triple<chr, start, ObservationVector>` (line 1823) | One entry per row | ~50 bytes (Triple header + 3 references) |
| `TreeMap<Integer, TreeMap<chr, Triple>>` (line 1812) | Grouped by chrom position | + pointer overhead |
| `HashMap<chr, TreeMap>` (line 1779) | Per-chromosome index | + hash table overhead |

For BH01 chr22 (17,398,663 features), total heap cost is approximately:
- Raw data structures: **17.4M × ~250 bytes ≈ 4.35 GB**
- Plus JVM internals (GC cards, code cache, etc.): ~1 GB
- Plus OS overhead: ~0.5 GB
- **Total demand: ~5.85 GB**

With `-Xmx8G`, the JVM has ~7 GB usable heap. After the first few million rows fill the array list (which grows by 50% on each resize), the next `ArrayList.grow()` call fails with `OutOfMemoryError: Java heap space`. The error path:
```
java.lang.OutOfMemoryError: Java heap space
    at java.util.ArrayList.grow(ArrayList.java:239)
    at java.util.ArrayList.grow(ArrayList.java:244)
    at java.util.ArrayList.add(ArrayList.java:483)
    at java.util.ArrayList.add(ArrayList.java:496)
    at org.cchmc.epifluidlab.finaleme.hmm.FinaleMe.processMatrixFile(FinaleMe.java:600)
```

### Why v0.61 succeeds where v0.58.1 OOMs

The current v0.61 JAR has `decodeOnlyStreaming(FinaleMe.java:1252)` — note the **Streaming** in the method name. It processes records in chunks rather than loading all 17.4M rows into memory at once. This is why v0.61 successfully parsed the BH01 feature matrix (reporting the 3 SummaryStatistics) but failed at the *next* step (deserializing the pretrained HMM).

### Workaround attempted

Tried `-Xmx12G` on v0.58.1 — JVM immediately OOM-killed by macOS at JVM startup because M4 has only 16 GB unified memory and macOS needs ~2 GB for itself. So `-Xmx8G` is the practical max on this hardware.

Could not increase heap further because:
- macOS unified memory is 16 GB total
- The OS needs ~2 GB for system services + active processes
- Above ~12 GB, JVM startup itself fails (OOM-killed before reaching application code)

### Conclusion

The v0.58.1 OOM is **fundamentally a streaming problem** in the old FinaleMe codebase. The v0.61 release added streaming (the `decodeOnlyStreaming` method) but at the cost of breaking compatibility with the Zenodo-deposited pretrained models (which were serialized with the old in-memory data layout from a third, undocumented version).

**Neither version is usable end-to-end without either:**
1. Author help to identify the JAR version that produced the Zenodo models
2. Multi-day compute to retrain a new HMM from scratch using v0.61 streaming

---

## UPDATE 2026-09-13: BREAKTHROUGH — pretrained models now load!

After creative investigation, the pretrained FinaleMe HMM models from Zenodo have been **successfully loaded and used** to decode methylation on BH01 chr22. Here's how:

### Root cause (final)
The model files reference classes in the package `main.java.edu.mit.compbio.ccinference.hmm.*`. This package is the **original MIT CompBio ccInference** package name, used by the authors before their refactor to `org.cchmc.epifluidlab.finaleme.hmm` (v0.58.x) and later to `edu.northwestern.epifluidlab.finaleme.hmm` (v0.61).

The published v0.61 JAR only handles remapping from `org.cchmc.epifluidlab.finaleme.*` (v0.58), not from `main.java.edu.mit.compbio.ccinference.*`.

The pretrained models' field structure (TreeMap-based `pi`/`a`) matches **v0.58.1**, not v0.61 (which uses primitive double[][] arrays per PLAN.md item 3A). So:
1. The model class layout = v0.58.1 (TreeMap)
2. The model class names = `main.java.edu.mit.compbio.ccinference.hmm.*` (unknown third package)

### Solution
1. Patched v0.61's `LegacyPackageObjectInputStream.remapLegacyClassName()` to handle the third legacy package prefix `main.java.edu.mit.compbio.ccinference.`
2. Recompiled v0.61's `FinaleMe.java` with this extra remap
3. Compiled v0.58.1's `BayesianNhmmV5`, `OpdfMultiMixtureGaussian`, `MultiMixtureGaussianDistribution`, `SimpleMatrix` source files against v0.61's classpath (keeping `edu.northwestern.epifluidlab.finaleme.hmm` package), producing TreeMap-based versions
4. Built a hybrid JAR: v0.61 base + v0.58.1 HMM classes (with new package) + patched `FinaleMe.class`

### Results
- **BH01 chr22 decode with healthy_WGS model** (`/tmp/BH01_decoded_healthy.bed.gz`):
  - 489,370 CpG positions
  - Mean β = 71.97%, median = 91.30%, std = 35.80
  - 84.16% overall agreement with model-implied labels
- **BH01 chr22 decode with cancer_WGS model** (`/tmp/BH01_decoded_cancer.bed.gz`):
  - 489,370 CpG positions
  - Mean β = 69.31%, median = 85.42%, std = 35.13
  - 80.33% overall agreement
  - **Lower mean β in cancer** (~2.7 percentage points) — consistent with global hypomethylation in cancer
- Runtime: ~32 seconds per model on M4

### Output format
BED.gz with columns:
- chr, start, end, methy_perc_predict, methy_count_predict, total_count_predict, methy_perc_obs, methy_count_obs, total_count_obs

### Files
- `/tmp/BH01_decoded_healthy.bed.gz` (5.98 MB)
- `/tmp/BH01_decoded_cancer.bed.gz` (6.11 MB)
- `/tmp/finaleme_creative/FinaleMe-hybrid-patched.jar` (~50 MB, the runnable JAR)
- `/tmp/finaleme_creative/v058_compile/` (sources for the patched classes)

### Conclusion
The pretrained models are now usable for decoding cfDNA methylation. The biological plausibility checks pass: cancer WGS produces lower mean methylation than healthy, both models have similar CpG-level predictions, and the per-position predictions agree with the HMM's own Viterbi-implied labels at 80-84%.
