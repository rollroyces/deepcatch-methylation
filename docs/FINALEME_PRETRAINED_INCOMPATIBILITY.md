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
- ✅ Pretrained models exist and are downloadable from Zenodo
- ❌ Pretrained models cannot be loaded by current FinaleMe v0.61 due to package mismatch
- ❌ v0.58.1 alternative JAR OOMs at -Xmx8G and has different package structure
- ❌ Final conclusion: cannot get pretrained model Step 3 working without author assistance

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
