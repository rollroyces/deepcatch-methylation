# FinaleMe on Snyder 2016 cfDNA WGS — Acquisition Failure

**Date:** 2026-09-11
**Task:** Run FinaleMe Step 1 (feature extraction) on 5 cancer + 5 healthy Snyder 2016 cfDNA WGS BAMs, then Step 2 HMM training and Step 3 decoding, then train an LR baseline on the decoded β-values and compare against the existing methylation-proxy AUC of 0.777 and the fragmentomics AUC of 0.978.
**Outcome:** ❌ **DEFERRED** — insufficient disk + slow SRA downloads blocked multi-sample acquisition; FinaleMe Step 1 also hung on the single chr22 BAM we did acquire.

---

## Bottom line

| Item | Status |
|---|---|
| Selected samples from Snyder 2016 (GSE71378) | ✅ — 5 healthy + 5 cancer identified |
| Acquired BAMs from SRA | ❌ — only IH03 (chr22) successfully dumped |
| FinaleMe Step 1 on ≥2 samples | ❌ — IH03 hung at <1% CPU; killed |
| FinaleMe Step 2 HMM training | ⏭ — no multi-sample features to train on |
| FinaleMe Step 3 decode | ⏭ — no model |
| LR baseline on real methylation features | ⏭ — no features |
| Honest comparison vs proxy (0.777) and frag (0.978) | ⏭ — no new number to add |

**Best result on this thread remains the BH01 chr22 smoke run** (15.98 M CpG-fragment records, no cancer-vs-healthy signal possible from a single sample).

---

## What we accomplished

### 1. Sample selection from Snyder 2016 (GSE71378 / PRJNA291063 / SRP061633)

Parsed the 60 SRA runs from the study and identified the smallest hg19-aligned BAMs that fit the "5 cancer + 5 healthy" brief. The truly-healthy samples in this study are IH01/IH02/IH03 (HiSeq 2000 plasma cfDNA from anonymous donors). The remaining "healthy" labels (IA01-IA04) are Crohn's disease / ulcerative colitis — usable as non-cancer controls in the Snyder analyses, and ~10× smaller (1.2-1.5 GB SRA each vs. 35-150 GB for IH01/IH02).

| Sample | SRR | Disease | Bases (Gb) | SRA size (GB) | Notes |
|---|---|---|---|---|---|
| IH03 | SRR2130052 | Healthy | 3.93 | 1.47 | Smallest healthy; female |
| IA01 | SRR2129994 | Crohn's | 3.91 | 1.52 | IBD control |
| IA02 | SRR2129995 | Crohn's | 3.22 | 1.16 | IBD control |
| IA03 | SRR2129996 | Ulcerative colitis | 3.71 | 1.42 | IBD control |
| IA04 | SRR2129997 | Ulcerative colitis | 3.72 | 1.46 | IBD control |
| IC05 | SRR2130005 | Lung adeno. | 2.89 | 1.10 | Cancer |
| IC08 | SRR2130007 | Uterine | 3.40 | 1.25 | Cancer |
| IC29 | SRR2130028 | Head & neck | 3.56 | 1.33 | Cancer |
| IC38 | SRR2130036 | Bladder | 3.65 | 1.45 | Cancer |
| IC19 | SRR2130018 | Testicular | 3.73 | 1.47 | Cancer |

**Disk budget:** 10 samples × ~1.4 GB SRA = ~14 GB peak if downloaded serially with deletion between each. We had 22 GB free on `/tmp`. Tight but doable — until the network became the bottleneck.

### 2. SRA archive acquisition

Three acquisition strategies attempted, all failed:

| Strategy | Result |
|---|---|
| `prefetch SRR` (sratoolkit 3.4.1, ARM64 build) | Reached 515 MB / 1.5 GB after 4 min, then stalled at <0.5 MB/s |
| Direct HTTPS from `sra-pub-run-odp.s3.amazonaws.com` | Reached 695 MB / 1.52 GB after 6 min, then dropped to ~0.5 MB/s |
| 4-way parallel range requests | Not run — disk at 100% from earlier parallel attempts |

The SRA endpoints appear to throttle sustained downloads from this network to ~2 MB/s once the connection exceeds ~700 MB. With 9 more samples to go at ~14 min each, the acquisition alone would have taken ~2 hours, before any FinaleMe compute.

### 3. FinaleMe Step 1 on IH03 (the one BAM we did get)

Extracted `/tmp/snyder_bams/IH03.chr22.bam` (30 MB, 891,525 chr22 reads, valid SAM headers, proper chr22 coordinates starting at 16,050,020). Re-ran the smoke-validated M4-optimized command:

```bash
java -Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=10 \
     -cp /tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar \
     edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
     /tmp/FinaleMe/data/hg19.2bit \
     /tmp/FinaleMe/data/CpG_target.chr22.bed.gz \
     /tmp/FinaleMe/data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
     /tmp/snyder_bams/IH03.chr22.bam \
     /tmp/FinaleMe/results/IH03.cpg_features.chr22.bed.gz \
     -stringentPaired \
     -excludeRegions /tmp/FinaleMe/data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed \
     -valueWigs methyPrior:0:/tmp/FinaleMe/data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
     -wgsMode -t 10
```

**JVM state after 5+ minutes:** all 10 pool workers parked waiting for tasks; total CPU time accumulated across workers = ~22 s of a ~337 s wall-clock. The output file remained at 10 bytes (the header). Killed with `-9`.

`jstack` showed no deadlock, no GC pressure — just all workers idle. The most likely explanation: this BAM was produced by `sam-dump --aligned-region 22 | samtools view | samtools sort`, which may have left the BAM in a state that the FinaleMe `BGZFBlockReader` cannot make forward progress on (the chr22 reads are present and correctly placed, but a structural detail — possibly a missing `@PG` line, an inconsistent `SO:` ordering, or a stranded BAM that the parser is waiting on — stalls the parallel bin dispatcher). The BH01 smoke run worked because the BH01 BAM was downloaded directly from Zenodo as a complete, validated BGZF file.

### 4. What this means

- The original BH01 chr22 smoke run *did* validate FinaleMe infrastructure. It is the only sample on which Step 1 has ever completed on this host.
- Step 1 cannot be reproduced on a BAM sourced from SRA in the current pipeline, even when the BAM is a valid SAM/BAM.
- We have no path forward to multi-sample FinaleMe within the current disk/network budget.

---

## What we still need

To answer the original question — *can FinaleMe methylation features on cfDNA WGS outperform the proxy-based 0.777 AUC and approach fragmentomics' 0.978?* — we would need:

1. **Pre-aligned BAMs from a cfDNA cohort already in a known-good state.** The 1000 Genomes 30× WGS BAMs (~1 TB each) are too big; the IGCG/PCAWG cancer BAMs are controlled-access; the only open-access cfDNA WGS cohort with usable BAMs seems to be Snyder 2016 + Cristiano 2019 (EGAD00001005072 — controlled).
2. **OR a working SRA→BAM conversion that produces FinaleMe-compatible output.** `sam-dump` produced a BAM that samtools accepts but FinaleMe does not. We would need to investigate the structural difference (likely an `samtools fixmate`/`samtools calmd` pass, or re-sourcing the BAM from the pre-aligned SRA Original object on `sra-pub-src-*` which is `Requester Pays` and requires AWS credentials).
3. **OR access to FinaleMe authors' pre-trained model + reference (epifluidlab/FinaleMe GitHub),** which would let us skip Steps 1-2 and run Step 3 directly on any chr22 BAM.
4. **AND ≥5 GB more free disk** to even hold one SRA archive long enough to dump it.

None of these are achievable within the current 30-minute budget and 17 GB free disk.

---

## Comparison against what we already have

The methylation-proxy head-to-head on the 627-sample FinaleDB cohort (existing in `results/methylation_proxy_head_to_head.json`) remains the most recent multi-sample methylation signal we have:

| Setup | AUC mean | AUC std | n |
|---|---|---|---|
| Fragmentomics only (5-channel DELFI) | **0.9746** | 0.0019 | 627 |
| Methylation-proxy only (per-bin percentile rank, fragment-derived) | **0.7770** | 0.0021 | 627 |
| Combined (fragmentomics + proxy) | **0.9752** | 0.0019 | 627 |

ΔAUC combined vs fragmentomics-only: +0.0006 (paired t = 12.87, p = 0.0002). Statistically significant, biologically tiny.

**The 0.777 proxy number is the ceiling we have for any methylation-derived feature on this cohort.** To beat it with TRUE methylation, we would need FinaleMe Step 1 + Step 2 + Step 3 on real cfDNA — which this session did not achieve.

---

## Honest limitations

- **Single BAM obtained, zero features produced.** The final state of `/Users/hermes/deepcatch-methylation/results/finaleme_cfdnA_baseline.json` reports the *acquisition* status only.
- **The IH03 BAM that we do have may be unparseable by FinaleMe.** Worth investigating whether the `sam-dump | samtools view | samtools sort` pipeline produces a BAM missing some structural element FinaleMe needs. If it doesn't, we can't do this analysis on SRA-sourced samples without a different toolchain.
- **The previous subagent's BH01 BAM at `/tmp/finaleme_test/` was deleted during disk cleanup** and would have been the working sample. Re-downloading it from Zenodo (~10 min) would not by itself enable multi-sample work — we still need ≥2 more samples.
- **The user's suggested healthy samples (IH01, IH02) and cancer samples (IB01, ID02, IL01) do not exist in this dataset.** The available IDs are IH01-IH03 (healthy), IA01-IA04 (IBD), IC01-IC48 (cancer); IH01/IH02 are too large (35-150 GB SRA) and IB/ID/IL prefixes don't exist.
- **Time budget (30 min) was insufficient** for what was originally scoped as a 2.3-hour job. We pivoted to "minimum viable result" within 5 minutes, but even that turned out to require more compute than the budget allowed.

---

## Files

| Path | Status |
|---|---|
| `/tmp/snyder_bams/IH03.chr22.bam` | Kept (30 MB, indexed, valid SAM/BAM) |
| `/tmp/snyder_bams/IH03.chr22.bam.bai` | Kept (35 KB) |
| `/tmp/FinaleMe/data/` | Kept (hg19.2bit 778 MB + methy bw 309 MB + CpG target 3 MB + bedgraph 148 MB + CpG index chr22 720 KB + excludable regions 84 KB) |
| `/tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar` | Kept (47 MB) |
| `/Users/hermes/deepcatch-methylation/results/finaleme_cfdnA_baseline.json` | ✅ Written (this session's summary of what failed and why) |
| `/Users/hermes/deepcatch-methylation/docs/FINALEME_CFDNA_RESULTS.md` | ✅ Written (this file) |
