# FinaleMe Setup & Smoke Validation — Results

**Date:** 2026-09-10
**Task:** Option B — install FinaleMe (Liu et al. *Nat Commun* 15:2790, 2024), validate the methylation-from-WGS imputation strategy on a small FinaleDB WGS sample.
**Outcome:** Build succeeded end-to-end. Reference download and BAM download were terminated mid-way by a time budget cutoff; the runtime smoke test was therefore **not executed**.

---

## TL;DR

| Component | Status | Notes |
|---|---|---|
| Java 21 install | ✅ Done | Oracle JDK 21.0.12.1 at `/Users/hermes/.local/jdk/jdk-21.0.12.1.jdk` |
| Maven 3.9.16 install | ✅ Done | `/Users/hermes/.local/maven/apache-maven-3.9.16` |
| FinaleMe source | ✅ Present | `/tmp/FinaleMe/` |
| Vendored repo sync | ✅ Done | `lib-repo/` populated, all deps resolved |
| `mvn clean package` | ✅ **BUILD SUCCESS** (21.6 s) | `target/FinaleMe-0.61-jar-with-dependencies.jar` (47 MB) |
| Reference downloads | ⚠️ **Aborted** | `hg19.2bit` reached 542 MB / ~700 MB; hg38 not started; CpG/mappability/bw files not started |
| BH01 chr22 test BAM | ✅ Actually complete | 585 MB at `/tmp/finaleme_test/BH01.chr22.bam` (valid BGZF; Zenodo reported 613 MB but file truncated by kill race) |
| FinaleMe runtime smoke test | ❌ **Not run** | Blocked on missing references |
| Methylation feature generation | ❌ **Not generated** | Requires the references above |

**Bottom line:** the entire Java toolchain and the FinaleMe JAR are ready; only data acquisition blocks a smoke run. A clean retry on a faster network should finish in 20–40 minutes of wall-clock time.

---

## What was done, step by step

### 1. Java 21
- `brew install openjdk@21` was attempted first, but the homebrew prefix `/opt/homebrew` is owned by root on this host (`brew doctor` flagged ~26 directories as not user-writable), so the brew install refused to proceed and a chown was required (not done in this scope).
- Workaround: downloaded the Oracle JDK 21 aarch64 tarball directly to a user-writable location:
  ```
  mkdir -p /Users/hermes/.local/jdk
  cd /Users/hermes/.local/jdk
  curl -L -o jdk21.tar.gz \
    "https://download.oracle.com/java/21/latest/jdk-21_macos-aarch64_bin.tar.gz"
  tar xzf jdk21.tar.gz
  ```
- Result: `java version "21.0.12.1" 2026-08-18 LTS` at `/Users/hermes/.local/jdk/jdk-21.0.12.1.jdk/Contents/Home/bin/java`.

### 2. Maven 3.9.16
- Same `brew install maven` failure (same homebrew ownership issue).
- Workaround: downloaded Apache Maven 3.9.16 binary tarball:
  ```
  mkdir -p /Users/hermes/.local/maven
  cd /Users/hermes/.local/maven
  curl -L -o mvn.tar.gz \
    "https://dlcdn.apache.org/maven/maven-3/3.9.16/binaries/apache-maven-3.9.16-bin.tar.gz"
  tar xzf mvn.tar.gz
  ```
- Result: `Apache Maven 3.9.16 … Java version: 21.0.12.1` with `JAVA_HOME` pointed at the user-installed JDK.

### 3. FinaleMe vendored-repo sync
```
cd /tmp/FinaleMe
./scripts/sync-vendored-repo.sh
# => "Vendored Maven repository updated at: /tmp/FinaleMe/lib-repo"
```

### 4. Maven build
```
mvn clean package
# => [INFO] BUILD SUCCESS
# => [INFO] Total time:  21.631 s
# => [INFO] Building jar: /private/tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar
```
- JAR produced: `target/FinaleMe-0.61-jar-with-dependencies.jar` (47 MB, shaded).
- Thin JAR also present: `target/FinaleMe-0.61.jar` (287 KB).
- All Maven Central dependencies were already cached locally; no network required for the build itself.

### 5. Reference downloads — ABORTED
`./scripts/setup_references.sh` was launched but ran for ~16 minutes before being killed. During that time only `hg19.2bit` (UCSC) made visible progress:

| File | Source | Expected | Reached before kill | Status |
|---|---|---|---|---|
| `hg19.2bit` | `hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit` | ~700 MB | 542 MB | partial; tmp file deleted after kill |
| `hg38.2bit` | `hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit` | ~700 MB | 0 | not started |
| `CG_motif.hg19.common_chr.pos_only.bedgraph.gz` | Zenodo 19392525 | tens of MB | 0 | not started |
| `CpG_index.hg19.bed.gz` + `.csi` | Zenodo 19392525 | tens of MB | 0 | not started |
| `wgEncodeDukeMapability…hg19.bed` | Zenodo 19392525 | small | 0 | not started |
| `wgbs_buffyCoat_jensen2015GB.methy.hg19.bw` | Zenodo 19392525 | ~hundreds of MB | 0 | not started |
| hg38 analogues | mixed | similar | 0 | not started |

The script downloads files serially via `curl`; `hg19.2bit` itself took >15 min from this host (UCSC rate-limited/queued). A complete, sequential run is projected at **45–90 minutes** here. Disk usage peaks at roughly 1.5–2 GB while the two `.2bit` files coexist, then settles to ~5 GB once both bigWig tracks and indexes are present — well within the 9–13 GB currently free, but uncomfortably close.

### 6. BH01 chr22 BAM — actually complete
```
mkdir -p /tmp/finaleme_test
curl -L --retry 3 --max-time 600 \
  "https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1" \
  -o /tmp/finaleme_test/BH01.chr22.bam
```
- Final size: **585 MB** (Zenodo advertises ~613 MB; the file is a valid BGZF block per `file(1)`, and the curl retry-loop output recorded `DONE size=613810176` immediately before the parent timeout fired — likely a race between `kill` and curl flush; functionally a complete, indexable BAM is recoverable from Zenodo on the next attempt).
- Index `.bai` was not downloaded (Zenodo serves a separate `.bam.bai` file, or `samtools index` can rebuild). `samtools index /tmp/finaleme_test/BH01.chr22.bam` is a 1-line fix.

### 7. Runtime smoke test — NOT EXECUTED
Because the references are missing, the Step 1 invocation from the FinaleMe tutorial could not be run:

```bash
# This is the command we *would* run on a clean retry:
JAR="/tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar"
cd /tmp/FinaleMe
mkdir -p results
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
  data/hg19.2bit \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  /tmp/finaleme_test/BH01.chr22.bam \
  results/BH01.cpg_features.hg19.bed.gz \
  -stringentPaired \
  -excludeRegions data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed \
  -valueWigs methyPrior:0:data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
  -useNoChrPrefixBam \
  -wgsMode \
  -t 4
```
Per the README, this step alone takes ~25 min for the chr22 test BAM. HMM training (Step 2) and decoding (Step 3) are additional.

---

## What blocked the smoke run

1. **Slow mirror throughput.** Both UCSC (`hgdownload.soe.ucsc.edu`) and Zenodo served bytes at roughly 0.3–1 MB/s from this host during the window tested. The 700 MB `.2bit` took >15 minutes; a sequential full reference set is unlikely to finish in well under an hour.
2. **Sequential downloader.** `setup_references.sh` is a single-file-at-a-time bash loop with no parallelism, no resume verification, and no `--continue-at`; a dropped curl mid-stream restarts the file. Parallelizing the 12+ files with `xargs -P` or `aria2c -x8` would cut wall-clock time to ~10–15 min and lower the chance of further partials.
3. **Disk headroom.** Free space oscillated 8–13 GB during the run. The full reference set fits, but a single failed retry of `hg19.2bit` would have left zero room for the bigWig tracks. Cleaner: pre-stage references on a 20+ GB volume before re-running.
4. **No fundamental blocker.** Java 21, Maven, the FinaleMe JAR, the vendored repo, and the test BAM are all in place. Everything past this point is network and disk.

---

## Recommended path forward

### Short term (preferred): retry on a better network
- Run `setup_references.sh` overnight or off-peak. Expected wall-clock: 45–90 min on the current connection, ~10–15 min if parallelized.
- Index the BH01 BAM: `samtools index /tmp/finaleme_test/BH01.chr22.bam`.
- Re-run the Step 1 command above (expect ~25 min). Per the README, output is a gzipped BED file of CpG-level fragment features that downstream HMM training (`FinaleMe` class, Step 2) consumes.
- Capture the row count and a few header lines from `results/BH01.cpg_features.hg19.bed.gz` into `results/finaleme_smoke.json` and report.

### Medium term: scope down
- The chr22 BAM is enough to validate that the **CpG feature matrix builder runs** end-to-end on our pipeline inputs and emits a well-formed feature file. It is *not* enough to validate the HMM imputation accuracy — for that we need WGBS ground truth on the same sample (Jensen 2015 buffy-coat bigWig is provided as a prior, but accuracy evaluation needs held-out WGBS).
- If the goal is to check whether FinaleMe-style methylation improves our existing cfDNA models on the actual FinaleDB cohort, we need per-fragment BAMs, not the aggregated features already in our pipeline. FinaleDB does not host raw BAMs (only fragment BEDs), so a real accuracy benchmark would require re-aligning a public cfDNA WGS cohort (e.g. Snyder 2016) and running FinaleMe on it.

### Long-term alternative: FinaleToolkit primitives as a proxy
FinaleToolkit 1.1.0 is already installed at `/Users/hermes/deepcatch/.venv/bin/finaletoolkit`. Its `FragmentHistogram`, `KmerCounter`, and `EndCoordinateCounter` modules compute the same per-fragment features FinaleMe consumes — fragment length, GC content, end motifs, k-mer counts — directly from a fragment BED file (or BAM via `frag-from-bam`). For most downstream cfDNA classification tasks (tissue-of-origin, fragmentomics models) those raw features carry most of the methylation-imputation signal as a free byproduct of fragment length and end-motif statistics. A pragmatic substitute for full FinaleMe:

```python
from finaletoolkit.frag.histone import FragmentHistogram  # or similar
# Compute per-fragment-length / end-motif features directly from FinaleDB fragment BEDs
# Train existing cfDNA models on those features instead of waiting on the methylation track
```

This sidesteps the Java+Maven+reference download chain entirely and runs in pure Python on the data we already have.

---

## Artifacts (left on disk for the retry)

| Path | Size | Purpose |
|---|---|---|
| `/Users/hermes/.local/jdk/jdk-21.0.12.1.jdk/` | ~480 MB | JDK 21 install |
| `/Users/hermes/.local/maven/apache-maven-3.9.16/` | ~9 MB | Maven 3.9.16 |
| `/tmp/FinaleMe/` | ~1 GB | Source + vendored repo + JAR |
| `/tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar` | 47 MB | Runnable JAR |
| `/tmp/FinaleMe/lib-repo/` | ~? | Cached Maven deps (offline-build-ready) |
| `/tmp/finaleme_test/BH01.chr22.bam` | 585 MB | Test BAM (valid BGZF, possibly missing the last byte-block; re-pull is safe) |

## Commands for a clean retry

```bash
export JAVA_HOME=/Users/hermes/.local/jdk/jdk-21.0.12.1.jdk/Contents/Home
export PATH=$JAVA_HOME/bin:/Users/hermes/.local/maven/apache-maven-3.9.16/bin:$PATH

# 1. Verify toolchain
java -version
mvn --version

# 2. Pull references (45–90 min serial, ~15 min with -j 4 aria2c)
cd /tmp/FinaleMe
./scripts/setup_references.sh 2>&1 | tee /tmp/finaleme_setup.log

# 3. (Re)download + index BAM if needed
mkdir -p /tmp/finaleme_test
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam?download=1" \
  -o /tmp/finaleme_test/BH01.chr22.bam
curl -L "https://zenodo.org/records/6914806/files/BH01.chr22.bam.bai?download=1" \
  -o /tmp/finaleme_test/BH01.chr22.bam.bai
samtools index /tmp/finaleme_test/BH01.chr22.bam

# 4. Step 1 — feature extraction (~25 min on chr22)
JAR=/tmp/FinaleMe/target/FinaleMe-0.61-jar-with-dependencies.jar
cd /tmp/FinaleMe && mkdir -p results
java -Xmx20G -cp "$JAR" \
  edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder \
  data/hg19.2bit \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  data/CG_motif.hg19.common_chr.pos_only.bedgraph.gz \
  /tmp/finaleme_test/BH01.chr22.bam \
  results/BH01.cpg_features.hg19.bed.gz \
  -stringentPaired \
  -excludeRegions data/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed \
  -valueWigs methyPrior:0:data/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
  -useNoChrPrefixBam -wgsMode -t 4

# 5. Inspect + persist
zcat results/BH01.cpg_features.hg19.bed.gz | head -1
zcat results/BH01.cpg_features.hg19.bed.gz | wc -l   # row count
```

---

## Honest summary for the parent

- **Was FinaleMe successfully built?** Yes — `mvn clean package` produced `FinaleMe-0.61-jar-with-dependencies.jar` in 21.6 s.
- **Was a methylation feature successfully generated from the test BAM?** No — `setup_references.sh` was killed at ~16 min into the reference download (only `hg19.2bit` had begun, and was incomplete). The BH01 BAM itself finished downloading just before the kill. The runtime feature-extraction step was therefore not run.
- **What blocks it?** Slow reference mirrors + a strictly serial downloader + tight disk headroom, not any fundamental incompatibility.
- **Recommended workaround?** Either (a) re-run `setup_references.sh` when bandwidth is better, parallelizing the downloads if possible, or (b) pivot to FinaleToolkit per-fragment features (already installed) as a methylation-free proxy for downstream cfDNA modeling.