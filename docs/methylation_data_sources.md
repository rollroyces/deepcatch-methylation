# Publicly Accessible Methylation Data Sources for cfDNA Cancer Detection

Survey to identify public methylation data that could extend the DeepCatch + cfdna-fragmentomics-pipeline project (currently fragmentomics only — FSD and DELFI 5 Mb / 100 kb ratios from FinaleDB WGS data) with CpG-methylation features for cfDNA cancer detection.

Compiled September 2026 from web search of primary literature, GEO, GDC, dbGaP, EGA, and the cfMethDB / FinaleDB catalogues. No code was written.

---

## Executive Summary

| Resource | Type | Cancer types | Cancer / Healthy | Assay | Open access? | Best use for the project |
|---|---|---|---|---|---|---|
| **TCGA (GDC)** | Tumor tissue + matched-normal | 33 | thousands (tissue) | 450 K (mostly) + some 27 K; few EPIC | Open (dbGaP) | Pan-cancer DMR discovery, paired normal control for in-silico mixture |
| **cfMethDB** | Aggregated cfDNA methylation | 7 | 4 828 datasets | mixed (WGBS, RRBS, EM-seq, 450 K, EPIC) | Open, no login | Curated biomarker catalogue + pre-computed DMCs |
| **GSE122126 (Moss 2018)** | cfDNA + 25 cell-type reference | colorectal, lung, breast | 11 cancer / 59 cfDNA total | Illumina 450 K / EPIC | Open via GEO | Cell-of-origin deconvolution NNLS atlas |
| **GSE186458 (Loyfer 2023)** | Healthy tissue WGBS reference | n/a (healthy) | 205 samples, 39 cell types | deep WGBS (~30×) | EGA DAC (study EGAS00001006791) | Stronger normal reference for deconvolution |
| **FinaleMe / FinaleToolkit** | Method (predicts CpG from WGS) | breast, prostate, HCC, others | 77 paired ULP-WGS / WGBS | WGBS + WGS (ULP and deep) | Code open (MIT); training dbGaP phs003287 | Lets us **recover CpG methylation from existing WGS** — no new assay needed |
| **FinaleDB** | cfDNA fragmentation | 23 conditions | 2 579 WGS samples | WGS only (NO methylation) | Open web portal | Already used. No methylation here. |
| **dbGaP phs000846 (Sun 2015/2018)** | cfDNA WGBS | HCC (small) | ~37 healthy + 8–9 HCC | WGBS | dbGaP DAC | WGBS paired tissue / plasma benchmark |
| **dbGaP phs001417 (Adalsteinsson 2017)** | cfDNA WGS/WES | metastatic breast | n/a — controls | WGS | dbGaP DAC | Large healthy-plasma WGS, methylation predictible via FinaleMe |
| **dbGaP phs003287 (FinaleMe)** | Newly generated cfDNA | breast, prostate | restricted | WGBS + ULP-WGS + ULP-WGBS | dbGaP DAC | Future re-use benchmark |
| **GRAIL CCGA / PATHFINDER** | Multi-cancer cfDNA methylation | 50+ cancer types | thousands | Targeted bisulfite | NO — controlled, no public deposit | Aspirational reference only |

**Bottom line:** There is no open WGBS cfDNA cohort with the scale of CCGA. The realistic public path for the project is **(a)** train on TCGA tissue-level 450 K DMRs, **(b)** validate on GSE122126 cfDNA arrays and cfMethDB-curated markers, and **(c)** **infer CpG methylation on the existing FinaleDB WGS samples using FinaleMe** rather than waiting for new methylation data.

---

## 1. TCGA Methylation Data

- **URL:** https://portal.gdc.cancer.gov/  |  https://gdc.cancer.gov/about-data/gdc-data-processing/resources-tcga-users
- **Pipeline:** GDC Methylation Array Harmonization Workflow (raw IDAT → Level 2 beta values → Level 3); supports HM27, HM450 (450 K), and EPIC. (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_Pipeline/)
- **Cancer types covered:** 33 tumor types in TCGA. 450 K was the workhorse for the bulk of the program; 27 K appears in the earliest samples; EPIC only appears in a small subset of projects / contributed data.
- **Sample counts (illustrative;** `TCGAbiolinks::GDCquery()` is the right way to get exact numbers for any single project):
  - TCGA-BRCA (breast): ~800 tumor / ~100 matched-normal on 450 K
  - TCGA-LUAD + LUSC (lung): ~900 tumor / ~70 matched-normal
  - TCGA-LIHC (liver / HCC): ~380 / 50
  - TCGA-COAD + READ (colorectal): ~520 / 80
  - TCGA-KIRC + KIRP + KICH (kidney): ~900 / ~250
  - TCGA-PAAD (pancreas): ~185 / 10
  - TCGA-OV (ovarian): ~600 / 0 matched
  - TCGA-STAD (gastric): ~450 / 2
  - TCGA-ESCA: ~200 / 16
  - TCGA-HNSC: ~580 / 50
  - TCGA-BLCA: ~450 / 20
  - TCGA-PRAD (prostate): ~550 / 50
  - TCGA-THCA: ~570 / 56
  - TCGA-SKCM: ~470 / 2
  (Numbers vary by query date; verify via `GDCquery()` or `TCGAbiolinks`.)
- **Matched tumor + normal?** YES — TCGA was specifically designed for matched pairs (blood-derived normal or adjacent-tissue normal). Good for tumor-vs-normal DMR calling.
- **Matched cfDNA + tumor?** NO. TCGA is tissue (mostly fresh-frozen primary tumor) plus germline normal. No plasma from the same patients.
- **Access requirements:** Open (no DAC) for tumor and normal 450 K IDATs / Level 3 beta values, after a free dbGaP / eRA Commons login for controlled-access tiers.
- **Estimated download size:** A single project's full Level 3 450 K matrix ≈ 5–25 GB compressed; a 30-project pull is ≈ 200 GB.
- **Companion resource:** PanCanAtlas publication archive at https://gdc.cancer.gov/about-data/publications/pancanatlas provides pre-processed merged 27 K + 450 K matrices.
- **How it fits the project:** Use TCGA as the gold-standard source of cancer-vs-normal DMRs. Probe sets are array-based, so they need to be lifted to WGBS coordinates (or you must restrict FinaleDB analysis to CpGs that overlap 450 K probes). For DELFI / FSD extension, the most direct use is to score each FinaleDB fragment by whether its CpGs fall inside a TCGA cancer-DMR.

---

## 2. GEO cfDNA Methylation Datasets

GEO is the most accessible deposit. Below are the accessions most relevant to cfDNA cancer detection.

### GSE122126 — Moss et al. 2018, Nature Communications (the original 25-tissue atlas)
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122126
- **Paper:** https://www.nature.com/articles/s41467-018-07466-6
- **Assay:** Illumina HumanMethylation 450 K + EPIC (the deconvolution reference is 450 K)
- **Samples:** 59 plasma cfDNA + 27 sorted cell populations + 15 technical controls. Includes 4 colorectal, 4 lung, 3 breast cancers and 4 CUP patients.
- **Cancer vs healthy:** 11 cancer, 4 healthy cfDNA directly; the rest is reference atlas.
- **Match type:** cfDNA only — no matched tumor tissue.
- **Access:** Open via `GEOquery` / `getGEO()`. No DAC.
- **Size:** ≈ 1 GB IDAT.
- **Use:** Standard cell-of-origin deconvolution benchmark. Already integrated in `nloyfer/meth_atlas` (Python) and `meth_atlas` R package. **Tier-3 of the existing pipeline.**

### GSE186458 — Loyfer et al. 2023, Nature (Loyfer Atlas)
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458   (Processed matrices only)
- **Raw FASTQ / BAM:** EGA study **EGAS00001006791** (DAC-controlled)
- **Paper:** https://www.nature.com/articles/s41586-022-05580-6
- **Assay:** Whole-genome bisulfite sequencing (WGBS), 30× paired-end NovaSeq 6000.
- **Samples:** 205 healthy-tissue samples, 39 sorted cell types.
- **Cancer vs healthy:** Healthy only — no cancer plasma.
- **Match type:** Tissue reference methylomes (no plasma in this study, but a smaller companion plasma panel exists at EGA).
- **Access:** Processed matrices open via GEO; raw FASTQs require an EGA Data Access Agreement (DAA template at https://www.cs.huji.ac.il/~tommy/DAA/DAA.Human_Methylation_Atlas.docx). Typically approved within weeks.
- **Size:** Raw ≈ 5–10 TB (205 samples × ~30× WGBS).
- **Use:** Stronger reference for NNLS-style deconvolution; supersedes Moss 2018 in cell-type resolution.

### GSE97923 / GSE97932 — colorectal / lung cancer cfDNA methylation
- **URLs:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97923 ; GSE97932
- **Paper:** https://www.nature.com/articles/s12276-023-01119-5 (Han et al., 2023)
- **Assay:** HM450 K arrays on plasma cfDNA + matched tumor + matched normal blood.
- **Samples:** ≈ 70 COAD pairs; similar for LUSC.
- **Access:** Open via GEO.

### GSE63775 / GSE70090 / GSE112221 / GSE55752 — HCC (Wen et al.)
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63775
- **Assay:** Targeted bisulfite sequencing (methylated CpG tandem amplification & sequencing); some WGBS of cfDNA + NAT (normal adjacent tissue) + plasma.
- **Samples:** 191 HCCs.
- **Access:** Open.

### GSE186458 / dbGaP phs000846 — Sun et al. plasma tissue mapping
- **URL:** https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000846.v1.p1
- **Assay:** Whole-genome bisulfite sequencing of cfDNA + gDNA buffy coat.
- **Samples:** 37 healthy + 8–9 HCC (pre-/post-surgery).
- **Access:** dbGaP DAC (typically approved for academic users in days–weeks).
- **Use:** The dataset FinaleMe was benchmarked against.

### dbGaP phs001417 — Adalsteinsson et al. 2017
- **URL:** https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001417.v1.p1
- **Assay:** cfDNA WGS / WES; ~16–39× coverage. Plus the FinaleMe cohort has 77 matched ULP-WGS + ULP-WGBS on breast, prostate, healthy.
- **Access:** dbGaP DAC.

### dbGaP phs003287 — FinaleMe generated WGBS / WGS
- **URL:** https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003287.v1.p1
- **Use:** Future re-use — required for cancer-specific methylation-on-WGS training.

### Smaller / less-aligned GEO sets (quick reconnaissance)
- **GSE248620** — UHRF1 MCF-7 WGBS (cell-line, low priority).
- **GSE124686, GSE81314** — Snyder / Ulz cfDNA WGS (no methylation; potential FinaleMe target).
- **GSE71378** — Snyder 2016 (BH01/IH01 healthy plasma WGS; no methylation but FinaleMe-ready).

### SRA accession patterns for cfDNA methylation sequencing
Most cfDNA methylation WGBS / RRBS / EM-seq deposits in SRA share the prefix `phs` (dbGaP) or simply appear as paired SRA runs (`SRR/ERR/DRR` numbers) under a `PRJNA` / `PRJEB` / `PRJCA` umbrella. Looking for cfDNA methylation data:

```
srapath `query --type study \
  --contains "cell-free DNA methylation" \
  --filter assay_type:"Bisulfite-Seq"`
```

or via the web: https://www.ebi.ac.uk/ena/browser/text-search?query=cell-free%20DNA%20methylation%20bisulfite

---

## 3. FinaleDB — Methylation Coverage?

**No.** FinaleDB hosts **fragmentation only**, never CpG-level methylation.

- **URL:** http://finaledb.research.cchmc.org/
- **Paper:** Zheng, Zhu & Liu, *Bioinformatics* 37(16):2502–2503 (2021) — DOI 10.1093/bioinformatics/btaa999.
- **Content:** 2 579 paired-end cfDNA WGS datasets across 23 pathological conditions (GEO / EGA / dbGaP-sourced); uniformly reprocessed with `finaledb_workflow`. Delivers fragment BED / bigWig / coverage / size profile tracks.
- **What it does NOT do:** methylation calling, CpG-level features, fragment-end motif methylation. The 'methylation' connection is indirect — `FinaleMe` (below) was built by the same CCHMC group to *predict* methylation from FinaleDB-style WGS.
- **Source code:** https://github.com/epifluidlab/finaledb_portal ; https://github.com/epifluidlab/finaledb_workflow

**For methylation features, FinaleDB data must be either (a) co-analyzed with externally downloaded methylation measurements, or (b) processed by `FinaleMe` to impute CpG methylation from the WGS fragments already available.** This is the cheapest path forward because no new raw data download is needed.

---

## 4. GRAIL CCGA / PATHFINDER — Public Availability?

**Public methylation data: NO.**

- Liu et al., *Annals of Oncology* 31:745–759 (2020) (CCGA substudy 1) describes a **targeted bisulfite sequencing panel (~17 000 CpGs)** across ~3 500 participants (1 600 non-cancer + 1 422 cancer covering >50 types). The training set draws on WGBS of tissue + plasma from CCGA plus commercially sourced tissue.
- CCGA / PATHFINDER / Galleri / G360 raw data and feature matrices are **proprietary** and not deposited in GEO / SRA / EGA / TCGA.
- Paper supplements contain classifier AUC figures and small marker lists, **not the per-sample methylation matrices**.
- The only public tie-in is the **CCGA2 preprint / DDW2020 posters** at https://grail.com/wp-content/uploads/2020/12/DDW_2020_Liu_CCGA2_Encore_POS_Final.pdf — no data.
- Final CCGA-3 multi-modality results are in *Cancer Cell* (2022) — press release at https://grail.com/press-releases/grail-ccga-discovery-results-published-in-cancer-cell-reveal-methylation-as-promising-dna-hallmark-for-multi-cancer-early-detection/ — again no data.

**What this means for the project:** Treat CCGA performance numbers as the **target benchmark**, not as a training resource. Build with TCGA + public cfDNA + FinaleMe-imputed WGS.

---

## 5. Other Public cfDNA Methylation Resources

### cfMethDB — comprehensive meta-database (Huazhong Agricultural U.)
- **URL:** https://cfmethdb.hzau.edu.cn/home
- **Paper:** Sun et al., *Genomics, Proteomics & Bioinformatics* (2025), DOI 10.1093/gpbjnl/qzaf092.
- **Content:** 4 828 curated cfDNA methylation datasets across 7 cancer types (hepatocellular, breast, colorectal, lung, gastric, pancreatic, ovarian). 54 cancer subtypes; 729 marker genes; 1 048 770 differentially methylated cytosines.
- **Caveat:** It is a *catalogue* — most entries point back to SRA / EGA / GEO for the raw data; only the curated DMCs and biomarker evaluations are first-party.
- **Access:** Open web query (gene, region, motif, biomarker evaluation). No login.
- **Use:** Quick lookup of DMRs / gene-level marker panels. Skip re-discovering biomarkers.

### DiseaseMeth / MethBank / MethAtlas — human disease methylation catalogues
- **methbank:** https://ngdc.cncb.ac.cn/methbank  (also hosts GSE186458 — Loyfer Atlas)
- **disease-meth:** http://bio-bigdata.hrbmu.edu.cn/disease_meth/
- **methatlas:** https://sun-lab.cn/methatlas (recent, multi-omic)
- **Use:** Lower-priority. Cross-validate DMRs from TCGA / cfMethDB.

### Shen et al. 2018 — Plasma cfDNA methylome immunoprecipitation
- **Paper:** https://www.nature.com/articles/s41586-018-0703-0 (Shen, Singhania et al., *Nature* 563:579)
- **Assay:** cfDNA methylome immunoprecipitation + sequencing (cfMeDIP-seq).
- **Samples:** Healthy + multiple tumor types (breast, colorectal, lung, ovarian, etc.).
- **Data availability:** Raw data partly in EGA under restricted access; protocol widely re-used by the field. Follow-up cfMeDIP-seq studies (BC, HCC) appear in GEO (e.g., GSE164138 is a related 5-formyluracil method, not this protocol).

### BASIS — breast cancer WGBS classifier
- **Paper:** https://www.nature.com/articles/s43018-022-00415-9
- **Sample count:** WGBS of breast PDX / cfDNA pairs (tens of samples).
- **Use:** Subtype-specific deconvolution reference.

### Misc targeted-bisulfite / RRBS sets (small but useful)
- **GSE100824** (leukocyte DNA methylation atlas, Hannum)
- **HCC epiLiver** — Wen et al., commercially marketed methylation assay, raw at GSE63775.
- **HCC cfMeDIP-seq** — recent deposit referenced in the 2025 ACS Brief (https://www.facs.org/for-medical-professionals/news-publications/news-and-articles/acs-brief/february-10-2026-issue/plasma-cell-free-dna-methylation-analysis-is-highly-sensitive-to-detecting-hepatocellular-carcinoma/) — verify GEO accession when integrating.

---

## FinaleMe — The Practical Path Forward

The single most important practical finding for the project:

> **FinaleMe (Liu et al., *Nature Communications* 15:2790, 2024) predicts single-CpG methylation status and tissue-of-origin from plain cfDNA WGS — no bisulfite conversion, no extra assay.**

- **Paper:** https://www.nature.com/articles/s41467-024-47196-6
- **Code:** https://github.com/epifluidlab/FinaleMe (MIT) ; DOI 10.5281/zenodo.7779198
- **Validation:** 80 paired deep-WGS / WGBS + 77 paired ULP-WGS / ULP-WGBS (breast, prostate, healthy) — auROC = 0.91 for fragments with ≥5 CpGs in CpG-rich regions.
- **Trained on:** Sun 2015 / Sun 2018 WGBS (dbGaP phs000846, phs001417); newly generated WGBS in phs003287.
- **Limitations:** Performance degrades in CpG-poor regions; relies on a Bayesian prior from healthy gDNA; some DMR false positives in cancer WGS.
- **Compatibility with FinaleDB:** Designed for the same fragment feature extraction that FinaleToolkit (https://github.com/epifluidlab/FinaleToolkit ; https://pypi.org/project/FinaleToolkit/) outputs. You can run FinaleMe on the same BED/bigWig fragments the project already ingests.

**Implication for the pipeline:** Existing WGS samples in FinaleDB can be re-analyzed with FinaleMe to get imputed CpG methylation. This is the most pragmatic way to add a methylation feature to the existing DELFI / FSD pipeline without a single new raw-data download.

---

## Recommended Data Strategy for the Project

1. **Tissue reference DMRs (training labels):** TCGA 450 K via `TCGAbiolinks::GDCquery()` for pan-cancer DMR discovery. Restrict to top cancer types of interest (BRCA, LUAD, LIHC, COAD, PAAD).
2. **Normal-cell reference:** Loyfer atlas (GSE186458) processed matrices (open) or Moss atlas (GSE122126).
3. **cfDNA array validation:** GSE122126, GSE97923, GSE97932.
4. **In-silico mixture ground truth:** TCGA tumor beta values mixed in known proportions with healthy buffy-coat reference — produces thousands of synthetic "cfDNA" profiles without sequencing a single plasma sample. (Pattern already demonstrated in KelyNorel/cfdna-methylation-atlas.)
5. **WGS → methylation imputation on FinaleDB samples:** Run FinaleMe on Snyder 2016 (GSE71378), Cristiano 2019, Jiang 2018, Sun 2019 — i.e., the existing WGS the project already pulls. This is the cheapest new feature source.
6. **Benchmark target:** CCGA AUCs (proprietary) as the upper bound to track against.

**Estimated total download budget for steps 1–5:** ≈ 1–3 TB raw WGBS data + ≈ 200 GB TCGA arrays. Workable on a single workstation with 2 TB SSD; the FinaleDB-side work is essentially zero new bytes.

---

## Appendix: Quick Links

- GDC: https://portal.gdc.cancer.gov/
- TCGA Wanderer (pre-aggregated): http://maplab.imppc.org/wanderer/
- cfMethDB: https://cfmethdb.hzau.edu.cn/home
- FinaleDB: http://finaledb.research.cchmc.org/
- FinaleMe code: https://github.com/epifluidlab/FinaleMe
- FinaleToolkit: https://github.com/epifluidlab/FinaleToolkit
- meth_atlas (Moss): https://github.com/nloyfer/meth_atlas
- Loyfer atlas processed (GEO): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458
- Loyfer atlas raw (EGA): https://ega-archive.org/studies/EGAS00001006791
- Weiss et al. CCGA-3 (CCGA Cancer Cell 2022): https://grail.com/press-releases/grail-ccga-discovery-results-published-in-cancer-cell-reveal-methylation-as-promising-dna-hallmark-for-multi-cancer-early-detection/

