# TB RNA-seq Analysis (Imperial NHLI): Summary & Reproducible Demo

**Role/Period:** Researcher, May–Aug 2024  
**Affiliation:** National Heart & Lung Institute, Imperial College London

This repository provides a **reproducible, CLI-driven downstream RNA-seq pipeline** starting from **gene-level count matrices** and **sample metadata**. It includes:
- **Differential expression (DESeq2)** with Benjamini–Hochberg FDR control and optional LFC shrinkage;
- **VST → PCA** with variance explained, loadings, and exposure-metric correlations;
- **GSEA dot-plots** from summary tables;
- **SCBD (CTV) visualisation** with heatmaps and nonparametric Wilcoxon tests.

> **Scope note.** Read alignment/quantification (e.g., FASTQ → STAR/Salmon → counts) is **out of scope**. Supply a gene-by-sample count matrix from your preferred upstream workflow.

**Input assumptions**
- Counts are **raw integer** gene counts (not TPM/FPKM/CPM) for DESeq2.
- Metadata contains `sample_id` and `group`. Common labels are **normalised** to `HC / Converter / Resister` (e.g., “Healthy control” → `HC`).
- Gene identifiers are consistent across counts and downstream summaries (GSEA/SCBD).

---

## Study snapshot (context)

- **Cohort:** High-exposure TB contacts, **n = 31** (IGRA-defined: 13 resisters, 8 converters, 10 controls).  
- **Objective:** Identify transcriptomic and pathway-level signals that distinguish **resisters** and relate to **exposure metrics**.  
- **Highlights:** Outlier-aware PCA; DESeq2 with BH–FDR (adj. *p* ≤ 0.05); GSEA on MSigDB KEGG/Hallmark v2023.2 (FDR *q* ≤ 0.25).

---

## Quick start

```bash
# 1) Clone
git clone https://github.com/JUNYUAN-GUAN/TB-rnaseq-pipeline.git
cd TB-rnaseq-pipeline

# 2) Install R packages (R ≥ 4.3 recommended)
R -q -e 'install.packages(c("optparse","readr","dplyr","tibble","ggplot2"), repos="https://cloud.r-project.org");
         if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager");
         BiocManager::install(c("DESeq2","pheatmap","apeglm"), ask=FALSE)'

# 3) One-click demo on toy data (writes to results/{dge,pca,gsea,ctv})
chmod +x scripts/run_all.sh
scripts/run_all.sh \
  -c example_data/counts.tsv \
  -m example_data/metadata.tsv \
  -o results \
  --ref HC --alt Resister