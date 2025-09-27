# TB RNA-seq Analysis (Imperial NHLI): Summary & Reproducible Demo

**Role/Period:** Researcher, May–Aug 2024  
**Affiliation:** National Heart & Lung Institute, Imperial College London

A **reproducible, CLI-driven downstream RNA-seq pipeline** starting from **gene-level count matrices** and **sample metadata**. It includes:
- **Differential expression (DESeq2)** with Benjamini–Hochberg FDR and optional LFC shrinkage;
- **VST → PCA** with variance explained, loadings, and exposure-metric correlations;
- **GSEA dot-plots** from summary tables;
- **SCBD (CTV) visualisation** with heatmaps and nonparametric Wilcoxon tests.

> **Scope.** Read alignment/quantification (FASTQ → STAR/Salmon → counts) is **out of scope**. Provide a gene-by-sample count matrix from your upstream workflow.

---

## Quick start

```bash
# 1) Clone
git clone https://github.com/JUNYUAN-GUAN/TB-rnaseq-pipeline.git
cd TB-rnaseq-pipeline

# 2) Install deps (R ≥ 4.3)
R -q -e 'install.packages(c("optparse","readr","dplyr","tibble","ggplot2"), repos="https://cloud.r-project.org");
         if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org");
         BiocManager::install(c("DESeq2","pheatmap","apeglm"), ask=FALSE, update=FALSE)'

# 3) One-click demo (writes to results/{dge,pca,gsea,ctv})
chmod +x scripts/run_all.sh
scripts/run_all.sh \
  -c example_data/counts.tsv \
  -m example_data/metadata.tsv \
  -o results \
  --ref HC --alt Resister
```

> Each module prints `sessionInfo()` for reproducibility.  
> Exclude an outlier in PCA/CTV: `--exclude TBRHC-015`.  
> Run a subset: `--only dge,pca` (see below).

---

## Input assumptions

- **Counts**: raw integer gene counts (not TPM/FPKM/CPM), genes × samples; first column header **`gene`**; remaining columns are **sample IDs**.  
- **Metadata**: columns **`sample_id`**, **`group`**. Common labels are normalised to **HC / Converter / Resister** (e.g., “Healthy control” → `HC`).  
- **GSEA tables**: columns (case-insensitive) **Genes** (set name), **NES**, **FDR_q_val**/`qval`/`padj`, **Size**.  
- **SCBD tables**: columns **Feature**, **SCBD** (synonyms tolerated). Heatmaps use **normalised counts**.

**Examples**
```
# counts.tsv (TSV/CSV)
gene    s1  s2  s3  s4  s5  s6
G1      10  12  11  9   5   7
G2      0   0   1   2   3   0
GAPDH   100 110 95  120 105 98
```
```
# metadata.tsv (CSV/TSV)
sample_id,group
s1,HC
s2,HC
s3,Converter
s4,Converter
s5,Resister
s6,Resister
```

---

## Modules & usage

> Each script supports `--help` for full options.

### 1) DGE (DESeq2)

```bash
Rscript scripts/DGE.R \
  --counts path/to/counts.tsv \
  --meta   path/to/metadata.tsv \
  --outdir results \
  --ref HC --alt Resister \
  --alpha 0.05
```

- Zero filtering via `--zero_frac` (default 0.95; relaxed automatically for tiny matrices).  
- Uses `fitType="mean"` for stability; **apeglm** shrinkage if available.

**Outputs → `results/dge/`**
- `tables/deseq2_results.tsv`  
- `tables/deseq2_ranked_genes.tsv`  
- `tables/normalized_counts.tsv`  
- `figures/volcano.png`, `figures/ma_plot.png`

---

### 2) PCA (VST → PCA)

```bash
Rscript scripts/PCA.R \
  --counts path/to/counts.tsv \
  --meta   path/to/metadata.tsv \
  --outdir results \
  --exclude TBRHC-015 \
  --pcs 1,2
```

- VST via DESeq2; robust fallback to `log1p(normalised)` for very small matrices.  
- Computes optional Spearman correlations vs exposure metrics if present.

**Outputs → `results/pca/`**
- `pca_scores.tsv`, `pca_loadings.tsv`, `variance_explained.tsv`, `spearman_correlations.tsv`  
- `figures/pca_PC1_PC2.png`, `figures/pca_scree.png`, `figures/pc1_loadings_top30.png`

---

### 3) GSEA (dot-plots)

```bash
Rscript scripts/GSEA.R \
  --tables results/gsea/gsea_merged6.tsv,results/gsea/gsea_legacy.tsv \
  --labels Merged_6,Legacy \
  --outdir results \
  --basename gsea_all
```

**Outputs → `results/gsea/`**
- Per input: `<label>_table.tsv`, `figures/<label>_dotplot.png`  
- Combined (≥2 inputs): `gsea_all_table.tsv`, `figures/gsea_all_dotplot.png`

---

### 4) CTV (SCBD visualisation)

```bash
Rscript scripts/CTV.R \
  --scbd_raw      path/to/out_SCBD_raw_values.tsv \
  --scbd_filtered path/to/out_SCBD_filtered_values.tsv \
  --counts_norm   results/dge/tables/normalized_counts.tsv \
  --meta          path/to/metadata.tsv \
  --outdir        results \
  --exclude       TBRHC-015 \
  --scbd_threshold 0.0002 \
  --top_k 50 \
  --top_heatmap 15 \
  --remove_globin TRUE
```

**Outputs → `results/ctv/`**
- `tables/` SCBD top lists, heatmap gene list, per-gene medians, expression sums, Wilcoxon tests  
- `figures/` SCBD barplots（full/top-K）, `heatmap_topN.png`, per-gene boxplots, expression-sum boxplot

---

## One-click runner

```bash
chmod +x scripts/run_all.sh
scripts/run_all.sh \
  -c path/to/counts.tsv \
  -m path/to/metadata.tsv \
  -o results \
  --ref HC --alt Resister \
  --exclude TBRHC-015 \
  --only dge,pca
# Use --skip gsea,ctv to skip modules; see scripts/run_all.sh --help
```

A combined log is written to `results/run.log`.

---

## Expected outputs

```
results/
  run.log
  dge/
    tables/  deseq2_results.tsv
             deseq2_ranked_genes.tsv
             normalized_counts.tsv
    figures/ volcano.png
             ma_plot.png
  pca/
    pca_scores.tsv
    pca_loadings.tsv
    variance_explained.tsv
    spearman_correlations.tsv
    figures/ pca_PC1_PC2.png
             pca_scree.png
             pc1_loadings_top30.png
  gsea/
    <label>_table.tsv
    gsea_all_table.tsv
    figures/ <label>_dotplot.png
             gsea_all_dotplot.png
  ctv/
    tables/  scbd_raw_top50.tsv
             scbd_filtered_top50.tsv
             heatmap_topN_gene_list.tsv
             expression_medians.tsv
             expression_sum_by_sample.tsv
             wilcoxon_pairwise.tsv
    figures/ scbd_raw_barplot.png
             scbd_raw_top50_barplot.png
             scbd_filtered_barplot.png
             scbd_filtered_top50_barplot.png
             heatmap_topN.png
             boxplots_per_gene.png
             expression_sum_boxplot.png
```

---

## Repository structure

```
scripts/
  DGE.R  PCA.R  GSEA.R  CTV.R  run_all.sh
example_data/
  counts.tsv  metadata.tsv
figures/               # (optional) small demo images for README
results/               # created on run; ignored by git
README.md  LICENSE  .gitignore
```

---

## Reproducibility

- Scripts print `sessionInfo()` (R & package versions).  
- Recommend **R ≥ 4.3**, current Bioconductor; consider `renv` to pin versions.  
- No stochastic steps by default; if you add sampling-based procedures, set seeds explicitly.

---

## Citations

- Love MI, Huber W, Anders S. 2014. DESeq2. *Genome Biol* 15:550. doi:10.1186/s13059-014-0550-8  
- Benjamini Y, Hochberg Y. 1995. Controlling the FDR. *J R Stat Soc B* 57:289–300.  
- Zhu A, Ibrahim JG, Love MI. 2019. apeglm. *Bioinformatics* 35:385–393. doi:10.1093/bioinformatics/bty895  
- Jolliffe IT, Cadima J. 2016. PCA review. *Philos Trans R Soc A* 374:20150202.  
- Subramanian A, et al. 2005. GSEA. *PNAS* 102:15545–15550.  
- Liberzon A, et al. 2015. MSigDB Hallmark. *Cell Syst* 1:417–425.  
- Kanehisa M, Goto S. 2000. KEGG. *Nucleic Acids Res* 28:27–30.  
- Kolde R. 2019. **pheatmap** v1.0.12. CRAN.  
- Wickham H. 2016. **ggplot2**. Springer.  
- Wickham H. et al. 2019. **tidyverse**. *JOSS* 4(43):1686.  
- Huber W, et al. 2015. Bioconductor. *Nat Methods* 12:115–121.  
- R Core Team. *R: A language and environment for statistical computing*. R Foundation, Vienna. https://www.R-project.org/

---

## License

Code is released under the **MIT License** (see `LICENSE`).  
Figures/example datasets: verify reuse terms before redistribution.

---

## Contact

Maintainer: **Junyuan Guan** (Imperial College London)  
Issues and feature requests: please open a GitHub Issue.