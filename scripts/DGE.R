#!/usr/bin/env Rscript
# Differential expression (DESeq2) with CLI parameters and robust defaults.
# Outputs are modularised under <outdir>/dge/. No hard-coded absolute paths.
#
# CLI inputs:
#   --counts     TSV/CSV (genes × samples). First column = gene id; other columns = sample IDs.
#   --meta       TSV/CSV with at least: sample_id, group.
#   --outdir     Output root (default: results) → module outputs to <outdir>/dge/
#   --ref        (optional) reference level of 'group' (e.g., "HC")
#   --alt        (optional) comparison level of 'group' (e.g., "Resister")
#   --alpha      FDR (BH) significance threshold [default: 0.05]
#   --zero_frac  Filter genes with proportion of zeros > zero_frac (relaxed for tiny toy data) [default: 0.95]
#
# Outputs (under <outdir>/dge/):
#   tables/: deseq2_results.tsv; deseq2_ranked_genes.tsv; normalized_counts.tsv
#   figures/: volcano.png; ma_plot.png
#
# Recommended citations (methods & software):
#   • Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data
#     with DESeq2. Genome Biology 15:550. doi:10.1186/s13059-014-0550-8
#   • Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach
#     to multiple testing. JRSS Series B 57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
#   • Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions and the apeglm method for RNA-seq.
#     Bioinformatics 35:385–393. doi:10.1093/bioinformatics/bty895
#   • Huber W, Carey VJ, Gentleman R, et al. 2015. Orchestrating high-throughput genomic analysis with
#     Bioconductor. Nature Methods 12:115–121. doi:10.1038/nmeth.3252
#   • Wickham H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer. doi:10.1007/978-3-319-24277-4
#
# Reproducibility: print sessionInfo() at the end to capture exact R/package versions.
# Author: Junyuan Guan (Imperial College London)

suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(readr); library(dplyr); library(tibble); library(ggplot2)
})

# CLI
option_list <- list(
  make_option(c("-c","--counts"),   type="character"),
  make_option(c("-m","--meta"),     type="character"),
  make_option(c("-o","--outdir"),   type="character", default="results"),
  make_option(     "--ref",         type="character", default=NULL, help="reference level of group"),
  make_option(     "--alt",         type="character", default=NULL, help="comparison level of group"),
  make_option(     "--alpha",       type="double",    default=0.05),
  make_option(     "--zero_frac",   type="double",    default=0.95)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$meta)) stop("Please provide --counts and --meta.", call.=FALSE)

# I/O dirs
outdir_dge <- file.path(opt$outdir, "dge")
figdir      <- file.path(outdir_dge, "figures")
tabdir      <- file.path(outdir_dge, "tables")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)

# read
counts <- readr::read_tsv(opt$counts, show_col_types = FALSE)
meta   <- readr::read_tsv(opt$meta,   show_col_types = FALSE)

# enforce first column is 'gene'
if (!("gene" %in% names(counts))) names(counts)[1] <- "gene"
stopifnot("gene" %in% names(counts))
stopifnot(all(c("sample_id","group") %in% names(meta)))

# align columns to metadata sample order
sample_ids <- meta$sample_id
stopifnot(all(sample_ids %in% colnames(counts)))
counts <- counts[, c("gene", sample_ids)]

# matrix
cts <- as.matrix(counts[,-1]); mode(cts) <- "integer"; rownames(cts) <- counts$gene

# normalise group labels to canonical levels: HC / Converter / Resister
meta$group <- as.character(meta$group)
meta$group <- tolower(meta$group)
meta$group <- dplyr::case_when(
  grepl("resist",   meta$group) ~ "Resister",
  grepl("convert",  meta$group) ~ "Converter",
  grepl("control|hc|healthy", meta$group) ~ "HC",
  TRUE ~ stringr::str_to_title(meta$group)
)
meta$group <- factor(meta$group)

# (optional) print levels to help debugging
message("group levels: ", paste(levels(meta$group), collapse=", "))

# optional ref/alt handling
if (!is.null(opt$ref) && !(opt$ref %in% levels(meta$group))) {
  stop(sprintf("--ref '%s' not found. Available levels: %s",
               opt$ref, paste(levels(meta$group), collapse=", ")), call.=FALSE)
}
if (!is.null(opt$alt) && !(opt$alt %in% levels(meta$group))) {
  stop(sprintf("--alt '%s' not found. Available levels: %s",
               opt$alt, paste(levels(meta$group), collapse=", ")), call.=FALSE)
}

# prefilter (relaxed for tiny data)
keep <- rowSums(cts==0)/ncol(cts) <= opt$zero_frac
if (sum(keep) >= 50) cts <- cts[keep, , drop=FALSE]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~ group)
dds <- DESeq(dds, fitType = "mean")  # stable on tiny toy sets

# normalized counts (size-factor) for reference/export
norm_counts <- counts(dds, normalized = TRUE)
readr::write_tsv(
  as.data.frame(norm_counts) |> tibble::rownames_to_column("gene"),
  file.path(tabdir, "normalized_counts.tsv")
)

# results: contrast if ref+alt provided; else default coefficient
get_results <- function(dds, alpha, ref=NULL, alt=NULL){
  if (!is.null(ref) && !is.null(alt)) {
    # alt vs ref (i.e., log2FC = alt - ref)
    results(dds, contrast = c("group", alt, ref), alpha = alpha)
  } else {
    results(dds, alpha = alpha)  # uses current reference; fine for 2-level factors
  }
}
res <- get_results(dds, alpha = opt$alpha, ref = opt$ref, alt = opt$alt)
res_df <- as.data.frame(res) |> tibble::rownames_to_column("gene")

# optional LFC shrinkage (apeglm if available and coefficients compatible)
shrunken <- NULL
suppressWarnings({
  if (requireNamespace("apeglm", quietly = TRUE)) {
    # try coef=2 (first non-intercept) unless contrast used
    if (is.null(opt$ref) || is.null(opt$alt)) {
      shrunken <- try(lfcShrink(dds, coef=2, type="apeglm"), silent=TRUE)
    } else {
      # when explicit contrast is used, specify via 'contrast' form
      shrunken <- try(lfcShrink(dds, contrast=c("group", opt$alt, opt$ref), type="apeglm"), silent=TRUE)
    }
  }
})
if (!inherits(shrunken, "try-error") && !is.null(shrunken)) {
  res_df <- as.data.frame(shrunken) |> tibble::rownames_to_column("gene")
}

# ranked list for GSEA
ranked <- res_df |> dplyr::filter(!is.na(stat)) |> dplyr::arrange(dplyr::desc(stat)) |> dplyr::select(gene, stat)

# write tables
readr::write_tsv(res_df, file.path(tabdir, "deseq2_results.tsv"))
readr::write_tsv(ranked, file.path(tabdir, "deseq2_ranked_genes.tsv"))

# plots
safe_neglog10 <- function(x){ y <- -log10(x); y[!is.finite(y)] <- NA_real_; y }
alpha <- opt$alpha

# Volcano
if (all(c("log2FoldChange","padj") %in% names(res_df))) {
  vdf <- res_df |> dplyr::filter(is.finite(log2FoldChange))
  vdf$neglog10padj <- safe_neglog10(vdf$padj)
  vdf$signif <- !is.na(vdf$padj) & (vdf$padj <= alpha)
  p_volc <- ggplot(vdf, aes(x=log2FoldChange, y=neglog10padj, color=signif)) +
    geom_point(alpha=0.7, size=1.6) +
    theme_minimal(base_size=12) +
    scale_color_manual(values=c("FALSE"="grey60","TRUE"="tomato")) +
    labs(title=sprintf("Volcano (FDR ≤ %.2f)", alpha),
         x="log2 fold change", y="-log10(padj)", color="Significant")
  ggsave(file.path(figdir, "volcano.png"), p_volc, width=7.5, height=5.5, dpi=300)
}

# MA plot
if (all(c("baseMean","log2FoldChange") %in% names(res_df))) {
  mdf <- res_df |> dplyr::mutate(A = log10(baseMean + 1))
  p_ma <- ggplot(mdf, aes(x=A, y=log2FoldChange)) +
    geom_point(alpha=0.6, size=1.4, color="grey50") +
    geom_hline(yintercept=0, linetype="dashed") +
    theme_minimal(base_size=12) +
    labs(title="MA plot", x="log10(mean normalised count + 1)", y="log2 fold change")
  ggsave(file.path(figdir, "ma_plot.png"), p_ma, width=7.5, height=5.5, dpi=300)
}

message("DGE complete. Outputs in: ", normalizePath(outdir_dge))
message("\n--- sessionInfo() ---")
print(utils::sessionInfo())