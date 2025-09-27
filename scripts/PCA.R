#!/usr/bin/env Rscript
# PCA for RNA-seq count data with DESeq2 variance-stabilizing transform (VST) and CLI-friendly I/O.
# No hard-coded absolute paths; all inputs/outputs are parameterized.
#
# Inputs:
#   --counts     TSV/CSV count matrix (genes × samples). First column = gene id; remaining columns = sample IDs.
#   --meta       TSV/CSV metadata with at least: sample_id, group. (Optional exposure metrics allowed.)
#   --outdir     Output directory (default: results). Figures go to <outdir>/pca/figures.
#   --exclude    Comma-separated sample IDs to exclude (optional), e.g., "TBRHC-015".
#   --zero_frac  Gene is filtered if proportion of zeros > zero_frac (default: 0.95; relaxed for tiny toy data).
#   --pcs        Principal components to plot for the main scatter (default: "1,2"; e.g., "3,4").
#
# Outputs (written under <outdir>/pca/):
#   Tables  : pca_scores.tsv  (sample × PC matrix + group)
#             pca_loadings.tsv (gene loadings)
#             variance_explained.tsv (PC, proportion of variance)
#             spearman_correlations.tsv (optional: PC1 vs exposure metrics)
#   Figures : pca_PC1_PC2.png, pca_PC3_PC4.png, pca_scree.png,
#             pc1_loadings_top30.png, pc2_loadings_top30.png
#
# Recommended citations (methods & software):
#   • Love MI, Huber W, Anders S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data
#     with DESeq2. Genome Biology 15:550. doi:10.1186/s13059-014-0550-8
#   • Jolliffe IT, Cadima J. (2016) Principal component analysis: a review and recent developments.
#     Phil. Trans. R. Soc. A 374:20150202. doi:10.1098/rsta.2015.0202
#   • R Core Team. (2025) R: A language and environment for statistical computing. R Foundation for
#     Statistical Computing, Vienna, Austria. https://www.R-project.org/
#   • Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer. doi:10.1007/978-3-319-24277-4
#   • (If used) Zeileis A, Hothorn T. (2002) Diagnostic checking in regression relationships. R News 2(3):7–10.
#     [For 'colors/plotting adjustments'—omit if not relevant.]
#
# Software citation practice:
#   Please also cite exact package versions via sessionInfo() captured at runtime (printed at script end).
#
# Author: Junyuan Guan (Imperial College London)

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(tibble)
  library(DESeq2);  library(ggplot2); library(stringr)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("-c","--counts"), type="character", help="Counts TSV/CSV (genes x samples). First column = gene id"),
  make_option(c("-m","--meta"),   type="character", help="Metadata TSV/CSV with columns: sample_id, group (+ optional exposures)"),
  make_option(c("-o","--outdir"), type="character", default="results", help="Output directory [default %default]"),
  make_option(c("-x","--exclude"),type="character", default=NULL, help="Comma-separated sample IDs to exclude (optional)"),
  make_option(     "--zero_frac", type="double",   default=0.95, help="Filter genes with >zero_frac zeros [default %default]"),
  make_option(     "--pcs",       type="character", default="1,2", help="PCs for main plot, e.g. '1,2' [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$meta)) stop("Please provide --counts and --meta.", call.=FALSE)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
outdir_pca <- file.path(opt$outdir, "pca")
dir.create(outdir_pca, showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(outdir_pca, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
read_table_auto <- function(path){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","tab","txt")) readr::read_tsv(path, show_col_types = FALSE)
  else if (ext %in% c("csv"))        readr::read_csv(path, show_col_types = FALSE)
  else {
    out <- try(readr::read_tsv(path, show_col_types = FALSE), silent=TRUE)
    if (inherits(out, "try-error")) out <- readr::read_csv(path, show_col_types = FALSE)
    out
  }
}
split_csv <- function(x){ if (is.null(x) || is.na(x) || x=="") character(0) else str_trim(unlist(strsplit(x, ","))) }

safe_vst <- function(dds){
  ng <- nrow(dds)
  # If we have a reasonable number of genes, use VST with mean fit and bounded nsub
  if (ng >= 50) {
    return(assay(varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean", nsub = min(1000, ng))))
  }
  # Fallback for tiny matrices: size-factor normalisation + log1p
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  message("Tiny matrix detected (nrow < 50). Using log1p(normalized counts) as a stable fallback.")
  return(log1p(norm_counts))
}

# ---------- read inputs ----------
counts_df <- read_table_auto(opt$counts)
meta_df   <- read_table_auto(opt$meta)

if (!("gene" %in% names(counts_df))) names(counts_df)[1] <- "gene"
stopifnot("gene" %in% names(counts_df))
stopifnot(all(c("sample_id","group") %in% names(meta_df)))

# select and order columns by metadata sample order
sample_ids <- meta_df$sample_id
counts_mat <- counts_df %>% dplyr::select(any_of(c("gene", sample_ids)))
stopifnot(all(sample_ids %in% colnames(counts_mat)))
counts_mat <- counts_mat %>% dplyr::select("gene", dplyr::all_of(sample_ids))

# exclude samples if requested
exc <- split_csv(opt$exclude)
if (length(exc) > 0) {
  keep_samples <- setdiff(sample_ids, exc)
  meta_df   <- meta_df   %>% dplyr::filter(sample_id %in% keep_samples)
  counts_mat<- counts_mat %>% dplyr::select("gene", dplyr::all_of(keep_samples))
  message(sprintf("Excluded samples: %s", paste(exc, collapse=", ")))
}

# ---------- build DESeqDataSet ----------
cts <- as.matrix(counts_mat[,-1]); mode(cts) <- "integer"; rownames(cts) <- counts_mat$gene
# zero filter (relaxed for tiny data)
zero_frac <- rowSums(cts==0)/ncol(cts)
keep <- zero_frac <= opt$zero_frac
if (sum(keep) < 50) { message("Low gene count detected; skipping strict zero filter to stabilise PCA."); keep[] <- TRUE }
cts <- cts[keep, , drop=FALSE]

# align meta rows to column order and coerce group to factor
meta_df <- meta_df[meta_df$sample_id %in% colnames(cts), , drop = FALSE]
meta_df <- meta_df[match(colnames(cts), meta_df$sample_id), , drop = FALSE]
meta_df$group <- as.factor(meta_df$group)
stopifnot(all(meta_df$sample_id == colnames(cts)))

dds <- DESeqDataSetFromMatrix(cts, colData = as.data.frame(meta_df), design = ~ group)

# ---------- transform for PCA (robust) ----------
mat <- safe_vst(dds)          # samples in columns
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

# ---------- variance explained ----------
var_expl <- (pca$sdev^2)/sum(pca$sdev^2)
var_df <- tibble(PC = paste0("PC", seq_along(var_expl)), variance_explained = var_expl)
readr::write_tsv(var_df,   file.path(outdir_pca, "variance_explained.tsv"))

# ---------- scores/loadings ----------
scores <- as.data.frame(pca$x) %>% tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(meta_df %>% dplyr::select(sample_id, group), by="sample_id")
readr::write_tsv(scores,   file.path(outdir_pca, "pca_scores.tsv"))

loadings <- as.data.frame(pca$rotation) %>% tibble::rownames_to_column("gene")
readr::write_tsv(loadings, file.path(outdir_pca, "pca_loadings.tsv"))

# ---------- plots ----------
pcs_main <- as.integer(strsplit(opt$pcs, ",")[[1]]); if (length(pcs_main)!=2) pcs_main <- c(1,2)
pcx <- pcs_main[1]; pcy <- pcs_main[2]
vx <- round(100*var_df$variance_explained[pcx], 1)
vy <- round(100*var_df$variance_explained[pcy], 1)

p_scatter <- ggplot(scores, aes_string(x=paste0("PC",pcx), y=paste0("PC",pcy), color="group")) +
  geom_point(size=2.8, alpha=0.85) +
  theme_minimal(base_size=12) +
  labs(title=sprintf("PCA: PC%d vs PC%d", pcx, pcy),
       x=sprintf("PC%d (%.1f%% var.)", pcx, vx),
       y=sprintf("PC%d (%.1f%% var.)", pcy, vy),
       color="Group")
ggsave(file.path(figdir, sprintf("pca_PC%d_PC%d.png", pcx, pcy)), p_scatter, width=7, height=5, dpi=300)

if (ncol(scores) >= 5) {
  vx34 <- round(100*var_df$variance_explained[3],1); vy34 <- round(100*var_df$variance_explained[4],1)
  p_scatter_34 <- ggplot(scores, aes(x=PC3, y=PC4, color=group)) +
    geom_point(size=2.5, alpha=0.85) + theme_minimal(base_size=12) +
    labs(title="PCA: PC3 vs PC4",
         x=sprintf("PC3 (%.1f%% var.)", vx34),
         y=sprintf("PC4 (%.1f%% var.)", vy34),
         color="Group")
  ggsave(file.path(figdir, "pca_PC3_PC4.png"), p_scatter_34, width=7, height=5, dpi=300)
}

topK <- min(10, length(var_expl))
scree_df <- var_df[1:topK, , drop=FALSE]
p_scree <- ggplot(scree_df, aes(x=factor(PC, levels=PC), y=variance_explained)) +
  geom_col() +
  geom_text(aes(label=sprintf("%.2f", variance_explained)), vjust=-0.3, size=3) +
  theme_minimal(base_size=12) +
  labs(title="Variance explained (top PCs)", x="PC", y="Proportion of variance")
ggsave(file.path(figdir, "pca_scree.png"), p_scree, width=7, height=5, dpi=300)

plot_top_loadings <- function(loadings_df, pc_col, k=30, outfile="pc_loadings.png"){
  ld <- loadings_df %>% dplyr::select(gene, dplyr::all_of(pc_col)) %>%
    dplyr::arrange(dplyr::desc(abs(.data[[pc_col]]))) %>%
    dplyr::slice(1:k) %>% dplyr::mutate(gene=factor(gene, levels=rev(gene)))
  p <- ggplot(ld, aes(x=gene, y=.data[[pc_col]])) +
    geom_col() + coord_flip() + theme_minimal(base_size=11) +
    labs(title=sprintf("Top %d loadings on %s", k, pc_col), x="Gene", y="Loading")
  ggsave(file.path(figdir, outfile), p, width=7, height=6, dpi=300)
}
if ("PC1" %in% names(loadings)) plot_top_loadings(loadings, "PC1", 30, "pc1_loadings_top30.png")
if ("PC2" %in% names(loadings)) plot_top_loadings(loadings, "PC2", 30, "pc2_loadings_top30.png")

# Optional: correlations of PC1 with numeric exposure metrics in metadata
numeric_cols <- names(meta_df)[sapply(meta_df, is.numeric)]
numeric_cols <- setdiff(numeric_cols, c("group"))
if (length(numeric_cols) > 0 && "PC1" %in% names(scores)) {
  merged <- dplyr::left_join(scores, meta_df, by="sample_id")
  cor_rows <- lapply(numeric_cols, function(col){
    x <- merged$PC1; y <- merged[[col]]
    if (all(is.finite(x)) && all(is.finite(y))) {
      ct <- suppressWarnings(cor.test(x, y, method="spearman"))
      tibble(metric = col, rho = unname(ct$estimate), p_value = ct$p.value)
    } else tibble(metric=col, rho=NA_real_, p_value=NA_real_)
  })
  cor_df <- dplyr::bind_rows(cor_rows)
  readr::write_tsv(cor_df, file.path(outdir_pca, "spearman_correlations.tsv"))
}

message("PCA completed. Outputs written to: ", normalizePath(opt$outdir))
message("\n--- sessionInfo() ---")
print(utils::sessionInfo())