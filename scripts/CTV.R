#!/usr/bin/env Rscript
# Coefficient/variability (SCBD) visualisation and heatmaps for RNA-seq signatures.
# CLI-friendly, no hard-coded absolute paths; all I/O is parameterised.
#
# Inputs (TSV/CSV; auto-detected by extension):
#   --scbd_raw, --scbd_filtered, --counts_norm, --meta, [--unique_counts, --unique_scbd]
# Key options:
#   --outdir <results>, --exclude, --scbd_threshold, --top_k, --top_heatmap, --remove_globin
#
# Outputs (under <outdir>/ctv/):
#   tables/: SCBD top lists, per-gene medians, expression sums, Wilcoxon tests, ...
#   figures/: barplots, heatmaps, faceted boxplots, expression-sum boxplots, ...
#
# Recommended citations (methods & software):
#   • Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer. doi:10.1007/978-3-319-24277-4
#   • Kolde R. (2019) pheatmap: Pretty Heatmaps. R package v1.0.12. https://CRAN.R-project.org/package=pheatmap
#   • Wickham H., et al. (2019) Welcome to the tidyverse. JOSS 4(43):1686. doi:10.21105/joss.01686
#   • Wilcoxon F. (1945) Individual comparisons by ranking methods. Biometrics Bull. 1:80–83. doi:10.2307/3001968
#     (equivalently: Mann H.B., Whitney D.R. (1947) Ann. Math. Stat. 18:50–60. doi:10.1214/aoms/1177730491)
#   • Murtagh F., Legendre P. (2014) Ward’s hierarchical agglomerative clustering method:
#     which algorithms implement Ward’s criterion? J. Classification 31:274–295. doi:10.1007/s00357-014-9161-z
#   • R Core Team (2025) R: A language and environment for statistical computing. https://www.R-project.org/
#
# Note: If your SCBD metric corresponds to a published definition, add its bibliographic reference here.
# Reproducibility: print sessionInfo() at the end to capture exact R/package versions.
# Author: Junyuan Guan (Imperial College London)

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
  library(stringr)
})

# ---------- CLI ----------
option_list <- list(
  make_option("--scbd_raw",      type="character", help="Path to SCBD table (raw)."),
  make_option("--scbd_filtered", type="character", default=NULL, help="Path to SCBD table (filtered)."),
  make_option("--counts_norm",   type="character", help="Path to normalised count matrix (genes x samples)."),
  make_option("--meta",          type="character", help="Path to metadata (sample IDs + group/condition)."),
  make_option("--unique_counts", type="character", default=NULL, help="(Optional) Unique matrix: normalised counts."),
  make_option("--unique_scbd",   type="character", default=NULL, help="(Optional) Unique matrix: SCBD table."),
  make_option("--outdir",        type="character", default="results", help="Output root [default %default]"),
  make_option("--exclude",       type="character", default=NULL, help="Comma-separated sample IDs to exclude"),
  make_option("--scbd_threshold",type="double",   default=0.0002, help="SCBD threshold [default %default]"),
  make_option("--top_k",         type="integer",  default=50,     help="Top-K features for barplots [default %default]"),
  make_option("--top_heatmap",   type="integer",  default=15,     help="Top-N features for heatmap [default %default]"),
  make_option("--remove_globin", type="logical",  default=TRUE,   help="Remove globin/MB genes [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$scbd_raw) || is.null(opt$counts_norm) || is.null(opt$meta)) {
  stop("Please provide at least --scbd_raw, --counts_norm, and --meta.", call.=FALSE)
}

# ---------- utilities ----------
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

split_csv <- function(x){
  if (is.null(x) || is.na(x) || x=="") character(0) else str_trim(unlist(strsplit(x, ",")))
}

standardise_scbd <- function(df){
  nms <- tolower(names(df))
  pick <- function(cands){ i <- which(nms %in% cands); if (length(i)) i[1] else NA_integer_ }
  i_feat <- pick(c("feature","gene","genes","symbol","id","name"))
  i_scbd <- pick(c("scbd","value","score","scbd_value","scbdscore"))
  if (!is.na(i_feat) && !is.na(i_scbd)) {
    return(tibble::tibble(Feature = as.character(df[[i_feat]]),
                          SCBD    = as.numeric(df[[i_scbd]])))
  }
  if (ncol(df) == 2) {
    is_num <- sapply(df, function(x){
      if (is.numeric(x)) return(TRUE)
      suppressWarnings(all(!is.na(as.numeric(as.character(head(x, 50))))))
    })
    scbd_idx  <- if (any(is_num)) which(is_num)[1] else 2
    feature_idx <- if (any(!is_num)) which(!is_num)[1] else setdiff(1:2, scbd_idx)[1]
    return(tibble::tibble(Feature = as.character(df[[feature_idx]]),
                          SCBD    = as.numeric(df[[scbd_idx]])))
  }
  is_num <- sapply(df, function(x){
    if (is.numeric(x)) return(TRUE)
    suppressWarnings(all(!is.na(as.numeric(as.character(head(x, 50))))))
  })
  if (any(is_num)) {
    scbd_idx <- which(is_num)[1]
    feature_idx <- if (any(!is_num)) which(!is_num)[1] else setdiff(seq_len(ncol(df)), scbd_idx)[1]
    return(tibble::tibble(Feature = as.character(df[[feature_idx]]),
                          SCBD    = as.numeric(df[[scbd_idx]])))
  }
  stop("SCBD table must have 'Feature' and 'SCBD' (or synonyms) or be a two-column table [feature, value].", call.=FALSE)
}

standardise_meta <- function(df){
  nms <- tolower(names(df))
  i_id <- which(nms %in% c("sample","sample_id","id"))[1]
  i_gp <- which(nms %in% c("group","condition","status"))[1]
  if (is.na(i_id) || is.na(i_gp)) stop("Meta must have 'Sample' (or sample_id) and 'group/condition'.", call.=FALSE)
  out <- tibble::tibble(sample_id = as.character(df[[i_id]]),
                        group     = as.character(df[[i_gp]]))
  out$group <- tolower(out$group)
  # map common synonyms to canonical labels
  out$group <- dplyr::case_when(
    grepl("resist", out$group) ~ "Resister",
    grepl("convert", out$group) ~ "Converter",
    grepl("control|hc|healthy", out$group) ~ "HC",
    TRUE ~ tools::toTitleCase(out$group)
  )
  out
}

fix_sample_ids <- function(x){
  # unify: replace dots with dashes to match IDs like TBRHC-015
  gsub("\\.", "-", x)
}

ensure_dir <- function(path){
  dir.create(path, showWarnings = FALSE, recursive = TRUE); path
}

# ---------- I/O dirs ----------
outdir_ctv <- file.path(opt$outdir, "ctv")
figdir     <- ensure_dir(file.path(outdir_ctv, "figures"))
tabdir     <- ensure_dir(file.path(outdir_ctv, "tables"))

# ---------- read inputs ----------
scbd_raw_df      <- standardise_scbd(read_table_auto(opt$scbd_raw))
scbd_filtered_df <- if (!is.null(opt$scbd_filtered)) standardise_scbd(read_table_auto(opt$scbd_filtered)) else NULL
counts_df        <- read_table_auto(opt$counts_norm)
meta_df_raw      <- read_table_auto(opt$meta)

# standardise metadata; fix sample IDs in both meta and counts columns
meta_df <- standardise_meta(meta_df_raw)
colnames(counts_df)[1] <- if (!tolower(colnames(counts_df)[1]) %in% c("gene","genes","feature","id","symbol")) "gene" else colnames(counts_df)[1]
if (tolower(colnames(counts_df)[1]) != "gene") names(counts_df)[1] <- "gene"

# rownames as gene ids; columns as samples
genes <- as.character(counts_df[[1]])
mat   <- as.matrix(counts_df[ , -1, drop=FALSE])
rownames(mat) <- genes
colnames(mat) <- fix_sample_ids(colnames(mat))
meta_df$sample_id <- fix_sample_ids(meta_df$sample_id)

# exclude samples if requested
exc <- split_csv(opt$exclude)
if (length(exc) > 0) {
  exc <- fix_sample_ids(exc)
  keep <- setdiff(colnames(mat), exc)
  mat  <- mat[ , keep, drop=FALSE]
  meta_df <- meta_df %>% dplyr::filter(sample_id %in% colnames(mat))
  message(sprintf("Excluded samples: %s", paste(exc, collapse=", ")))
}

# make sure metadata order matches matrix columns
meta_df <- meta_df %>% dplyr::filter(sample_id %in% colnames(mat))
meta_df <- meta_df[match(colnames(mat), meta_df$sample_id), , drop=FALSE]
stopifnot(all(meta_df$sample_id == colnames(mat)))

# globin removal (applied to any SCBD table we use for plotting)
globin_genes <- c("CYGB","HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ","MB")
maybe_rm_globin <- function(df){
  if (isTRUE(opt$remove_globin)) df <- df %>% dplyr::filter(!(Feature %in% globin_genes))
  df
}

# ---------- plotting helpers ----------
plot_scbd_bars <- function(df, title, outfile, add_hline=NA){
  p <- ggplot(df, aes(x = reorder(Feature, -SCBD), y = SCBD)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "Feature", y = "SCBD") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (is.finite(add_hline)) p <- p + geom_hline(yintercept = add_hline)
  ggsave(outfile, p, width = 10, height = 6, dpi = 300)
}

make_heatmap <- function(submat, annot_df, title, outfile, method="complete"){
  # pheatmap expects annotation rownames to match columns of submat
  ann <- data.frame(condition = factor(annot_df$group, levels = c("HC","Converter","Resister")))
  rownames(ann) <- annot_df$sample_id
  stopifnot(all(colnames(submat) %in% rownames(ann)))
  suppressWarnings(
    pheatmap::pheatmap(
      submat,
      annotation_col = ann[colnames(submat),, drop=FALSE],
      show_rownames = TRUE,
      show_colnames = TRUE,
      scale = "row",
      clustering_method = method,
      main = title,
      fontsize_col = 9,
      fontsize_row = 8,
      filename = outfile,
      width = 10, height = 8
    )
  )
}

# ---------- 1) SCBD raw plots ----------
scbd_raw_df <- maybe_rm_globin(scbd_raw_df) %>% dplyr::arrange(dplyr::desc(SCBD))
readr::write_tsv(head(scbd_raw_df, opt$top_k), file.path(tabdir, "scbd_raw_top50.tsv"))
plot_scbd_bars(scbd_raw_df, "SCBD Values (raw)", file.path(figdir, "scbd_raw_barplot.png"), add_hline = opt$scbd_threshold)
plot_scbd_bars(head(scbd_raw_df, opt$top_k), sprintf("Top %d SCBD Values (raw)", opt$top_k), file.path(figdir, "scbd_raw_top50_barplot.png"), add_hline = opt$scbd_threshold)

# ---------- 2) SCBD filtered plots (if provided) ----------
if (!is.null(scbd_filtered_df)) {
  scbd_filtered_df <- maybe_rm_globin(scbd_filtered_df) %>% dplyr::arrange(dplyr::desc(SCBD))
  readr::write_tsv(head(scbd_filtered_df, opt$top_k), file.path(tabdir, "scbd_filtered_top50.tsv"))
  plot_scbd_bars(scbd_filtered_df, "SCBD Values (filtered)", file.path(figdir, "scbd_filtered_barplot.png"))
  plot_scbd_bars(head(scbd_filtered_df, opt$top_k), sprintf("Top %d SCBD Values (filtered)", opt$top_k), file.path(figdir, "scbd_filtered_top50_barplot.png"))
}

# choose source for heatmap/features: prefer filtered if provided, else raw
scbd_for_heatmap <- if (!is.null(scbd_filtered_df)) scbd_filtered_df else scbd_raw_df
topN <- head(scbd_for_heatmap, opt$top_heatmap)
readr::write_tsv(topN %>% dplyr::select(Feature, SCBD), file.path(tabdir, "heatmap_topN_gene_list.tsv"))

# extract submatrix for topN features
top_genes <- intersect(topN$Feature, rownames(mat))
submat <- mat[top_genes, , drop=FALSE]
# remove excluded again (safety)
if (length(exc) > 0) submat <- submat[, setdiff(colnames(submat), exc), drop=FALSE]

# Heatmap (complete linkage for the main)
make_heatmap(submat, meta_df, sprintf("Heatmap of Top %d Highly Variable Genes", nrow(submat)), file.path(figdir, "heatmap_topN.png"), method="complete")

# ---------- 3) Per-gene distributions & medians by group ----------
long_df <- submat %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="Gene") %>%
  tidyr::pivot_longer(-Gene, names_to="Sample", values_to="Expression") %>%
  dplyr::left_join(meta_df, by=c("Sample"="sample_id"))

# sanity: drop NAs if any
long_df <- long_df %>% dplyr::filter(!is.na(group))
# order facets by gene list
long_df$Gene <- factor(long_df$Gene, levels = top_genes)

medians_df <- long_df %>%
  dplyr::group_by(Gene, group) %>%
  dplyr::summarise(Median_Expression = median(Expression, na.rm = TRUE), .groups="drop")
readr::write_tsv(medians_df, file.path(tabdir, "expression_medians.tsv"))

# per-gene faceted boxplots
p_box <- ggplot(long_df, aes(x = group, y = Expression, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free") +
  theme_minimal(base_size = 11) +
  labs(title = "Distribution of Gene Signatures Between Groups", x = "Group", y = "Expression")
ggsave(file.path(figdir, "boxplots_per_gene.png"), p_box, width = 12, height = 9, dpi = 300)

# ---------- 4) Sum across top genes + Wilcoxon tests ----------
sum_by_sample <- data.frame(
  Sample = colnames(submat),
  ExpressionSum = colSums(submat, na.rm=TRUE),
  stringsAsFactors = FALSE
) %>% dplyr::left_join(meta_df, by=c("Sample"="sample_id"))
readr::write_tsv(sum_by_sample, file.path(tabdir, "expression_sum_by_sample.tsv"))

p_sum <- ggplot(sum_by_sample, aes(x = group, y = ExpressionSum, fill = group)) +
  geom_boxplot() +
  theme_minimal(base_size = 12) +
  labs(title = "Gene Expression Sums by Group", x = "Group", y = "Sum of expression")
ggsave(file.path(figdir, "expression_sum_boxplot.png"), p_sum, width = 7.5, height = 5.5, dpi = 300)

# pairwise Wilcoxon among groups present
groups_present <- unique(sum_by_sample$group)
pairs <- combn(groups_present, 2, simplify = FALSE)
wilc_rows <- lapply(pairs, function(pr){
  a <- pr[1]; b <- pr[2]
  sub <- sum_by_sample %>% dplyr::filter(group %in% c(a,b))
  if (length(unique(sub$group))==2) {
    wt <- suppressWarnings(wilcox.test(ExpressionSum ~ group, data = sub, exact = FALSE))
    tibble::tibble(group_A=a, group_B=b, p_value = wt$p.value, n_A = sum(sub$group==a), n_B = sum(sub$group==b))
  } else {
    tibble::tibble(group_A=a, group_B=b, p_value = NA_real_, n_A = NA_integer_, n_B = NA_integer_)
  }
})
wilc_df <- dplyr::bind_rows(wilc_rows)
readr::write_tsv(wilc_df, file.path(tabdir, "wilcoxon_pairwise.tsv"))

# ---------- 5) "Unique" matrix branch (optional) ----------
if (!is.null(opt$unique_counts) && !is.null(opt$unique_scbd)) {
  u_counts_df <- read_table_auto(opt$unique_counts)
  u_scbd_df   <- standardise_scbd(read_table_auto(opt$unique_scbd)) %>% maybe_rm_globin() %>% dplyr::arrange(dplyr::desc(SCBD))
  
  # normalised unique counts matrix
  colnames(u_counts_df)[1] <- if (!tolower(colnames(u_counts_df)[1]) %in% c("gene","genes","feature","id","symbol")) "gene" else colnames(u_counts_df)[1]
  if (tolower(colnames(u_counts_df)[1]) != "gene") names(u_counts_df)[1] <- "gene"
  u_genes <- as.character(u_counts_df[[1]])
  u_mat   <- as.matrix(u_counts_df[,-1, drop=FALSE]); rownames(u_mat) <- u_genes
  colnames(u_mat) <- fix_sample_ids(colnames(u_mat))
  
  # align and exclude
  keep_s <- intersect(colnames(u_mat), meta_df$sample_id)
  u_mat  <- u_mat[, keep_s, drop=FALSE]
  meta_u <- meta_df %>% dplyr::filter(sample_id %in% keep_s)
  if (length(exc) > 0) {
    keep_s2 <- setdiff(colnames(u_mat), exc)
    u_mat <- u_mat[, keep_s2, drop=FALSE]
    meta_u <- meta_u %>% dplyr::filter(sample_id %in% keep_s2)
  }
  
  # barplots: unique SCBD (full & top50)
  readr::write_tsv(head(u_scbd_df, opt$top_k), file.path(tabdir, "unique_scbd_top50.tsv"))
  plot_scbd_bars(u_scbd_df, "Unique SCBD Values", file.path(figdir, "unique_scbd_barplot.png"), add_hline = opt$scbd_threshold)
  plot_scbd_bars(head(u_scbd_df, opt$top_k), sprintf("Top %d Unique SCBD Values", opt$top_k), file.path(figdir, "unique_scbd_top50_barplot.png"), add_hline = opt$scbd_threshold)
  
  # threshold-based selection & heatmaps
  u_pass <- u_scbd_df %>% dplyr::filter(SCBD >= opt$scbd_threshold)
  if (nrow(u_pass) > 1) {
    u_sub <- u_mat[intersect(u_pass$Feature, rownames(u_mat)), , drop=FALSE]
    # optional exclusion of specific sample (your old code excluded IMTB004); here use --exclude instead
    make_heatmap(u_sub, meta_u, sprintf("Heatmap of %d Highly Variable Genes (unique ≥ %.4f)", nrow(u_sub), opt$scbd_threshold),
                 file.path(figdir, "unique_heatmap_threshold.png"), method="ward.D2")
  }
  
  # heatmap of unique top-50 (if available)
  u_top50 <- head(u_scbd_df, opt$top_k)
  u_topmat <- u_mat[intersect(u_top50$Feature, rownames(u_mat)), , drop=FALSE]
  if (nrow(u_topmat) >= 2) {
    make_heatmap(u_topmat, meta_u, sprintf("Heatmap of Top %d Unique SCBD Genes", nrow(u_topmat)),
                 file.path(figdir, "unique_heatmap_top50.png"), method="complete")
  }
  
  # sum across top50 + wilcoxon
  if (nrow(u_topmat) >= 1) {
    u_sum <- data.frame(Sample = colnames(u_topmat), ExpressionSum = colSums(u_topmat, na.rm=TRUE)) %>%
      dplyr::left_join(meta_u, by=c("Sample"="sample_id"))
    readr::write_tsv(u_sum, file.path(tabdir, "unique_expression_sum_by_sample.tsv"))
    p_u_sum <- ggplot(u_sum, aes(x = group, y = ExpressionSum, fill = group)) +
      geom_boxplot() + theme_minimal(base_size = 12) +
      labs(title = "Gene Expression Sums by Group (unique)", x = "Group", y = "Sum of expression")
    ggsave(file.path(figdir, "unique_expression_sum_boxplot.png"), p_u_sum, width = 7.5, height = 5.5, dpi = 300)
    
    u_groups <- unique(u_sum$group)
    u_pairs <- combn(u_groups, 2, simplify = FALSE)
    u_rows <- lapply(u_pairs, function(pr){
      a <- pr[1]; b <- pr[2]
      sub <- u_sum %>% dplyr::filter(group %in% c(a,b))
      if (length(unique(sub$group))==2) {
        wt <- suppressWarnings(wilcox.test(ExpressionSum ~ group, data = sub, exact = FALSE))
        tibble::tibble(group_A=a, group_B=b, p_value = wt$p.value, n_A = sum(sub$group==a), n_B = sum(sub$group==b))
      } else {
        tibble::tibble(group_A=a, group_B=b, p_value = NA_real_, n_A = NA_integer_, n_B = NA_integer_)
      }
    })
    readr::write_tsv(dplyr::bind_rows(u_rows), file.path(tabdir, "unique_wilcoxon_pairwise.tsv"))
  }
}

message("CTV module complete. Outputs in: ", normalizePath(outdir_ctv))
message("\n--- sessionInfo() ---")
print(utils::sessionInfo())