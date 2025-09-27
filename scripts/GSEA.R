#!/usr/bin/env Rscript
# GSEA dot-plot renderer (CLI-friendly, path-agnostic).
# Reads one or multiple GSEA summary tables and produces bubble plots.
#
# Expected columns (case-insensitive; synonyms handled):
#   - Geneset name : Genes | pathway | term | geneset | set | gs
#   - NES         : NES | nes
#   - FDR q-value : FDR_q_val | qval | q_value | q | padj (use the most appropriate)
#   - Size        : Size | size | n | k
#
# Inputs:
#   --tables   Comma-separated paths to TSV/CSV, e.g. "gsea_merged6.tsv,gsea_medicus.tsv,gsea_legacy.tsv"
#   --labels   Optional comma-separated labels for titles/filenames (must match table count)
#   --outdir   Output root (default: results); module outputs go to <outdir>/gsea/
#   --basename Optional base name for combined outputs (default: "gsea_all")
#
# Outputs (under <outdir>/gsea/):
#   - <label>_table.tsv           : cleaned/sorted table per input
#   - figures/<label>_dotplot.png : bubble plot per input
#   - gsea_all_table.tsv          : concatenated table if ≥2 inputs
#   - figures/gsea_all_dotplot.png: faceted plot if ≥2 inputs
#
# Recommended citations (methods & resources):
#   • Subramanian A, et al. (2005) Gene set enrichment analysis. PNAS 102:15545–15550. doi:10.1073/pnas.0506580102
#   • Liberzon A, et al. (2015) The MSigDB hallmark gene set collection. Cell Systems 1:417–425. doi:10.1016/j.cels.2015.12.004
#   • Sergushichev A. (2016) Fast preranked GSEA via cumulative statistic. bioRxiv 060012. doi:10.1101/060012  (if fgsea used upstream)
#   • Kanehisa M, Goto S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res 28:27–30. doi:10.1093/nar/28.1.27
#   • Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer. doi:10.1007/978-3-319-24277-4
#   • Wickham H., et al. (2019) Welcome to the tidyverse. JOSS 4(43):1686. doi:10.21105/joss.01686
#
# Reproducibility: print sessionInfo() at the end to capture exact R/package versions.
# Note: If using MSigDB/KEGG gene-set names, ensure compliance with their terms of use.
# Author: Junyuan Guan (Imperial College London)

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("-t","--tables"),  type="character", help="Comma-separated TSV/CSV paths for GSEA summaries"),
  make_option(c("-l","--labels"),  type="character", default=NULL, help="Comma-separated labels (same length as --tables)"),
  make_option(c("-o","--outdir"),  type="character", default="results", help="Output root [default %default]"),
  make_option(     "--basename",   type="character", default="gsea_all", help="Base name for combined outputs [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$tables)) stop("Please provide --tables with one or more TSV/CSV paths.", call.=FALSE)

# ---------- utils ----------
split_csv <- function(x){
  if (is.null(x) || is.na(x) || x=="") character(0) else str_trim(unlist(strsplit(x, ",")))
}
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
standardise_cols <- function(df){
  # lower-case for matching
  nms <- tolower(names(df))
  std <- function(targets){
    i <- which(nms %in% targets)
    if (length(i)==0) return(NA_integer_) else return(i[1])
  }
  i_name <- std(c("genes","pathway","term","geneset","set","gs","description","name","id","ids"))
  i_nes  <- std(c("nes","normalized_enrichment_score"))
  i_q    <- std(c("fdr_q_val","qval","q_value","q","fdr","padj","fdr.q.val"))
  i_size <- std(c("size","n","k"))
  if (is.na(i_name) || is.na(i_nes) || is.na(i_q) || is.na(i_size)) {
    stop("Could not find required columns (Genes/NES/FDR_q_val/Size or synonyms).", call.=FALSE)
  }
  out <- tibble::tibble(
    Genes = df[[i_name]],
    NES   = as.numeric(df[[i_nes]]),
    FDR_q_val = as.numeric(df[[i_q]]),
    Size  = as.numeric(df[[i_size]])
  )
  out
}

sanitise <- function(x){
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

plot_dot <- function(df, title, outfile){
  p <- ggplot(df, aes(x = NES, y = reorder(Genes, NES), size = Size, color = FDR_q_val)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "NES", y = "", color = "FDR q") +
    theme(axis.text.y = element_text(size = 9), axis.title.x = element_text(size = 11))
  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
}

# ---------- I/O setup ----------
outdir_gsea <- file.path(opt$outdir, "gsea")
dir.create(outdir_gsea, showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(outdir_gsea, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# ---------- main ----------
paths  <- split_csv(opt$tables)
labels <- split_csv(opt$labels)
if (length(labels)==0) labels <- basename(paths)
if (length(labels) != length(paths)) stop("--labels length must match --tables.", call.=FALSE)

all_list <- list()
for (k in seq_along(paths)) {
  path <- paths[k]
  lab  <- labels[k]
  df_raw <- read_table_auto(path)
  df <- standardise_cols(df_raw) %>%
    dplyr::filter(is.finite(NES), is.finite(FDR_q_val), is.finite(Size)) %>%
    dplyr::arrange(NES)
  # write cleaned/sorted table
  lab_s <- sanitise(lab)
  out_tab <- file.path(outdir_gsea, sprintf("%s_table.tsv", lab_s))
  readr::write_tsv(df, out_tab)
  # plot
  out_png <- file.path(figdir, sprintf("%s_dotplot.png", lab_s))
  plot_dot(df, sprintf("GSEA — %s", lab), out_png)
  # collect
  df$Dataset <- lab
  all_list[[k]] <- df
}

# Combined facet (if multiple)
if (length(all_list) >= 2) {
  combined <- dplyr::bind_rows(all_list)
  readr::write_tsv(combined, file.path(outdir_gsea, sprintf("%s_table.tsv", sanitise(opt$basename))))
  p_all <- ggplot(combined, aes(x = NES, y = reorder(Genes, NES), size = Size, color = FDR_q_val)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal(base_size = 12) +
    labs(title = "GSEA (all datasets)", x = "NES", y = "", color = "FDR q") +
    theme(axis.text.y = element_text(size = 7)) +
    facet_wrap(~ Dataset, scales = "free_y")
  ggsave(file.path(figdir, sprintf("%s_dotplot.png", sanitise(opt$basename))), p_all, width = 10, height = 7, dpi = 300)
}

message("GSEA plots written to: ", normalizePath(figdir))