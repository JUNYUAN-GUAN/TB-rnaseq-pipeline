#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# RNA-seq pipeline runner
# Modules: DGE → PCA → (GSEA plot) → (CTV)
# -----------------------------

usage() {
  cat <<'USAGE'
Usage:
  scripts/run_all.sh [OPTIONS]

Core inputs:
  -c, --counts FILE           Counts matrix (genes x samples). Default: example_data/counts.tsv
  -m, --meta   FILE           Metadata with columns: sample_id, group. Default: example_data/metadata.tsv
  -o, --outdir DIR            Output root. Default: results

DGE options:
      --ref STR               Reference group label (after normalization, e.g. "HC")
      --alt STR               Comparison group label (e.g. "Resister")
      --alpha FLOAT           FDR threshold. Default: 0.05
      --zero-frac FLOAT       Zero proportion filter. Default: 0.95

PCA options:
      --exclude IDS           Comma-separated sample IDs to exclude (e.g. TBRHC-015)
      --pcs STR               PCs for main plot, e.g. "1,2". Default: 1,2

GSEA (dot-plot) options:
      --gsea-tables LIST      Comma-separated TSV/CSV paths
      --gsea-labels LIST      Comma-separated labels (same length as tables)
      --gsea-basename STR     Basename for combined outputs. Default: gsea_all

CTV options:
      --scbd-raw FILE         SCBD table (raw). Two cols: Feature, SCBD (or synonyms)
      --scbd-filtered FILE    SCBD table (filtered)
      --counts-norm FILE      Normalized counts for CTV. Default: results/dge/tables/normalized_counts.tsv
      --unique-counts FILE    (optional) Unique normalized counts
      --unique-scbd FILE      (optional) Unique SCBD table
      --scbd-threshold FLOAT  SCBD threshold. Default: 0.0002
      --top-k INT             Top-K features for barplots. Default: 50
      --top-heatmap INT       Top-N features for heatmap. Default: 15
      --no-remove-globin      Do NOT remove globin/MB genes (default is remove)

Control:
      --only LIST             Run only these modules (comma list from: dge,pca,gsea,ctv)
      --skip LIST             Skip these modules (comma list)
  -h,  --help                 Show this help

Examples:
  # Full run (DGE→PCA, and if GSEA/CTV inputs provided, run them too)
  scripts/run_all.sh -c example_data/counts.tsv -m example_data/metadata.tsv -o results --ref HC --alt Resister

  # Run only DGE + PCA with an exclusion
  scripts/run_all.sh -c counts.tsv -m metadata.tsv --exclude TBRHC-015 --only dge,pca

  # Add GSEA dot-plots
  scripts/run_all.sh --gsea-tables results/gsea/gsea_merged6.tsv,results/gsea/gsea_legacy.tsv --gsea-labels Merged6,Legacy
USAGE
}

# -------- Defaults --------
COUNTS="example_data/counts.tsv"
META="example_data/metadata.tsv"
OUTDIR="results"

REF=""
ALT=""
ALPHA="0.05"
ZERO_FRAC="0.95"

EXCLUDE=""
PCS="1,2"

GSEA_TABLES=""
GSEA_LABELS=""
GSEA_BASENAME="gsea_all"

SCBD_RAW=""
SCBD_FILTERED=""
COUNTS_NORM=""
UNIQUE_COUNTS=""
UNIQUE_SCBD=""
SCBD_THRESHOLD="0.0002"
TOP_K="50"
TOP_HEATMAP="15"
REMOVE_GLOBIN="TRUE"  # default remove

ONLY=""
SKIP=""

# -------- Parse args --------
args=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--counts)        COUNTS="$2"; shift 2;;
    -m|--meta)          META="$2"; shift 2;;
    -o|--outdir)        OUTDIR="$2"; shift 2;;
    --ref)              REF="${2:-}"; shift 2;;
    --alt)              ALT="${2:-}"; shift 2;;
    --alpha)            ALPHA="$2"; shift 2;;
    --zero-frac)        ZERO_FRAC="$2"; shift 2;;
    --exclude)          EXCLUDE="$2"; shift 2;;
    --pcs)              PCS="$2"; shift 2;;
    --gsea-tables)      GSEA_TABLES="$2"; shift 2;;
    --gsea-labels)      GSEA_LABELS="${2:-}"; shift 2;;
    --gsea-basename)    GSEA_BASENAME="$2"; shift 2;;
    --scbd-raw)         SCBD_RAW="$2"; shift 2;;
    --scbd-filtered)    SCBD_FILTERED="$2"; shift 2;;
    --counts-norm)      COUNTS_NORM="$2"; shift 2;;
    --unique-counts)    UNIQUE_COUNTS="$2"; shift 2;;
    --unique-scbd)      UNIQUE_SCBD="$2"; shift 2;;
    --scbd-threshold)   SCBD_THRESHOLD="$2"; shift 2;;
    --top-k)            TOP_K="$2"; shift 2;;
    --top-heatmap)      TOP_HEATMAP="$2"; shift 2;;
    --no-remove-globin) REMOVE_GLOBIN="FALSE"; shift 1;;
    --only)             ONLY="$2"; shift 2;;
    --skip)             SKIP="$2"; shift 2;;
    -h|--help)          usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

# -------- Helpers --------
die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }
have() { [[ -n "${1:-}" ]]; }
in_list() { [[ ",$1," == *",$2,"* ]]; }   # in_list "a,b,c" "b"

# Ensure deps
need Rscript

mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
exec > >(tee -a "$LOG") 2>&1
echo ">>> $(date) – Starting pipeline"
echo "Args: ${args[*]}"

# Resolve normalized counts for CTV (default from DGE)
if [[ -z "$COUNTS_NORM" ]]; then
  COUNTS_NORM="$OUTDIR/dge/tables/normalized_counts.tsv"
fi

# Which modules to run?
RUN_DGE=1; RUN_PCA=1; RUN_GSEA=0; RUN_CTV=0
# Enable GSEA/CTV only if inputs provided
have "$GSEA_TABLES" && RUN_GSEA=1
have "$SCBD_RAW" && have "$COUNTS_NORM" && have "$META" && RUN_CTV=1

# Apply --only/--skip
if have "$ONLY"; then
  RUN_DGE=0; RUN_PCA=0; RUN_GSEA=0; RUN_CTV=0
  IFS=',' read -r -a arr <<<"$ONLY"
  for m in "${arr[@]}"; do
    case "$m" in
      dge)  RUN_DGE=1;;
      pca)  RUN_PCA=1;;
      gsea) RUN_GSEA=1;;
      ctv)  RUN_CTV=1;;
      *) die "--only contains unknown module: $m (valid: dge,pca,gsea,ctv)";;
    esac
  done
fi
if have "$SKIP"; then
  IFS=',' read -r -a arr <<<"$SKIP"
  for m in "${arr[@]}"; do
    case "$m" in
      dge)  RUN_DGE=0;;
      pca)  RUN_PCA=0;;
      gsea) RUN_GSEA=0;;
      ctv)  RUN_CTV=0;;
      *) die "--skip contains unknown module: $m (valid: dge,pca,gsea,ctv)";;
    esac
  done
fi

# Check files exist
[[ -f "$COUNTS" ]] || die "Counts not found: $COUNTS"
[[ -f "$META"   ]] || die "Meta not found: $META"

# -------- Module runners --------
run_dge() {
  [[ "$RUN_DGE" -eq 1 ]] || { echo "[DGE] skipped"; return; }
  [[ -f scripts/DGE.R ]] || { echo "[DGE] scripts/DGE.R not found, skipping"; return; }
  echo ">>> $(date) – DGE"
  Rscript scripts/DGE.R \
    --counts "$COUNTS" \
    --meta   "$META" \
    --outdir "$OUTDIR" \
    ${REF:+--ref "$REF"} \
    ${ALT:+--alt "$ALT"} \
    --alpha "$ALPHA" \
    --zero_frac "$ZERO_FRAC"
}

run_pca() {
  [[ "$RUN_PCA" -eq 1 ]] || { echo "[PCA] skipped"; return; }
  [[ -f scripts/PCA.R ]] || { echo "[PCA] scripts/PCA.R not found, skipping"; return; }
  echo ">>> $(date) – PCA"
  Rscript scripts/PCA.R \
    --counts "$COUNTS" \
    --meta   "$META" \
    --outdir "$OUTDIR" \
    ${EXCLUDE:+--exclude "$EXCLUDE"} \
    --pcs "$PCS"
}

run_gsea() {
  [[ "$RUN_GSEA" -eq 1 ]] || { echo "[GSEA] skipped"; return; }
  [[ -f scripts/GSEA.R ]] || { echo "[GSEA] scripts/GSEA.R not found, skipping"; return; }
  echo ">>> $(date) – GSEA (dot-plots)"
  Rscript scripts/GSEA.R \
    --tables "$GSEA_TABLES" \
    ${GSEA_LABELS:+--labels "$GSEA_LABELS"} \
    --outdir "$OUTDIR" \
    --basename "$GSEA_BASENAME"
}

run_ctv() {
  [[ "$RUN_CTV" -eq 1 ]] || { echo "[CTV] skipped"; return; }
  [[ -f scripts/CTV.R ]] || { echo "[CTV] scripts/CTV.R not found, skipping"; return; }

  # If normalized counts not yet produced (no DGE run and default path), warn
  if [[ "$COUNTS_NORM" == "$OUTDIR/dge/tables/normalized_counts.tsv" && ! -f "$COUNTS_NORM" ]]; then
    echo "[CTV] normalized_counts not found at $COUNTS_NORM; you can provide --counts-norm to override."
  fi

  echo ">>> $(date) – CTV"
  Rscript scripts/CTV.R \
    --scbd_raw      "$SCBD_RAW" \
    ${SCBD_FILTERED:+--scbd_filtered "$SCBD_FILTERED"} \
    --counts_norm   "$COUNTS_NORM" \
    --meta          "$META" \
    ${UNIQUE_COUNTS:+--unique_counts "$UNIQUE_COUNTS"} \
    ${UNIQUE_SCBD:+--unique_scbd "$UNIQUE_SCBD"} \
    --outdir        "$OUTDIR" \
    ${EXCLUDE:+--exclude "$EXCLUDE"} \
    --scbd_threshold "$SCBD_THRESHOLD" \
    --top_k          "$TOP_K" \
    --top_heatmap    "$TOP_HEATMAP" \
    --remove_globin  "$REMOVE_GLOBIN"
}

# -------- Execute --------
run_dge
run_pca
run_gsea
run_ctv

echo ">>> $(date) – Done."
echo "Results under: $OUTDIR/"