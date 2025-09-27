#!/usr/bin/env bash
# Orchestrator for the TB RNA-seq pipeline (CLI entry point).
# Runs DESeq2 (DGE) and PCA modules; GSEA/CTV typically require extra inputs and
# are meant to be run via their R scripts (see README).
#
# Example:
#   ./scripts/run_all.sh \
#     --counts example_data/counts.tsv \
#     --meta   example_data/metadata.tsv \
#     --outdir results \
#     --ref HC --alt Resister \
#     --exclude TBRHC-015 \
#     --only dge,pca
#
# Tips:
#   - Use --only to run a subset (e.g., dge,pca). Available tokens: dge, pca.
#   - Use --skip to skip modules (e.g., --skip pca).
#   - GSEA and CTV require additional user-provided inputs; run their R scripts
#     explicitly as documented in the README.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_all.sh --counts <file> --meta <file> [--outdir <dir>]
             [--ref <level>] [--alt <level>]
             [--exclude <sample_ids>] [--only <modules>] [--skip <modules>]

Required:
  --counts   Path to counts matrix (genes x samples). First column 'gene'.
  --meta     Path to sample metadata (must contain: sample_id, group).

Optional:
  --outdir   Output root directory (default: results).
  --ref      Reference level for DE contrast (e.g., HC).
  --alt      Comparison level for DE contrast (e.g., Resister).
  --exclude  Comma-separated sample IDs to exclude for PCA (e.g., TBRHC-015).
  --only     Comma-separated modules to run (subset of: dge,pca).
  --skip     Comma-separated modules to skip.

Notes:
  - DGE and PCA are supported here. For GSEA/CTV, call their R scripts directly.
  - See README for full I/O specifications and examples.

Examples:
  ./scripts/run_all.sh --counts example_data/counts.tsv --meta example_data/metadata.tsv --outdir results --ref HC --alt Resister
  ./scripts/run_all.sh -c example_data/counts.tsv -m example_data/metadata.tsv -o results --exclude TBRHC-015 --only pca
USAGE
}

# -------------------------
# Parse CLI arguments
# -------------------------
COUNTS=""
META=""
OUTDIR="results"
REF=""
ALT=""
EXCLUDE=""
ONLY=""
SKIP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--counts)   COUNTS="${2:-}"; shift 2 ;;
    -m|--meta)     META="${2:-}";   shift 2 ;;
    -o|--outdir)   OUTDIR="${2:-}"; shift 2 ;;
    --ref)         REF="${2:-}";    shift 2 ;;
    --alt)         ALT="${2:-}";    shift 2 ;;
    --exclude)     EXCLUDE="${2:-}"; shift 2 ;;
    --only)        ONLY="${2:-}";    shift 2 ;;
    --skip)        SKIP="${2:-}";    shift 2 ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

# Basic validation
if [[ -z "${COUNTS}" || -z "${META}" ]]; then
  echo "[FATAL] --counts and --meta are required." >&2
  usage
  exit 2
fi

# Resolve repository root (works when invoked from project root or elsewhere)
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Create output root and start a log
mkdir -p "${OUTDIR}"
LOG="${OUTDIR%/}/run.log"
{
  echo "===== TB RNA-seq pipeline ====="
  echo "Start time: $(date)"
  echo "ROOT_DIR  : ${ROOT_DIR}"
  echo "COUNTS    : ${COUNTS}"
  echo "META      : ${META}"
  echo "OUTDIR    : ${OUTDIR}"
  echo "REF/ALT   : ${REF:+$REF}/${ALT:+$ALT}"
  echo "EXCLUDE   : ${EXCLUDE:-<none>}"
  echo "ONLY/SKIP : ${ONLY:-<none>}/${SKIP:-<none>}"
} | tee "${LOG}"

# Helper: determine whether a module should run considering --only/--skip
want() {
  local module="$1"                    # e.g., dge or pca
  if [[ -n "${ONLY}" ]]; then
    [[ ",${ONLY}," == *",${module},"* ]] || return 1
  fi
  if [[ -n "${SKIP}" ]]; then
    [[ ",${SKIP}," == *",${module},"* ]] && return 1
  fi
  return 0
}

# Trap for friendly error messages
trap 'echo "[ERROR] Aborted at line $LINENO"; exit 1' ERR

# -------------------------
# DGE (DESeq2)
# -------------------------
if want "dge"; then
  echo "[DGE] Running DESeq2..." | tee -a "${LOG}"
  Rscript "${ROOT_DIR}/scripts/DGE.R" \
    --counts "${COUNTS}" \
    --meta   "${META}" \
    --outdir "${OUTDIR}" \
    ${REF:+--ref "${REF}"} \
    ${ALT:+--alt "${ALT}"} \
    | tee -a "${LOG}"
  echo "[DGE] Done." | tee -a "${LOG}"
else
  echo "[DGE] Skipped." | tee -a "${LOG}"
fi

# -------------------------
# PCA (VST → PCA)
# -------------------------
if want "pca"; then
  echo "[PCA] Running VST → PCA..." | tee -a "${LOG}"
  Rscript "${ROOT_DIR}/scripts/PCA.R" \
    --counts "${COUNTS}" \
    --meta   "${META}" \
    --outdir "${OUTDIR}" \
    ${EXCLUDE:+--exclude "${EXCLUDE}"} \
    --pcs 1,2 \
    | tee -a "${LOG}"
  echo "[PCA] Done." | tee -a "${LOG}"
else
  echo "[PCA] Skipped." | tee -a "${LOG}"
fi

# -------------------------
# Final message
# -------------------------
{
  echo "Finished: $(date)"
  echo "Outputs are under: ${OUTDIR%/}/"
  echo "  - DGE : ${OUTDIR%/}/dge/{tables,figures}"
  echo "  - PCA : ${OUTDIR%/}/pca/{figures,*.tsv}"
  echo
  echo "Note: GSEA and CTV are not run here; invoke their R scripts explicitly (see README)."
} | tee -a "${LOG}"