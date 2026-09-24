#!/usr/bin/env sh
# ============================================================================
# run_all.sh — one-shot pipeline for the GSE221492 validation module.
#
# Order matters: the two implementations are run first, then asserted to agree,
# then the manuscript's reported values are asserted against both, then the
# repository is checked for artefacts that must not be published.
#
# Any failing assertion exits non-zero.
#
# Usage:  sh run_all.sh
# ============================================================================
set -eu

ROOT=$(cd "$(dirname "$0")" && pwd)
export LPS_ROOT="$ROOT"
cd "$ROOT"

echo "=============================================================="
echo "LipidScore_DDLPS — full pipeline"
echo "  repo root: $ROOT"
echo "=============================================================="

DATA="$ROOT/data/GSE221492_bulk/raw/GSE221492_bulk_LPS_atlas_raw_counts_matrix.tsv.gz"
if [ ! -f "$DATA" ]; then
  echo
  echo "❌ source data not found:"
  echo "     $DATA"
  echo
  echo "   The data/ tree is excluded from version control (third-party data)."
  echo "   Download GSE221492 as described in DATA_SOURCES.md, then re-run."
  exit 1
fi

run() {
  echo
  echo "--------------------------------------------------------------"
  echo ">>> $*"
  echo "--------------------------------------------------------------"
  "$@"
}

# 1) R implementation (independent re-derivation)
run Rscript scripts/07_bulk_validation/01_lipid_score_bulk.R

# 2) Python implementation (origin of the reported numbers)
run python scripts/07_bulk_validation/01_lipid_score_gse221492.py

# 3) cross-validation — implementations must agree, and must match the manuscript
run python scripts/09_qc/assert_manuscript_numbers.py --root "$ROOT"

# 4) repository hygiene — nothing unpublished may be tracked
run python scripts/00_repo_hygiene/check_hygiene.py --root "$ROOT"

echo
echo "=============================================================="
echo "✅ pipeline complete — all assertions passed"
echo "   generated: results/tables/, results/objects/, results/figures/"
echo "   (not tracked: see .gitignore)"
echo "=============================================================="
