#!/bin/bash
set -euo pipefail

# REGENIE step 1 template.
# Set required variables before running, for example:
#   TRAIT=trait_name DATA_DIR=/path/to/analysis GENO_PREFIX=/path/to/genotypes SNP_LIST=/path/to/snps.txt bash regenie_step1_template.sh

: "${TRAIT:?Set TRAIT to the phenotype column name}"
: "${DATA_DIR:?Set DATA_DIR to the analysis working directory}"
: "${GENO_PREFIX:?Set GENO_PREFIX to the PLINK bed/bim/fam prefix}"
: "${SNP_LIST:?Set SNP_LIST to the variant list used for step 1}"

PHENO_FILE="${PHENO_FILE:-${DATA_DIR}/pheno.txt}"
COVAR_FILE="${COVAR_FILE:-${DATA_DIR}/cov.txt}"
KEEP_FILE="${KEEP_FILE:-${DATA_DIR}/keepIDs.txt}"

OUT_BASE="${OUT_BASE:-${DATA_DIR}/regenie_results/${TRAIT}}"
OUT_DIR="${OUT_DIR:-${OUT_BASE}/step1}"
TMP_DIR="${TMP_DIR:-${OUT_DIR}/tmp}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"

COVAR_LIST="${COVAR_LIST:-}"
CAT_COVAR_LIST="${CAT_COVAR_LIST:-}"
MAX_CAT_LEVELS="${MAX_CAT_LEVELS:-110}"
TRAIT_TYPE="${TRAIT_TYPE:-quantitative}"
B_SIZE="${B_SIZE:-1000}"

mkdir -p "${TMP_DIR}" "${LOG_DIR}"

RUN_TIME=$(date +%Y%m%d_%H%M%S)
CONSOLE_LOG="${LOG_DIR}/${TRAIT}_step1_${RUN_TIME}.console.log"
exec > >(tee -a "${CONSOLE_LOG}") 2>&1

require_file() {
  local file_path="$1"
  if [[ ! -s "${file_path}" ]]; then
    echo "ERROR: missing or empty file: ${file_path}"
    exit 1
  fi
}

echo "===== REGENIE Step 1 ====="
echo "Start time: $(date)"
echo "Trait: ${TRAIT}"
echo "REGENIE executable: $(command -v regenie)"
regenie --version

require_file "${GENO_PREFIX}.bed"
require_file "${GENO_PREFIX}.bim"
require_file "${GENO_PREFIX}.fam"
require_file "${SNP_LIST}"
require_file "${PHENO_FILE}"
require_file "${COVAR_FILE}"
require_file "${KEEP_FILE}"

REGENIE_ARGS=(
  --step 1
  --bed "${GENO_PREFIX}"
  --extract "${SNP_LIST}"
  --keep "${KEEP_FILE}"
  --phenoFile "${PHENO_FILE}"
  --phenoCol "${TRAIT}"
  --covarFile "${COVAR_FILE}"
  --bsize "${B_SIZE}"
  --lowmem
  --lowmem-prefix "${TMP_DIR}/regenie_l0"
  --out "${OUT_DIR}/${TRAIT}_step1"
)

if [[ -n "${COVAR_LIST}" ]]; then
  REGENIE_ARGS+=(--covarColList "${COVAR_LIST}")
fi

if [[ -n "${CAT_COVAR_LIST}" ]]; then
  REGENIE_ARGS+=(--catCovarList "${CAT_COVAR_LIST}" --maxCatLevels "${MAX_CAT_LEVELS}")
fi

case "${TRAIT_TYPE}" in
  quantitative)
    REGENIE_ARGS+=(--qt)
    ;;
  binary)
    REGENIE_ARGS+=(--bt)
    ;;
  none)
    ;;
  *)
    echo "ERROR: TRAIT_TYPE must be quantitative, binary, or none"
    exit 1
    ;;
esac

echo "Starting REGENIE Step 1..."
regenie "${REGENIE_ARGS[@]}"

echo "REGENIE Step 1 completed."
echo "End time: $(date)"
echo "Prediction list: ${OUT_DIR}/${TRAIT}_step1_pred.list"
