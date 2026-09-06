#!/bin/bash
set -euo pipefail

# REGENIE step 2 template.
# GENO_PREFIX_PATTERN must contain {CHR}, for example: /path/to/genotypes_chr{CHR}

: "${TRAIT:?Set TRAIT to the phenotype column name}"
: "${DATA_DIR:?Set DATA_DIR to the analysis working directory}"
: "${GENO_PREFIX_PATTERN:?Set GENO_PREFIX_PATTERN to the chromosome-specific PLINK prefix pattern}"

PHENO_FILE="${PHENO_FILE:-${DATA_DIR}/pheno.txt}"
COVAR_FILE="${COVAR_FILE:-${DATA_DIR}/cov.txt}"
KEEP_FILE="${KEEP_FILE:-${DATA_DIR}/keepIDs.txt}"
PRED_FILE="${PRED_FILE:-${DATA_DIR}/regenie_results/${TRAIT}/step1/${TRAIT}_step1_pred.list}"

STEP2_DIR="${STEP2_DIR:-${DATA_DIR}/regenie_results/${TRAIT}/step2}"
LOG_DIR="${LOG_DIR:-${STEP2_DIR}/logs}"

COVAR_LIST="${COVAR_LIST:-}"
CAT_COVAR_LIST="${CAT_COVAR_LIST:-}"
MAX_CAT_LEVELS="${MAX_CAT_LEVELS:-110}"
TRAIT_TYPE="${TRAIT_TYPE:-quantitative}"
B_SIZE="${B_SIZE:-1000}"
CHROMOSOMES="${CHROMOSOMES:-$(seq 1 22)}"

mkdir -p "${STEP2_DIR}" "${LOG_DIR}"

require_file() {
  local file_path="$1"
  if [[ ! -s "${file_path}" ]]; then
    echo "ERROR: missing or empty file: ${file_path}"
    exit 1
  fi
}

chrom_prefix() {
  local pattern="$1"
  local chr="$2"
  printf "%s\n" "${pattern//\{CHR\}/${chr}}"
}

echo "===== REGENIE Step 2 ====="
echo "Start time: $(date)"
echo "Trait: ${TRAIT}"
echo "REGENIE executable: $(command -v regenie)"
regenie --version

require_file "${PHENO_FILE}"
require_file "${COVAR_FILE}"
require_file "${KEEP_FILE}"
require_file "${PRED_FILE}"

for chr in ${CHROMOSOMES}; do
  GENO_PREFIX="$(chrom_prefix "${GENO_PREFIX_PATTERN}" "${chr}")"
  CHR_OUT_DIR="${STEP2_DIR}/chr${chr}"
  OUT_PREFIX="${CHR_OUT_DIR}/chr${chr}"
  CONSOLE_LOG="${LOG_DIR}/${TRAIT}_chr${chr}.console.log"
  RESULT_FILE="${OUT_PREFIX}_${TRAIT}.regenie"

  mkdir -p "${CHR_OUT_DIR}"

  require_file "${GENO_PREFIX}.bed"
  require_file "${GENO_PREFIX}.bim"
  require_file "${GENO_PREFIX}.fam"

  REGENIE_ARGS=(
    --step 2
    --chr "${chr}"
    --bed "${GENO_PREFIX}"
    --keep "${KEEP_FILE}"
    --phenoFile "${PHENO_FILE}"
    --phenoCol "${TRAIT}"
    --covarFile "${COVAR_FILE}"
    --pred "${PRED_FILE}"
    --bsize "${B_SIZE}"
    --out "${OUT_PREFIX}"
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

  echo "Starting chromosome ${chr}: $(date)"
  regenie "${REGENIE_ARGS[@]}" 2>&1 | tee "${CONSOLE_LOG}"

  if [[ ! -s "${RESULT_FILE}" ]]; then
    echo "ERROR: chromosome ${chr} finished without producing: ${RESULT_FILE}"
    exit 1
  fi

  echo "Chromosome ${chr} completed: $(date)"
done

echo "REGENIE Step 2 completed."
echo "End time: $(date)"
