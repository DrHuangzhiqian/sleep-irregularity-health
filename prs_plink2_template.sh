#!/bin/bash
set -euo pipefail

# PLINK2 PRS scoring template.
# GENO_PREFIX_PATTERN and WEIGHT_PATTERN must contain {CHR}.

: "${PLINK2:?Set PLINK2 to the plink2 executable}"
: "${GENO_PREFIX_PATTERN:?Set GENO_PREFIX_PATTERN to the chromosome-specific PLINK prefix pattern}"
: "${KEEP_FILE:?Set KEEP_FILE to the target sample keep file}"
: "${WEIGHT_PATTERN:?Set WEIGHT_PATTERN to the chromosome-specific weight file pattern}"
: "${SCORE_DIR:?Set SCORE_DIR to the score output directory}"

CHROMOSOMES="${CHROMOSOMES:-$(seq 1 22)}"
SCORE_NAME="${SCORE_NAME:-prs_score}"
LOG_DIR="${LOG_DIR:-${SCORE_DIR}/logs}"

mkdir -p "${SCORE_DIR}" "${LOG_DIR}"

RUN_TIME=$(date +%Y%m%d_%H%M%S)
CONSOLE_LOG="${LOG_DIR}/${SCORE_NAME}_${RUN_TIME}.console.log"
exec > >(tee -a "${CONSOLE_LOG}") 2>&1

require_file() {
  local file_path="$1"
  if [[ ! -s "${file_path}" ]]; then
    echo "ERROR: missing or empty file: ${file_path}"
    exit 1
  fi
}

require_executable() {
  local file_path="$1"
  if [[ ! -x "${file_path}" ]]; then
    echo "ERROR: not executable: ${file_path}"
    exit 1
  fi
}

chrom_path() {
  local pattern="$1"
  local chr="$2"
  printf "%s\n" "${pattern//\{CHR\}/${chr}}"
}

require_executable "${PLINK2}"
require_file "${KEEP_FILE}"

echo "===== PLINK2 PRS scoring ====="
echo "Start time: $(date)"
"${PLINK2}" --version

SCORE_FILES=()

for chr in ${CHROMOSOMES}; do
  GENO_PREFIX="$(chrom_path "${GENO_PREFIX_PATTERN}" "${chr}")"
  WEIGHT_FILE="$(chrom_path "${WEIGHT_PATTERN}" "${chr}")"
  VARIANT_LIST="${SCORE_DIR}/${SCORE_NAME}_chr${chr}.score_ids.txt"
  OUT_PREFIX="${SCORE_DIR}/${SCORE_NAME}_chr${chr}"

  require_file "${GENO_PREFIX}.bed"
  require_file "${GENO_PREFIX}.bim"
  require_file "${GENO_PREFIX}.fam"
  require_file "${WEIGHT_FILE}"

  awk 'NR > 1 {print $1}' "${WEIGHT_FILE}" > "${VARIANT_LIST}"
  require_file "${VARIANT_LIST}"

  echo "Scoring chromosome ${chr}: $(date)"
  "${PLINK2}" \
    --bfile "${GENO_PREFIX}" \
    --keep "${KEEP_FILE}" \
    --extract "${VARIANT_LIST}" \
    --score "${WEIGHT_FILE}" 1 2 3 header-read cols=fid,nallele,scoresums \
    --out "${OUT_PREFIX}"

  SSCORE="${OUT_PREFIX}.sscore"
  require_file "${SSCORE}"
  SCORE_FILES+=("${SSCORE}")
done

FINAL_SCORE="${SCORE_DIR}/${SCORE_NAME}_PRS_raw.txt"

awk '
BEGIN {
  OFS = "\t"
}

FNR == 1 {
  fid_col = 0
  iid_col = 0
  score_col = 0

  for (i = 1; i <= NF; i++) {
    if ($i == "#FID" || $i == "FID")
      fid_col = i
    if ($i == "IID")
      iid_col = i
    if ($i ~ /_SUM$/)
      score_col = i
  }

  if (fid_col == 0 || iid_col == 0 || score_col == 0) {
    print "ERROR: required columns not found in " FILENAME > "/dev/stderr"
    exit 1
  }

  next
}

{
  fid = $fid_col
  iid = $iid_col
  score = $score_col
  key = fid SUBSEP iid

  if (!(key in observed)) {
    observed[key] = 1
    order[++sample_count] = key
    saved_fid[key] = fid
    saved_iid[key] = iid
  }

  prs_sum[key] += score
  chromosome_count[key]++
}

END {
  print "FID", "IID", "PRS_raw"

  for (i = 1; i <= sample_count; i++) {
    key = order[i]

    if (chromosome_count[key] != 22) {
      print "ERROR: incomplete chromosome count for sample " saved_fid[key], saved_iid[key] > "/dev/stderr"
      error_found = 1
      continue
    }

    print saved_fid[key], saved_iid[key], prs_sum[key]
  }

  if (error_found)
    exit 1
}
' "${SCORE_FILES[@]}" > "${FINAL_SCORE}"

require_file "${FINAL_SCORE}"

echo "PRS scoring completed."
echo "End time: $(date)"
echo "Final PRS file: ${FINAL_SCORE}"
echo "Console log: ${CONSOLE_LOG}"
