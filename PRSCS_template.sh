#!/bin/bash
#SBATCH --job-name=prscs_auto
#SBATCH --array=1-22
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --output=logs/prscs_%A_chr%a.out
#SBATCH --error=logs/prscs_%A_chr%a.err

set -euo pipefail

# PRS-CS-auto template.
# Create the logs directory before sbatch, or adjust the SBATCH output paths.

: "${PRSCS_PY:?Set PRSCS_PY to PRScs.py}"
: "${REF_DIR:?Set REF_DIR to the LD reference directory}"
: "${SUMSTATS:?Set SUMSTATS to the cleaned summary statistics file}"
: "${BIM_PREFIX:?Set BIM_PREFIX to the target PLINK bim prefix}"
: "${OUT_DIR:?Set OUT_DIR to the PRS-CS output directory}"
: "${N_GWAS:?Set N_GWAS to the GWAS sample size}"

CHR="${SLURM_ARRAY_TASK_ID:?Run as a Slurm array job or set SLURM_ARRAY_TASK_ID}"
OUT_NAME="${OUT_NAME:-prs_trait}"
SEED_BASE="${SEED_BASE:-100000}"
SNPINFO_FILE="${SNPINFO_FILE:-}"
LD_FILE_PATTERN="${LD_FILE_PATTERN:-}"

if [[ -n "${MODULE_NAME:-}" ]]; then
  module load "${MODULE_NAME}"
fi

if [[ -n "${CONDA_ENV:-}" ]]; then
  conda activate "${CONDA_ENV}"
fi

require_file() {
  local file_path="$1"
  if [[ ! -s "${file_path}" ]]; then
    echo "ERROR: missing or empty file: ${file_path}"
    exit 1
  fi
}

chrom_path() {
  local pattern="$1"
  local chr="$2"
  printf "%s\n" "${pattern//\{CHR\}/${chr}}"
}

echo "===== PRS-CS-auto ====="
echo "Start time: $(date)"
echo "Chromosome: ${CHR}"
echo "Python: $(command -v python)"
python --version

require_file "${PRSCS_PY}"
require_file "${SUMSTATS}"
require_file "${BIM_PREFIX}.bim"

if [[ -n "${SNPINFO_FILE}" ]]; then
  require_file "${SNPINFO_FILE}"
fi

if [[ -n "${LD_FILE_PATTERN}" ]]; then
  require_file "$(chrom_path "${LD_FILE_PATTERN}" "${CHR}")"
fi

mkdir -p "${OUT_DIR}"

python -c '
import h5py
import numpy
import scipy
print("PRS-CS Python dependencies OK")
'

python "${PRSCS_PY}" \
  --ref_dir="${REF_DIR}" \
  --bim_prefix="${BIM_PREFIX}" \
  --sst_file="${SUMSTATS}" \
  --n_gwas="${N_GWAS}" \
  --chrom="${CHR}" \
  --seed="$((SEED_BASE + CHR))" \
  --out_dir="${OUT_DIR}/${OUT_NAME}"

echo "Chromosome ${CHR} completed."
echo "End time: $(date)"
