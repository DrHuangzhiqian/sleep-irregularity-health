#!/usr/bin/env bash
# Batch template for running get_sleep on CWA files with Slurm.
#
# Required manifest column:
#   cwa_path
#
# Optional manifest columns:
#   file_id, filename, asleep_id

#SBATCH --job-name=asleep_batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/asleep_%A_%a.out
#SBATCH --error=logs/asleep_%A_%a.err

set -euo pipefail

: "${MANIFEST:?Set MANIFEST to the input CSV manifest}"
: "${OUT_ROOT:?Set OUT_ROOT to the output directory}"
: "${LOG_ROOT:?Set LOG_ROOT to the log directory}"

DEVICE="${DEVICE:-cuda:0}"
MIN_WEAR="${MIN_WEAR:-22}"
PYTHON_BIN="${PYTHON_BIN:-python}"
GET_SLEEP_BIN="${GET_SLEEP_BIN:-get_sleep}"
LOCAL_REPO_PATH="${LOCAL_REPO_PATH:-}"
MODEL_WEIGHT_PATH="${MODEL_WEIGHT_PATH:-}"
WORK_ROOT="${WORK_ROOT:-${SLURM_TMPDIR:-/tmp}/asleep_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-1}}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUT_ROOT}" "${LOG_ROOT}" "${WORK_ROOT}"

cleanup() {
  rm -rf "${WORK_ROOT}"
}
trap cleanup EXIT

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

ROW_JSON="$("${PYTHON_BIN}" - "${MANIFEST}" "${TASK_ID}" <<'PY'
import csv
import json
import sys

manifest = sys.argv[1]
task_id = int(sys.argv[2])

with open(manifest, newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

if task_id < 1 or task_id > len(rows):
    raise SystemExit(f"task_id {task_id} outside manifest rows 1..{len(rows)}")

row = rows[task_id - 1]
if not row.get("cwa_path"):
    raise SystemExit("manifest must include a non-empty cwa_path column")

print(json.dumps(row))
PY
)"

CWA_PATH="$("${PYTHON_BIN}" - <<PY
import json
row = json.loads('''${ROW_JSON}''')
print(row["cwa_path"])
PY
)"

CWA_STEM="$("${PYTHON_BIN}" - <<PY
from pathlib import Path

path = Path(r'''${CWA_PATH}''')
name = path.name

if name.endswith(".cwa.gz"):
    name = name[:-7]
elif name.endswith(".cwa"):
    name = name[:-4]
else:
    name = path.stem

print(name)
PY
)"

ASLEEP_ID="$("${PYTHON_BIN}" - <<PY
import json

row = json.loads('''${ROW_JSON}''')
value = row.get("asleep_id") or row.get("file_id") or row.get("filename") or "${CWA_STEM}"

if value.endswith(".cwa.gz"):
    value = value[:-7]
elif value.endswith(".cwa"):
    value = value[:-4]

print(value)
PY
)"

FINAL_DIR="${OUT_ROOT}/${ASLEEP_ID}"
RUN_LOG="${LOG_ROOT}/${ASLEEP_ID}.log"

if [[ ! -f "${CWA_PATH}" ]]; then
  echo "ERROR: missing CWA file: ${CWA_PATH}" >&2
  exit 10
fi

if [[ -s "${FINAL_DIR}/predictions.csv" ]]; then
  echo "Already finished: ${ASLEEP_ID}"
  exit 0
fi

EXTRA_ARGS=()
if [[ -n "${LOCAL_REPO_PATH}" ]]; then
  EXTRA_ARGS+=(--local_repo_path "${LOCAL_REPO_PATH}")
fi
if [[ -n "${MODEL_WEIGHT_PATH}" ]]; then
  EXTRA_ARGS+=(--model_weight_path "${MODEL_WEIGHT_PATH}")
fi

echo "Running get_sleep for ${ASLEEP_ID}"
echo "Scratch work root: ${WORK_ROOT}"
echo "Final output: ${FINAL_DIR}"
echo "Log: ${RUN_LOG}"

"${GET_SLEEP_BIN}" "${CWA_PATH}" \
  --outdir "${WORK_ROOT}" \
  --pytorch_device "${DEVICE}" \
  --min_wear "${MIN_WEAR}" \
  --remove_intermediate_files \
  "${EXTRA_ARGS[@]}" \
  > "${RUN_LOG}" 2>&1

SCRATCH_DIR="${WORK_ROOT}/${CWA_STEM}"

required_files=(
  predictions.csv
  sleep_block.csv
  day_summary.csv
  summary.json
  info.json
)

for file_name in "${required_files[@]}"; do
  if [[ ! -s "${SCRATCH_DIR}/${file_name}" ]]; then
    echo "ERROR: required output missing or empty: ${SCRATCH_DIR}/${file_name}" >&2
    exit 20
  fi
done

PARTIAL_DIR="${OUT_ROOT}/${ASLEEP_ID}.partial.${SLURM_JOB_ID:-manual}.${SLURM_ARRAY_TASK_ID:-1}"
rm -rf "${PARTIAL_DIR}"
mkdir -p "${PARTIAL_DIR}"

for file_name in "${required_files[@]}"; do
  cp -f "${SCRATCH_DIR}/${file_name}" "${PARTIAL_DIR}/${file_name}"
done

if [[ -e "${FINAL_DIR}" && ! -d "${FINAL_DIR}" ]]; then
  echo "ERROR: final output path exists but is not a directory: ${FINAL_DIR}" >&2
  exit 30
fi

rm -rf "${FINAL_DIR}"
mv "${PARTIAL_DIR}" "${FINAL_DIR}"

echo "Finished: ${ASLEEP_ID}"
