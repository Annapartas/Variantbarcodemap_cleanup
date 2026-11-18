#!/usr/bin/env bash
#SBATCH --job-name=BC_FC_filter
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=BC_FC_filter_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# Under SLURM, use the directory where sbatch was submitted.
# Outside SLURM, use this script's directory.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    cd "${SLURM_SUBMIT_DIR:?Submit from the repository directory}"
else
    cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
fi

# Select a parameter file with PARAMS (sbatch) or an argument (bash).
PARAM_FILE="${1:-${PARAMS:-barcoded_variant_filter_params.json}}"

# Optional Conda activation when CONDA_BASE is supplied.
if [[ -n "${CONDA_BASE:-}" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${BC_FILTER_ENV:-BC_FC_variant_filter}"
fi

if [[ ! -f BC_FC_variant_filter.py ]]; then
    echo "ERROR: BC_FC_variant_filter.py not found in $(pwd)" >&2
    exit 1
fi
if [[ ! -f "$PARAM_FILE" ]]; then
    echo "ERROR: parameter file not found: $PARAM_FILE" >&2
    exit 1
fi

PYTHON="${PYTHON:-python3}"
command -v "$PYTHON" >/dev/null || { echo "ERROR: $PYTHON not found" >&2; exit 1; }

echo "Job: ${SLURM_JOB_ID:-local}; directory: $(pwd)"
echo "Parameters: $PARAM_FILE"
"$PYTHON" BC_FC_variant_filter.py --parameters "$PARAM_FILE"
echo "BC/FC variant filtering completed successfully."
