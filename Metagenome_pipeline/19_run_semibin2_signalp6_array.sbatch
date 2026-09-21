#!/usr/bin/env bash
#SBATCH --job-name=SignalP6_SemiBin2
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-187%30
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/signalp/logs/signalp6_%A_%a.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/signalp/logs/signalp6_%A_%a.err

set -Eeo pipefail

###############################################################################
## User settings
###############################################################################

PROJECT_DIR="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

INPUT_DIR="${PROJECT_DIR}/functional_analysis/FeGenie_SemiBin2_full_run_20260918/fegenie_results/ORF_calls"

SIGNALP_WORK_DIR="${PROJECT_DIR}/functional_analysis/signalp"

## Reuse the exact SignalP 6 installation/models from the validated
## 631-genome workflow.
SIGNALP_PACKAGE_DIR="$HOME/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/02_signalp/signalP/signalp6_slow_sequential/signalp-6-package"

MODEL_DIR="${SIGNALP_PACKAGE_DIR}/models"

OUTPUT_ROOT="${SIGNALP_WORK_DIR}/results_by_genome"
LOG_DIR="${SIGNALP_WORK_DIR}/logs"

CONDA_ENV="signalp6"

SIGNALP_MODE="slow-sequential"

ORGANISM="other"

###############################################################################
## Helper functions
###############################################################################

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

die() {
    echo "[$(timestamp)] ERROR: $*" >&2
    exit 1
}

###############################################################################
## Initial checks
###############################################################################

echo "[$(timestamp)] Starting SignalP 6 array task."
echo "SLURM job ID:       ${SLURM_JOB_ID:-NA}"
echo "SLURM array task:   ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Host:               $(hostname)"
echo "Requested CPUs:     ${SLURM_CPUS_PER_TASK:-NA}"
echo "SignalP mode:       ${SIGNALP_MODE}"
echo

[[ -d "$INPUT_DIR" ]] ||
    die "Input directory does not exist: $INPUT_DIR"

[[ -d "$SIGNALP_PACKAGE_DIR" ]] ||
    die "SignalP package directory does not exist: $SIGNALP_PACKAGE_DIR"

[[ -d "$MODEL_DIR" ]] ||
    die "SignalP model directory does not exist: $MODEL_DIR"

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR"

###############################################################################
## Activate SignalP environment
###############################################################################

if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    die "Could not locate conda.sh under miniforge3 or miniconda3."
fi

conda activate "$CONDA_ENV"

command -v signalp6 >/dev/null 2>&1 ||
    die "signalp6 is not available after activating environment: $CONDA_ENV"

echo "SignalP executable: $(command -v signalp6)"
echo "Python executable:  $(command -v python)"
echo

###############################################################################
## Select one FeGenie proteome for this array task
###############################################################################

mapfile -t INPUT_FILES < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*-proteins.faa' \
        -size +0c \
        -print |
    LC_ALL=C sort
)

N_INPUTS="${#INPUT_FILES[@]}"

[[ "$N_INPUTS" -gt 0 ]] ||
    die "No non-empty *-proteins.faa files found in: $INPUT_DIR"

[[ "$N_INPUTS" -eq 187 ]] ||
    die "Expected 187 proteomes but detected ${N_INPUTS}."

echo "Detected proteomes: $N_INPUTS"

TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if (( TASK_INDEX < 0 || TASK_INDEX >= N_INPUTS )); then
    die "Array task ${SLURM_ARRAY_TASK_ID} has no corresponding input. Detected only ${N_INPUTS} files."
fi

INPUT_FILE="${INPUT_FILES[$TASK_INDEX]}"
INPUT_BASENAME="$(basename "$INPUT_FILE")"

GENOME_ID="${INPUT_BASENAME%-proteins.faa}"
GENOME_ID="${GENOME_ID%.fa}"
GENOME_ID="${GENOME_ID%.fna}"
GENOME_ID="${GENOME_ID%.fasta}"

FINAL_OUT="${OUTPUT_ROOT}/${GENOME_ID}"
COMPLETE_MARKER="${FINAL_OUT}/SIGNALP6_COMPLETE.txt"

echo "Selected input:      $INPUT_FILE"
echo "Genome ID:           $GENOME_ID"
echo "Final output:        $FINAL_OUT"
echo

###############################################################################
## Resume-safe behaviour
###############################################################################

if [[ -s "$COMPLETE_MARKER" ]]; then
    echo "[$(timestamp)] Complete output already exists. Skipping:"
    echo "$COMPLETE_MARKER"
    exit 0
fi

TMP_PARENT="${SLURM_TMPDIR:-${SIGNALP_WORK_DIR}/tmp}"
TMP_OUT="${TMP_PARENT}/signalp6_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${GENOME_ID}"

rm -rf "$TMP_OUT"
mkdir -p "$TMP_OUT"

cleanup() {
    status=$?
    if [[ "$status" -ne 0 ]]; then
        echo "[$(timestamp)] Task failed with exit status $status." >&2
        echo "Temporary output retained at: $TMP_OUT" >&2
    fi
}
trap cleanup EXIT

###############################################################################
## Run SignalP 6
###############################################################################

echo "[$(timestamp)] Running SignalP 6."

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

signalp6 \
    --fastafile "$INPUT_FILE" \
    --organism "$ORGANISM" \
    --output_dir "$TMP_OUT" \
    --format none \
    --mode "$SIGNALP_MODE" \
    --model_dir "$MODEL_DIR"

###############################################################################
## Validate and publish output
###############################################################################

PREDICTION_FILE="${TMP_OUT}/prediction_results.txt"

[[ -s "$PREDICTION_FILE" ]] ||
    die "SignalP finished without creating a non-empty prediction_results.txt"

N_INPUT_PROTEINS="$(
    grep -c '^>' "$INPUT_FILE" || true
)"

N_PREDICTIONS="$(
    awk '
        BEGIN { n = 0 }
        /^[[:space:]]*#/ { next }
        NF > 0 { n++ }
        END { print n }
    ' "$PREDICTION_FILE"
)"

if [[ "$N_PREDICTIONS" -ne "$N_INPUT_PROTEINS" ]]; then
    die "Prediction-row count (${N_PREDICTIONS}) does not match input-protein count (${N_INPUT_PROTEINS})."
fi

rm -rf "$FINAL_OUT"
mkdir -p "$(dirname "$FINAL_OUT")"
mv "$TMP_OUT" "$FINAL_OUT"

cat > "$COMPLETE_MARKER" <<EOF
SignalP 6 task completed successfully.
completion_time=$(timestamp)
job_id=${SLURM_JOB_ID}
array_task_id=${SLURM_ARRAY_TASK_ID}
host=$(hostname)
input_file=${INPUT_FILE}
genome_id=${GENOME_ID}
organism=${ORGANISM}
mode=${SIGNALP_MODE}
model_dir=${MODEL_DIR}
input_proteins=${N_INPUT_PROTEINS}
prediction_rows=${N_PREDICTIONS}
EOF

trap - EXIT

echo
echo "[$(timestamp)] SignalP 6 completed successfully."
echo "Proteins processed:  $N_INPUT_PROTEINS"
echo "Prediction rows:     $N_PREDICTIONS"
echo "Output directory:    $FINAL_OUT"
