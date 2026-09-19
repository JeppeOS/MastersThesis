#!/usr/bin/env bash
#SBATCH --job-name=DeepTMHMM_candidates
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --array=1-631%30
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/04_deeptmhmm/logs/deeptmhmm_%A_%a.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/04_deeptmhmm/logs/deeptmhmm_%A_%a.err

set -Eeo pipefail

PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
WORK_DIR="${PROJECT}/04_deeptmhmm"
INPUT_DIR="${WORK_DIR}/candidates_by_genome"
OUTPUT_DIR="${WORK_DIR}/results_by_genome"
LOG_DIR="${WORK_DIR}/logs"
DEEPTMHMM_DIR="/faststorage/project/methanotrophs/methanotroph_project/jeppe/DeepTMHMM-Academic-License-v1.0"
DEEPTMHMM="${DEEPTMHMM_DIR}/predict.py"
CONDA_ENV="DeepTMHMM"

timestamp(){ date '+%Y-%m-%d %H:%M:%S'; }
die(){ echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
[[ -d "$INPUT_DIR" ]] || die "Candidate FASTA directory not found: $INPUT_DIR"
[[ -s "$DEEPTMHMM" ]] || die "DeepTMHMM predict.py not found: $DEEPTMHMM"

if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    die "Could not find conda.sh"
fi
conda activate "$CONDA_ENV"

mapfile -t FASTAS < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.faa' -print | LC_ALL=C sort)
N_FASTAS="${#FASTAS[@]}"
[[ "$N_FASTAS" -eq 631 ]] || echo "WARNING: expected 631 candidate FASTAs but detected $N_FASTAS."
idx=$((SLURM_ARRAY_TASK_ID - 1))
if (( idx < 0 || idx >= N_FASTAS )); then
    echo "[$(timestamp)] No input assigned to this array task. Exiting."
    exit 0
fi

faa="${FASTAS[$idx]}"
genome="$(basename "$faa" .faa)"
genome_out="${OUTPUT_DIR}/${genome}"

echo "[$(timestamp)] Starting DeepTMHMM candidate task."
echo "Job ID:             ${SLURM_JOB_ID:-NA}"
echo "Array task:         ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Host:               $(hostname)"
echo "Genome:             $genome"
echo "Candidate FASTA:    $faa"
echo "Python:             $(command -v python3)"
echo "DeepTMHMM:          $DEEPTMHMM"

if [[ ! -s "$faa" ]]; then
    mkdir -p "$genome_out"
    printf 'No FeGenie/FindMeHemes candidates for %s\n' "$genome" > "${genome_out}/DEEPTMHMM_NO_CANDIDATES.txt"
    echo "[$(timestamp)] No candidates. Nothing to run."
    exit 0
fi

n_candidates="$(grep -c '^>' "$faa" || true)"
echo "Candidate proteins:  $n_candidates"

complete_file="${genome_out}/predicted_topologies.3line"
if [[ -s "$complete_file" ]]; then
    echo "[$(timestamp)] Complete output already exists. Skipping."
    exit 0
fi
if [[ -d "$genome_out" ]]; then
    echo "[$(timestamp)] Removing incomplete previous output: $genome_out"
    rm -rf "$genome_out"
fi

cd "$DEEPTMHMM_DIR"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

python3 "$DEEPTMHMM" --fasta "$faa" --output "$genome_out"
[[ -s "$complete_file" ]] || die "DeepTMHMM finished but predicted_topologies.3line is missing/empty."

echo "[$(timestamp)] DeepTMHMM completed successfully."
echo "Genome:             $genome"
echo "Candidates:         $n_candidates"
echo "Output:             $genome_out"
