#!/usr/bin/env bash

#SBATCH --job-name=DeepTMHMM_SemiBin2

#SBATCH --partition=normal

#SBATCH --cpus-per-task=4

#SBATCH --mem=10G

#SBATCH --time=12:00:00

#SBATCH --array=1-187%30

#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/deeptmhmm/logs/deeptmhmm_%A_%a.out

#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/deeptmhmm/logs/deeptmhmm_%A_%a.err


set -Eeo pipefail


###############################################################################
## Paths
###############################################################################

PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"


WORK_DIR="${PROJECT}/functional_analysis/deeptmhmm"


INPUT_DIR="${WORK_DIR}/candidates_by_genome"


OUTPUT_DIR="${WORK_DIR}/results_by_genome"


LOG_DIR="${WORK_DIR}/logs"


## Reuse the exact same licensed DeepTMHMM installation.
DEEPTMHMM_DIR="/faststorage/project/methanotrophs/methanotroph_project/jeppe/DeepTMHMM-Academic-License-v1.0"


DEEPTMHMM="${DEEPTMHMM_DIR}/predict.py"


CONDA_ENV="DeepTMHMM"


###############################################################################
## Helpers
###############################################################################

timestamp(){

    date '+%Y-%m-%d %H:%M:%S'

}


die(){

    echo "[$(timestamp)] ERROR: $*" >&2

    exit 1

}


###############################################################################
## Initial checks
###############################################################################

mkdir -p \
    "$OUTPUT_DIR" \
    "$LOG_DIR"


[[ -d "$INPUT_DIR" ]] \
    || die \
    "Candidate FASTA directory not found: $INPUT_DIR"


[[ -s "$DEEPTMHMM" ]] \
    || die \
    "DeepTMHMM predict.py not found: $DEEPTMHMM"


###############################################################################
## Activate environment
###############################################################################

if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then

    source "$HOME/miniforge3/etc/profile.d/conda.sh"


elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then

    source "$HOME/miniconda3/etc/profile.d/conda.sh"


else

    die "Could not find conda.sh"

fi


conda activate "$CONDA_ENV"


###############################################################################
## Discover candidate FASTAs
###############################################################################

mapfile -t FASTAS < <(

    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.faa' \
        -print \
        | LC_ALL=C sort

)


N_FASTAS="${#FASTAS[@]}"


[[ "$N_FASTAS" -eq 187 ]] \
    || echo \
    "WARNING: expected 187 candidate FASTAs but detected $N_FASTAS."


###############################################################################
## Select array input
###############################################################################

idx=$((SLURM_ARRAY_TASK_ID - 1))


if (( idx < 0 || idx >= N_FASTAS )); then

    echo \
        "[$(timestamp)] No input assigned to this array task. Exiting."

    exit 0

fi


faa="${FASTAS[$idx]}"


genome="$(
    basename "$faa" .faa
)"


genome_out="${OUTPUT_DIR}/${genome}"


###############################################################################
## Report
###############################################################################

echo "[$(timestamp)] Starting DeepTMHMM candidate task."

echo "Job ID:             ${SLURM_JOB_ID:-NA}"

echo "Array task:         ${SLURM_ARRAY_TASK_ID:-NA}"

echo "Host:               $(hostname)"

echo "Genome:             $genome"

echo "Candidate FASTA:    $faa"

echo "Python:             $(command -v python3)"

echo "DeepTMHMM:          $DEEPTMHMM"


###############################################################################
## Handle genomes with zero candidates
###############################################################################

if [[ ! -s "$faa" ]]; then

    mkdir -p "$genome_out"


    printf \
        'No FeGenie/FindMeHemes candidates for %s\n' \
        "$genome" \
        > "${genome_out}/DEEPTMHMM_NO_CANDIDATES.txt"


    echo \
        "[$(timestamp)] No candidates. Nothing to run."


    exit 0

fi


###############################################################################
## Count candidates
###############################################################################

n_candidates="$(
    grep -c '^>' "$faa" \
    || true
)"


echo \
    "Candidate proteins: $n_candidates"


###############################################################################
## Resume-safe behaviour
###############################################################################

complete_file="${genome_out}/predicted_topologies.3line"


if [[ -s "$complete_file" ]]; then

    echo \
        "[$(timestamp)] Complete output already exists. Skipping."

    exit 0

fi


if [[ -d "$genome_out" ]]; then

    echo \
        "[$(timestamp)] Removing incomplete previous output: $genome_out"


    rm -rf "$genome_out"

fi


###############################################################################
## Run DeepTMHMM
###############################################################################

cd "$DEEPTMHMM_DIR"


export PYTHONUNBUFFERED=1


export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"


export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"


export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"


python3 "$DEEPTMHMM" \
    --fasta "$faa" \
    --output "$genome_out"


###############################################################################
## Validate output
###############################################################################

[[ -s "$complete_file" ]] \
    || die \
    "DeepTMHMM finished but predicted_topologies.3line is missing/empty."


###############################################################################
## Finish
###############################################################################

echo \
    "[$(timestamp)] DeepTMHMM completed successfully."


echo \
    "Genome:             $genome"


echo \
    "Candidates:         $n_candidates"


echo \
    "Output:             $genome_out"
