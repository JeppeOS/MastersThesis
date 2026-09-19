#!/usr/bin/env bash
#SBATCH --job-name=MMseqs2_exported_heme
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/06_clustering/logs/mmseqs2_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/06_clustering/logs/mmseqs2_%j.err

set -Eeo pipefail

PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
WORK_DIR="${PROJECT}/06_clustering"

INPUT_FASTA="${WORK_DIR}/exported_heme_candidates.faa"

RUN_NAME="mmseqs40_cov80"
RESULT_PREFIX="${WORK_DIR}/${RUN_NAME}"
TMP_DIR="${SLURM_TMPDIR:-${WORK_DIR}/tmp_${SLURM_JOB_ID}}"

MIN_SEQ_ID="0.40"
COVERAGE="0.80"
COV_MODE="0"

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

die() {
    echo "[$(timestamp)] ERROR: $*" >&2
    exit 1
}

mkdir -p "${WORK_DIR}/logs"

[[ -s "$INPUT_FASTA" ]] ||
    die "Input FASTA missing or empty: $INPUT_FASTA"

# Use an already available MMseqs2 if present. Otherwise try a dedicated
# conda environment named mmseqs2.
if ! command -v mmseqs >/dev/null 2>&1; then
    if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniforge3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    else
        die "mmseqs is not in PATH and conda.sh could not be found."
    fi

    conda activate mmseqs2 ||
        die "Could not activate conda environment 'mmseqs2'."
fi

command -v mmseqs >/dev/null 2>&1 ||
    die "mmseqs is not available."

echo "[$(timestamp)] MMseqs2 clustering starting."
echo "Executable:       $(command -v mmseqs)"
echo "Version:"
mmseqs version
echo "Input FASTA:      $INPUT_FASTA"
echo "Protein count:    $(grep -c '^>' "$INPUT_FASTA")"
echo "Minimum identity: $MIN_SEQ_ID"
echo "Coverage:         $COVERAGE"
echo "Coverage mode:    $COV_MODE (bidirectional)"
echo "Threads:          ${SLURM_CPUS_PER_TASK:-16}"
echo "Temporary dir:    $TMP_DIR"
echo

# Refuse to silently overwrite a completed clustering.
if [[ -s "${RESULT_PREFIX}_cluster.tsv" ]]; then
    die "Existing cluster TSV found: ${RESULT_PREFIX}_cluster.tsv"
fi

rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

export MMSEQS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

mmseqs easy-cluster \
    "$INPUT_FASTA" \
    "$RESULT_PREFIX" \
    "$TMP_DIR" \
    --min-seq-id "$MIN_SEQ_ID" \
    -c "$COVERAGE" \
    --cov-mode "$COV_MODE" \
    --threads "${SLURM_CPUS_PER_TASK:-16}"

[[ -s "${RESULT_PREFIX}_cluster.tsv" ]] ||
    die "MMseqs2 finished but cluster TSV is missing/empty."

[[ -s "${RESULT_PREFIX}_rep_seq.fasta" ]] ||
    die "MMseqs2 finished but representative FASTA is missing/empty."

n_input="$(grep -c '^>' "$INPUT_FASTA")"
n_memberships="$(wc -l < "${RESULT_PREFIX}_cluster.tsv")"
n_clusters="$(cut -f1 "${RESULT_PREFIX}_cluster.tsv" | sort -u | wc -l)"

if [[ "$n_memberships" -ne "$n_input" ]]; then
    die "Cluster membership rows ($n_memberships) != input proteins ($n_input)."
fi

cat > "${WORK_DIR}/${RUN_NAME}_RUN_COMPLETE.txt" <<EOF
MMseqs2 clustering completed successfully.
completion_time=$(timestamp)
job_id=${SLURM_JOB_ID}
input_fasta=${INPUT_FASTA}
input_proteins=${n_input}
clusters=${n_clusters}
min_seq_id=${MIN_SEQ_ID}
coverage=${COVERAGE}
cov_mode=${COV_MODE}
threads=${SLURM_CPUS_PER_TASK:-16}
mmseqs_executable=$(command -v mmseqs)
EOF

echo
echo "[$(timestamp)] MMseqs2 clustering complete."
echo "Input proteins: $n_input"
echo "Clusters:       $n_clusters"
echo "Cluster TSV:    ${RESULT_PREFIX}_cluster.tsv"
echo "Rep FASTA:      ${RESULT_PREFIX}_rep_seq.fasta"
