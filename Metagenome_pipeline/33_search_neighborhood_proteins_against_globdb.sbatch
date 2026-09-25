#!/usr/bin/env bash

#SBATCH --job-name=NBH_GlobDB_annot
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/methylococcales_neighborhoods/globdb_annotation_transfer/logs/mmseqs_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/methylococcales_neighborhoods/globdb_annotation_transfer/logs/mmseqs_%j.err

set -Eeo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

WORK_DIR="${PROJECT}/functional_analysis/methylococcales_neighborhoods"

TRANSFER_DIR="${WORK_DIR}/globdb_annotation_transfer"


QUERY_FASTA="${WORK_DIR}/priority_neighborhood_proteins.faa"

TARGET_FASTA="${TRANSFER_DIR}/GlobDB631_all_proteins.faa"

RESULTS="${TRANSFER_DIR}/neighborhood_vs_GlobDB631.tsv"

TMP_DIR="${SLURM_TMPDIR:-${TRANSFER_DIR}/tmp_${SLURM_JOB_ID}}"


timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


die() {
    echo "[$(timestamp)] ERROR: $*" >&2
    exit 1
}


for path in \
    "$QUERY_FASTA" \
    "$TARGET_FASTA"
do
    [[ -s "$path" ]] ||
        die "Missing or empty input: $path"
done


if ! command -v mmseqs >/dev/null 2>&1; then

    if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniforge3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    else
        die "Cannot locate conda.sh"
    fi

    conda activate MMseqs2 ||
        die "Could not activate mmseqs2"
fi


rm -rf "$TMP_DIR"

mkdir -p "$TMP_DIR"


echo "[$(timestamp)] Searching neighborhood proteins against GlobDB631."
echo "Queries: $(grep -c '^>' "$QUERY_FASTA")"
echo "Targets: $(grep -c '^>' "$TARGET_FASTA")"


mmseqs easy-search \
    "$QUERY_FASTA" \
    "$TARGET_FASTA" \
    "$RESULTS" \
    "$TMP_DIR" \
    -s 7.5 \
    -e 1e-5 \
    --max-seqs 50 \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --format-output \
"query,target,pident,alnlen,qcov,tcov,evalue,bits"


[[ -f "$RESULTS" ]] ||
    die "MMseqs did not produce results."


echo
echo "[$(timestamp)] Search complete."
echo "Result rows: $(wc -l < "$RESULTS")"
echo "Results: $RESULTS"
