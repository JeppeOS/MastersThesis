#!/usr/bin/env bash

#SBATCH --job-name=MG_local_clusters
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/methylococcales_neighborhoods/local_clustering/logs/mmseqs_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/methylococcales_neighborhoods/local_clustering/logs/mmseqs_%j.err

set -Eeo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

WORK_DIR="${PROJECT}/functional_analysis/methylococcales_neighborhoods"

OUT_DIR="${WORK_DIR}/local_clustering"

INPUT="${WORK_DIR}/priority_neighborhood_proteins.faa"

PREFIX="${OUT_DIR}/MG_local40_cov80"

TMP="${SLURM_TMPDIR:-${OUT_DIR}/tmp_${SLURM_JOB_ID}}"


mkdir -p \
    "$OUT_DIR" \
    "${OUT_DIR}/logs"


timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


die() {
    echo "[$(timestamp)] ERROR: $*" >&2
    exit 1
}


[[ -s "$INPUT" ]] ||
    die "Missing or empty input FASTA: $INPUT"


## -------------------------------------------------------------- ##
## Activate MMseqs2 if it is not already available.
## -------------------------------------------------------------- ##

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


N_INPUT=$(grep -c '^>' "$INPUT")


if [[ "$N_INPUT" -ne 845 ]]; then
    die "Expected 845 input proteins but found $N_INPUT"
fi


rm -rf "$TMP"

mkdir -p "$TMP"


echo "============================================================"
echo "STAGE 57 - LOCAL NEIGHBORHOOD PROTEIN CLUSTERING"
echo "============================================================"
echo
echo "Input proteins:       $N_INPUT"
echo "Minimum identity:     40%"
echo "Minimum coverage:     80%"
echo "Coverage mode:        0"
echo


## -------------------------------------------------------------- ##
## Same sequence-family definition used in the previous local
## architecture workflow:
##
##     >=40% amino-acid identity
##     >=80% coverage
##
## --cov-mode 0 means coverage is evaluated using the shorter
## sequence denominator as used in the existing workflow.
## -------------------------------------------------------------- ##

mmseqs easy-cluster \
    "$INPUT" \
    "$PREFIX" \
    "$TMP" \
    --min-seq-id 0.40 \
    -c 0.80 \
    --cov-mode 0 \
    --threads "${SLURM_CPUS_PER_TASK:-8}"


CLUSTER_TSV="${PREFIX}_cluster.tsv"

[[ -s "$CLUSTER_TSV" ]] ||
    die "MMseqs cluster TSV was not created: $CLUSTER_TSV"


N_ASSIGNMENTS=$(wc -l < "$CLUSTER_TSV")

N_CLUSTERS=$(cut -f1 "$CLUSTER_TSV" | sort -u | wc -l)


if [[ "$N_ASSIGNMENTS" -ne "$N_INPUT" ]]; then
    die \
        "Expected one cluster assignment per protein. " \
        "Input=$N_INPUT assignments=$N_ASSIGNMENTS"
fi


echo
echo "============================================================"
echo "STAGE 57 COMPLETE"
echo "============================================================"
echo
echo "Input proteins:       $N_INPUT"
echo "Cluster assignments:  $N_ASSIGNMENTS"
echo "Local families:       $N_CLUSTERS"
echo
echo "Membership:"
echo "  $CLUSTER_TSV"
echo
echo "Representative FASTA:"
echo "  ${PREFIX}_rep_seq.fasta"
