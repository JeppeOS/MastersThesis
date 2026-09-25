#!/usr/bin/env bash

#SBATCH --job-name=C121_local
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/comparative_analysis/Cluster_00121/local_clustering/logs/mmseqs_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/comparative_analysis/Cluster_00121/local_clustering/logs/mmseqs_%j.err

set -Eeo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

WORK="${PROJECT}/comparative_analysis/Cluster_00121"

OUT="${WORK}/local_clustering"

INPUT="${WORK}/Cluster_00121_combined_neighborhood_proteins.faa"

PREFIX="${OUT}/C121_local40_cov80"

TMP="${SLURM_TMPDIR:-${OUT}/tmp_${SLURM_JOB_ID}}"


mkdir -p \
    "$OUT" \
    "${OUT}/logs"


fail() {
    echo "ERROR: $*" >&2
    exit 1
}


[[ -s "$INPUT" ]] ||
    fail "Input FASTA missing: $INPUT"


command -v mmseqs >/dev/null 2>&1 ||
    fail "mmseqs is not available."


N_INPUT=$(grep -c '^>' "$INPUT")


echo "============================================================"
echo "STAGE 64 - CLUSTER_00121 SHARED LOCAL FAMILIES"
echo "============================================================"
echo
echo "Neighborhood proteins: $N_INPUT"
echo "Minimum identity:      40%"
echo "Minimum coverage:      80%"
echo "Coverage mode:         0"
echo


rm -rf "$TMP"

mkdir -p "$TMP"


mmseqs easy-cluster \
    "$INPUT" \
    "$PREFIX" \
    "$TMP" \
    --min-seq-id 0.40 \
    -c 0.80 \
    --cov-mode 0 \
    --threads "${SLURM_CPUS_PER_TASK:-8}"


CLUSTERS="${PREFIX}_cluster.tsv"


[[ -s "$CLUSTERS" ]] ||
    fail "MMseqs cluster table was not produced."


N_ASSIGNMENTS=$(wc -l < "$CLUSTERS")

N_FAMILIES=$(cut -f1 "$CLUSTERS" | sort -u | wc -l)


[[ "$N_ASSIGNMENTS" -eq "$N_INPUT" ]] ||
    fail \
        "Expected $N_INPUT assignments; found $N_ASSIGNMENTS"


echo
echo "============================================================"
echo "STAGE 64 COMPLETE"
echo "============================================================"
echo
echo "Protein assignments: $N_ASSIGNMENTS"
echo "Raw local families:  $N_FAMILIES"
echo
echo "Membership:"
echo "  $CLUSTERS"
