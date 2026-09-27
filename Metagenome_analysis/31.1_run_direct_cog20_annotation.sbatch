#!/usr/bin/env bash

#SBATCH --job-name=MG_COG20
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/methylococcales_neighborhoods/direct_cog20/logs/cog20_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/methylococcales_neighborhoods/direct_cog20/logs/cog20_%j.err

set -Eeo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

WORK_DIR="${PROJECT}/functional_analysis/methylococcales_neighborhoods"

OUT_DIR="${WORK_DIR}/direct_cog20"

INPUT="${WORK_DIR}/priority_neighborhood_proteins.faa"

OUTPUT="${OUT_DIR}/priority_neighborhood_proteins_COG20.tsv"

COG_DIR="${PROJECT}/reference_data/anvio_cogs"

EXPECTED=845


mkdir -p "$OUT_DIR"


## Activate isolated anvi'o environment. ##

source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate anvio-9


command -v anvi-run-ncbi-cogs >/dev/null 2>&1 || {
    echo "ERROR: anvi-run-ncbi-cogs not found."
    exit 1
}


[[ -s "$INPUT" ]] || {
    echo "ERROR: input FASTA missing."
    exit 1
}


N_INPUT=$(grep -c '^>' "$INPUT")


[[ "$N_INPUT" -eq "$EXPECTED" ]] || {
    echo "ERROR: expected $EXPECTED proteins; found $N_INPUT."
    exit 1
}


echo "============================================================"
echo "STAGE 59 - DIRECT COG20 ANNOTATION"
echo "============================================================"
echo
echo "Proteins:    $N_INPUT"
echo "COG version: COG20"
echo "Threads:     ${SLURM_CPUS_PER_TASK:-16}"
echo


anvi-run-ncbi-cogs \
    --fasta-file "$INPUT" \
    --output-file "$OUTPUT" \
    --cog-version COG20 \
    --cog-data-dir "$COG_DIR" \
    --search-with diamond \
    --num-threads "${SLURM_CPUS_PER_TASK:-16}"


[[ -s "$OUTPUT" ]] || {
    echo "ERROR: COG20 output was not created."
    exit 1
}


echo
echo "============================================================"
echo "STAGE 59 COMPLETE"
echo "============================================================"
echo
echo "Output:"
echo "  $OUTPUT"
echo
echo "Lines:"
wc -l "$OUTPUT"
echo
echo "Header / first records:"
head -n 10 "$OUTPUT"
