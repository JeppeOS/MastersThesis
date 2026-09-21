#!/bin/bash

#SBATCH --job-name=checkm2_mags
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/checkm2_%j.out
#SBATCH --error=logs/checkm2_%j.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate checkm2_cpu

set -euo pipefail

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome
OUTDIR="${PROJECT}/checkm2_results"

mkdir -p "${OUTDIR}"

mapfile -t BINS < <(
    find "${PROJECT}/metabat2_results" \
        -regextype posix-extended \
        -type f \
        -regex '.*/[^/]+_bin\.[0-9]+\.fa' \
        | sort
)

echo "Number of bins: ${#BINS[@]}"

checkm2 predict \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --input "${BINS[@]}" \
    --output-directory "${OUTDIR}"

echo "CheckM2 completed successfully"