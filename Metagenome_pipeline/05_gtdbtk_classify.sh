#!/bin/bash

#SBATCH --job-name=gtdbtk_mags
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gtdbtk_%j.out
#SBATCH --error=logs/gtdbtk_%j.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate gtdbtk_2.6.1

set -euo pipefail

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome
BATCH="${PROJECT}/gtdbtk_bins.tsv"
OUTDIR="${PROJECT}/gtdbtk_results"

mkdir -p "${OUTDIR}"

echo "GTDB-Tk version:"
gtdbtk --version

echo "Reference database:"
echo "${GTDBTK_DATA_PATH}"

echo "Number of bins:"
wc -l "${BATCH}"

gtdbtk classify_wf \
    --batchfile "${BATCH}" \
    --out_dir "${OUTDIR}" \
    --cpus "${SLURM_CPUS_PER_TASK}"

echo "GTDB-Tk completed successfully."