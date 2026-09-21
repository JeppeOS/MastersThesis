#!/bin/bash

#SBATCH --job-name=bin_checkm2
#SBATCH --partition=normal
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bin_checkm2_%j.out
#SBATCH --error=logs/bin_checkm2_%j.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate checkm2_cpu

set -euo pipefail

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

INPUT="${PROJECT}/binner_comparison/input_bins"

OUTDIR="${PROJECT}/binner_comparison/checkm2"


N=$(find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)


echo "CheckM2 input bins: ${N}"

if [[ "${N}" -ne 612 ]]; then
    echo "ERROR: Expected 612 bins."
    exit 1
fi


rm -rf "${OUTDIR}"


checkm2 predict \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --input "${INPUT}" \
    --output-directory "${OUTDIR}" \
    -x fa


REPORT="${OUTDIR}/quality_report.tsv"

N_RESULTS=$(awk '
    NR > 1 {n++}
    END {print n+0}
' "${REPORT}")


echo
echo "CheckM2 complete"
echo "Results: ${N_RESULTS}"

if [[ "${N_RESULTS}" -ne 612 ]]; then
    echo "WARNING: Expected 612 results."
fi
