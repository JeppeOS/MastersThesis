#!/bin/bash

#SBATCH --job-name=das_checkm2
#SBATCH --partition=normal
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/das_checkm2_%j.out
#SBATCH --error=logs/das_checkm2_%j.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate checkm2_cpu

set -euo pipefail


PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

INPUT="${PROJECT}/dastool_evaluation/input_bins"

OUTDIR="${PROJECT}/dastool_evaluation/checkm2"


N=$(find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)


echo "DAS Tool bins: ${N}"


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
echo "============================================================"
echo "CheckM2 complete"
echo "Input bins: ${N}"
echo "Results:    ${N_RESULTS}"
echo "============================================================"
