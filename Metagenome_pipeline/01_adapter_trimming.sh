#!/bin/bash

#SBATCH --job-name=porechop
#SBATCH --partition=normal
#SBATCH --array=1-8%5
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/porechop_%A_%a.out
#SBATCH --error=logs/porechop_%A_%a.err
#SBATCH --account=methanotrophs




source ~/.bashrc

conda activate ont_metagenome

set -euo pipefail

## Get files and make outputs
PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome
MANIFEST="${PROJECT}/metagenome_fastqs.txt"
OUTDIR="${PROJECT}/porechop_results"

mkdir -p "${OUTDIR}"

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")

if [[ -z "${FASTQ}" ]]; then
    echo "ERROR: No FASTQ found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BASENAME=$(basename "${FASTQ}")
SAMPLE=${BASENAME%.fastq.gz}

OUTPUT="${OUTDIR}/${SAMPLE}_trimmed.fastq.gz"

echo "========================================"
echo "Task:    ${SLURM_ARRAY_TASK_ID}"
echo "Sample:  ${SAMPLE}"
echo "Input:   ${FASTQ}"
echo "Output:  ${OUTPUT}"
echo "========================================"

porechop_abi \
    -i "${FASTQ}" \
    -o "${OUTPUT}" \
    --threads "${SLURM_CPUS_PER_TASK}"

echo "Porechop completed for ${SAMPLE}"

