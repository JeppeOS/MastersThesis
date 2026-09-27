#!/bin/bash

#SBATCH --job-name=flye_meta
#SBATCH --partition=normal
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/flye_%A_%a.out
#SBATCH --error=logs/flye_%A_%a.err
#SBATCH --account=methanotrophs

source ~/.bashrc

conda activate ont_metagenome

set -euo pipefail

## Get files and make outputs
PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome
MANIFEST="${PROJECT}/trimmed_fastqs.txt"
OUTDIR="${PROJECT}/flye_results"

mkdir -p "${OUTDIR}"

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")

if [[ -z "${FASTQ}" ]]; then
    echo "ERROR: No FASTQ found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BASENAME=$(basename "${FASTQ}")
SAMPLE=${BASENAME%_trimmed.fastq.gz}

echo "========================================"
echo "SLURM job:       ${SLURM_JOB_ID}"
echo "Array task:      ${SLURM_ARRAY_TASK_ID}"
echo "Sample:          ${SAMPLE}"
echo "Input FASTQ:     ${FASTQ}"
echo "Output:          ${OUTDIR}/${SAMPLE}"
echo "Threads:         ${SLURM_CPUS_PER_TASK}"
echo "========================================"

flye --nano-hq "${FASTQ}" \
    --meta \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --out-dir "${OUTDIR}/${SAMPLE}"

echo "Flye competed for ${SAMPLE}"