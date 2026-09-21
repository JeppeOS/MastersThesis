#!/bin/bash

#SBATCH --job-name=map_meta
#SBATCH --partition=normal
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/map_%A_%a.out
#SBATCH --error=logs/map_%A_%a.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate ont_metagenome

set -euo pipefail

## Paths ##
PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome
MANIFEST="${PROJECT}/trimmed_fastqs.txt"
FLYEDIR="${PROJECT}/flye_results"
OUTDIR="${PROJECT}/mapping"

mkdir -p "${OUTDIR}"

## Get FASTQ for this task ##
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")

if [[ -z "${FASTQ}" ]]; then
    echo "ERROR: No FASTQ found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

## Sample name ##
BASENAME=$(basename "${FASTQ}")
SAMPLE=${BASENAME%_trimmed.fastq.gz}

## Corresponding Flye assembly ##
ASSEMBLY="${FLYEDIR}/${SAMPLE}/assembly.fasta"

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found:"
    echo "${ASSEMBLY}"
    exit 1
fi

SAMPLEDIR="${OUTDIR}/${SAMPLE}"
mkdir -p "${SAMPLEDIR}"

BAM="${SAMPLEDIR}/${SAMPLE}.sorted.bam"

echo "========================================"
echo "Task:       ${SLURM_ARRAY_TASK_ID}"
echo "Sample:     ${SAMPLE}"
echo "Reads:      ${FASTQ}"
echo "Assembly:   ${ASSEMBLY}"
echo "Output BAM: ${BAM}"
echo "Threads:    ${SLURM_CPUS_PER_TASK}"
echo "========================================"

## Map reads back to assembly and coordinate-sort ##
minimap2 \
    -ax map-ont \
    -t "${SLURM_CPUS_PER_TASK}" \
    "${ASSEMBLY}" \
    "${FASTQ}" \
    | samtools sort \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -o "${BAM}" \
        -

## Index BAM ##
samtools index \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${BAM}"

## Basic mapping statistics ##
samtools flagstat \
    -@ "${SLURM_CPUS_PER_TASK}" \
    "${BAM}" \
    > "${SAMPLEDIR}/${SAMPLE}.flagstat.txt"

echo "Mapping completed for ${SAMPLE}"