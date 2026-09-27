#!/bin/bash

#SBATCH --job-name=sb2_depth
#SBATCH --partition=normal
#SBATCH --array=1-7%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sb2_depth_%A_%a.out
#SBATCH --error=logs/sb2_depth_%A_%a.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc
conda activate ont_metagenome

set -euo pipefail


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

OUTDIR="${PROJECT}/semibin2_analysis/depth"

mkdir -p "${OUTDIR}"


## ------------------------------------------------------------
## Samples
## ------------------------------------------------------------

SAMPLES=(
    barcode10_seqs
    barcode11_seqs
    barcode12_seqs
    barcode13_seqs
    barcode14_seqs
    barcode15_seqs
    barcode16_seqs
)

SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"


ASSEMBLY="${PROJECT}/flye_results/${SAMPLE}/assembly.fasta"

BAM="${PROJECT}/mapping/${SAMPLE}/${SAMPLE}.sorted.bam"

DEPTH="${OUTDIR}/${SAMPLE}_depth.txt"


## ------------------------------------------------------------
## Validate
## ------------------------------------------------------------

if [[ ! -s "${ASSEMBLY}" ]]; then

    echo "ERROR: Assembly missing:"
    echo "${ASSEMBLY}"
    exit 1

fi


if [[ ! -s "${BAM}" ]]; then

    echo "ERROR: BAM missing:"
    echo "${BAM}"
    exit 1

fi


## ------------------------------------------------------------
## Generate depth table
##
## Same 97% mapping identity threshold used previously.
## ------------------------------------------------------------

jgi_summarize_bam_contig_depths \
    --percentIdentity 97 \
    --referenceFasta "${ASSEMBLY}" \
    --outputDepth "${DEPTH}" \
    "${BAM}"


echo "============================================================"
echo "Sample: ${SAMPLE}"
echo "Depth:  ${DEPTH}"
echo "============================================================"
