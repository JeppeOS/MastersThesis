#!/bin/bash

#SBATCH --job-name=semibin2
#SBATCH --partition=normal
#SBATCH --array=1-14%3
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/semibin2_%A_%a.out
#SBATCH --error=logs/semibin2_%A_%a.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc
conda activate semibin2

set -euo pipefail


## ------------------------------------------------------------
## Project
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

OUTBASE="${PROJECT}/semibin2_results"

mkdir -p "${OUTBASE}/flye"
mkdir -p "${OUTBASE}/strainberry"


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


## ------------------------------------------------------------
## Determine sample and assembly type
##
## Tasks:
##
##   1-7   Flye
##   8-14  Strainberry
## ------------------------------------------------------------

TASK=${SLURM_ARRAY_TASK_ID}

if [[ "${TASK}" -le 7 ]]; then

    MODE="flye"

    IDX=$((TASK - 1))

else

    MODE="strainberry"

    IDX=$((TASK - 8))

fi


SAMPLE="${SAMPLES[$IDX]}"


## ------------------------------------------------------------
## Input paths
## ------------------------------------------------------------

if [[ "${MODE}" == "flye" ]]; then

    ASSEMBLY="${PROJECT}/flye_results/${SAMPLE}/assembly.fasta"

    BAM="${PROJECT}/mapping/${SAMPLE}/${SAMPLE}.sorted.bam"

else

    ASSEMBLY="${PROJECT}/strainberry_results/${SAMPLE}/assembly.scaffolds.fa"

    BAM="${PROJECT}/strainberry_results/${SAMPLE}/assembly.scaffolds.bam"

fi


OUTDIR="${OUTBASE}/${MODE}/${SAMPLE}"


## ------------------------------------------------------------
## Validate inputs
## ------------------------------------------------------------

if [[ ! -s "${ASSEMBLY}" ]]; then

    echo "ERROR: Missing assembly:"
    echo "${ASSEMBLY}"

    exit 1

fi


if [[ ! -s "${BAM}" ]]; then

    echo "ERROR: Missing BAM:"
    echo "${BAM}"

    exit 1

fi


## ------------------------------------------------------------
## Prevent overwrite
## ------------------------------------------------------------

if [[ -d "${OUTDIR}" ]] && \
   [[ -n "$(find "${OUTDIR}" -mindepth 1 -print -quit 2>/dev/null)" ]]
then

    echo "ERROR: Output directory already exists and is not empty:"
    echo "${OUTDIR}"

    exit 1

fi


mkdir -p "${OUTDIR}"


## ------------------------------------------------------------
## Summary
## ------------------------------------------------------------

echo "============================================================"
echo "SemiBin2"
echo
echo "Task:       ${TASK}"
echo "Mode:       ${MODE}"
echo "Sample:     ${SAMPLE}"
echo "Assembly:   ${ASSEMBLY}"
echo "BAM:        ${BAM}"
echo "Output:     ${OUTDIR}"
echo "Threads:    ${SLURM_CPUS_PER_TASK}"
echo "Min contig: 2500 bp"
echo "============================================================"


## ------------------------------------------------------------
## Run SemiBin2
##
## self-supervised:
##   train a sample-specific model rather than assuming a
##   predefined habitat.
##
## long_read:
##   appropriate mode for the ONT assemblies.
##
## min-len 2500:
##   same lower contig threshold used for MetaBAT2.
##
## random-seed 42:
##   reproducible clustering/training.
## ------------------------------------------------------------

SemiBin2 single_easy_bin \
    --self-supervised \
    --sequencing-type long_read \
    --input-fasta "${ASSEMBLY}" \
    --input-bam "${BAM}" \
    --output "${OUTDIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --min-len 2500 \
    --random-seed 42


## ------------------------------------------------------------
## Count output bins
## ------------------------------------------------------------

BINDIR="${OUTDIR}/output_bins"

if [[ ! -d "${BINDIR}" ]]; then

    echo
    echo "ERROR: SemiBin2 output_bins directory not found:"
    echo "${BINDIR}"

    exit 1

fi


N_BINS=$(find "${BINDIR}" \
    -maxdepth 1 \
    -type f \
    \( -name '*.fa' -o -name '*.fasta' -o -name '*.fa.gz' \) \
    | wc -l)


echo
echo "============================================================"
echo "SemiBin2 complete"
echo "Mode:   ${MODE}"
echo "Sample: ${SAMPLE}"
echo "Bins:   ${N_BINS}"
echo "============================================================"
