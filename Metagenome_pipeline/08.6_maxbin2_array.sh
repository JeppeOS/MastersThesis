#!/bin/bash

#SBATCH --job-name=maxbin2
#SBATCH --partition=normal
#SBATCH --array=1-14%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/maxbin2_%A_%a.out
#SBATCH --error=logs/maxbin2_%A_%a.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## General setup
## ------------------------------------------------------------

source ~/.bashrc

set -euo pipefail


## ------------------------------------------------------------
## Project
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

OUTBASE="${PROJECT}/maxbin2_results"

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
## Inputs
## ------------------------------------------------------------

if [[ "${MODE}" == "flye" ]]; then

    ASSEMBLY="${PROJECT}/flye_results/${SAMPLE}/assembly.fasta"

    BAM="${PROJECT}/mapping/${SAMPLE}/${SAMPLE}.sorted.bam"

else

    ASSEMBLY="${PROJECT}/strainberry_results/${SAMPLE}/assembly.scaffolds.fa"

    BAM="${PROJECT}/strainberry_results/${SAMPLE}/assembly.scaffolds.bam"

fi


OUTDIR="${OUTBASE}/${MODE}/${SAMPLE}"

DEPTH="${OUTDIR}/${SAMPLE}_${MODE}_depth.tsv"

ABUNDANCE="${OUTDIR}/${SAMPLE}_${MODE}_abundance.tsv"

BINPREFIX="${OUTDIR}/${SAMPLE}_${MODE}_maxbin"


## ------------------------------------------------------------
## Validate
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
echo "MaxBin2"
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
## Generate abundance information
##
## Use the same MetaBAT2/JGI depth calculation and the same
## 97% identity requirement used in our previous binning.
## ------------------------------------------------------------

conda activate ont_metagenome


echo
echo "Generating contig depths..."


jgi_summarize_bam_contig_depths \
    --percentIdentity 97 \
    --referenceFasta "${ASSEMBLY}" \
    --outputDepth "${DEPTH}" \
    "${BAM}"


## ------------------------------------------------------------
## Convert JGI depth table into MaxBin abundance format
##
## JGI:
##
##   column 1 = contigName
##   column 2 = contigLen
##   column 3 = totalAvgDepth
##
## MaxBin:
##
##   contig<TAB>abundance
## ------------------------------------------------------------

awk -F'\t' '
    NR == 1 {
        if ($1 != "contigName" || $3 != "totalAvgDepth") {
            print "ERROR: Unexpected JGI depth header." > "/dev/stderr"
            exit 1
        }

        next
    }

    {
        print $1 "\t" $3
    }
' "${DEPTH}" > "${ABUNDANCE}"


if [[ ! -s "${ABUNDANCE}" ]]; then

    echo "ERROR: MaxBin2 abundance file was not generated."

    exit 1

fi


echo
echo "Abundance entries:"
wc -l "${ABUNDANCE}"


## ------------------------------------------------------------
## Run MaxBin2
## ------------------------------------------------------------

conda activate maxbin2


echo
echo "Running MaxBin2..."


run_MaxBin.pl \
    -contig "${ASSEMBLY}" \
    -abund "${ABUNDANCE}" \
    -out "${BINPREFIX}" \
    -thread "${SLURM_CPUS_PER_TASK}" \
    -min_contig_length 2500


## ------------------------------------------------------------
## Count MaxBin bins
##
## Typical outputs:
##
##   prefix.001.fasta
##   prefix.002.fasta
##   ...
## ------------------------------------------------------------

N_BINS=$(find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -regextype posix-extended \
    -regex ".*/$(basename "${BINPREFIX}")\.[0-9]+\.fasta" \
    | wc -l)


echo
echo "============================================================"
echo "MaxBin2 complete"
echo "Mode:   ${MODE}"
echo "Sample: ${SAMPLE}"
echo "Bins:   ${N_BINS}"
echo "============================================================"
