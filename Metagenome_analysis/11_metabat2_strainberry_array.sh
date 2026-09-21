#!/bin/bash

#SBATCH --job-name=sb_metabat
#SBATCH --partition=normal
#SBATCH --array=1-7%3
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/sb_metabat_%A_%a.out
#SBATCH --error=logs/sb_metabat_%A_%a.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc

## MetaBAT2 + samtools were already available here ##
conda activate ont_metagenome

set -euo pipefail


## ------------------------------------------------------------
## Project paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

STRAINDIR="${PROJECT}/strainberry_results"
OUTBASE="${PROJECT}/strainberry_metabat2"

mkdir -p "${OUTBASE}"


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


## ------------------------------------------------------------
## Inputs
## ------------------------------------------------------------

INDIR="${STRAINDIR}/${SAMPLE}"

ASSEMBLY="${INDIR}/assembly.scaffolds.fa"
BAM="${INDIR}/assembly.scaffolds.bam"

OUTDIR="${OUTBASE}/${SAMPLE}"

DEPTH="${OUTDIR}/${SAMPLE}_depth.txt"

BINPREFIX="${OUTDIR}/${SAMPLE}_bin"


## ------------------------------------------------------------
## Create output directory
## ------------------------------------------------------------

mkdir -p "${OUTDIR}"


## ------------------------------------------------------------
## Validate inputs
## ------------------------------------------------------------

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: Missing Strainberry assembly:"
    echo "${ASSEMBLY}"
    exit 1
fi

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: Missing Strainberry BAM:"
    echo "${BAM}"
    exit 1
fi


samtools quickcheck -v "${BAM}"


## ------------------------------------------------------------
## Summary
## ------------------------------------------------------------

echo "============================================================"
echo "Sample:       ${SAMPLE}"
echo "Assembly:     ${ASSEMBLY}"
echo "BAM:          ${BAM}"
echo "Output:       ${OUTDIR}"
echo "============================================================"


## ------------------------------------------------------------
## Generate MetaBAT2 depth file
##
## Same 97% alignment identity cutoff used for the original
## Flye/MetaBAT2 workflow.
## ------------------------------------------------------------

echo
echo "Generating contig depth table..."

jgi_summarize_bam_contig_depths \
    --percentIdentity 97 \
    --referenceFasta "${ASSEMBLY}" \
    --outputDepth "${DEPTH}" \
    "${BAM}"


## ------------------------------------------------------------
## Run MetaBAT2
##
## Same settings as original workflow:
##
##   minimum contig length = 2500 bp
##   seed                  = 42
## ------------------------------------------------------------

echo
echo "Running MetaBAT2..."

metabat2 \
    -i "${ASSEMBLY}" \
    -a "${DEPTH}" \
    -o "${BINPREFIX}" \
    -m 2500 \
    -t "${SLURM_CPUS_PER_TASK}" \
    --seed 42 \
    --saveCls \
    --unbinned


## ------------------------------------------------------------
## Count actual numeric bins
## ------------------------------------------------------------

N_BINS=$(find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -regextype posix-extended \
    -regex ".*/${SAMPLE}_bin\.[0-9]+\.fa" \
    | wc -l)


echo
echo "============================================================"
echo "MetaBAT2 completed"
echo "Sample: ${SAMPLE}"
echo "Numeric bins: ${N_BINS}"
echo "============================================================"
