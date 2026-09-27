#!/bin/bash

#SBATCH --job-name=sb_checkm2
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sb_checkm2_%j.out
#SBATCH --error=logs/sb_checkm2_%j.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc
conda activate checkm2_cpu

set -euo pipefail


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

BINROOT="${PROJECT}/strainberry_metabat2"

OUTDIR="${PROJECT}/strainberry_checkm2_results"

INPUTDIR="${OUTDIR}/input_bins"


## ------------------------------------------------------------
## Prepare clean input directory
## ------------------------------------------------------------

rm -rf "${OUTDIR}"

mkdir -p "${INPUTDIR}"


## ------------------------------------------------------------
## Link numbered MetaBAT2 bins only
## ------------------------------------------------------------

find "${BINROOT}" \
    -type f \
    -regextype posix-extended \
    -regex '.*/barcode(10|11|12|13|14|15|16)_seqs_bin\.[0-9]+\.fa' \
    -print0 \
| while IFS= read -r -d '' f
do

    ln -s "${f}" \
        "${INPUTDIR}/$(basename "${f}")"

done


## ------------------------------------------------------------
## Validate input count
## ------------------------------------------------------------

N_BINS=$(find "${INPUTDIR}" \
    -maxdepth 1 \
    -name '*.fa' \
    | wc -l)


echo "============================================================"
echo "Strainberry CheckM2"
echo "Input directory: ${INPUTDIR}"
echo "Output:          ${OUTDIR}"
echo "Bins:            ${N_BINS}"
echo "CPUs:            ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


if [[ "${N_BINS}" -ne 214 ]]; then

    echo "ERROR: Expected 214 numbered bins but found ${N_BINS}."
    exit 1

fi


## ------------------------------------------------------------
## Run CheckM2
##
## IMPORTANT:
## -x fa tells CheckM2 that genome files end in .fa
## ------------------------------------------------------------

checkm2 predict \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --input "${INPUTDIR}" \
    --output-directory "${OUTDIR}/checkm2" \
    -x fa


## ------------------------------------------------------------
## Validate result
## ------------------------------------------------------------

REPORT="${OUTDIR}/checkm2/quality_report.tsv"

if [[ ! -s "${REPORT}" ]]; then

    echo "ERROR: CheckM2 quality_report.tsv not found."
    exit 1

fi


N_RESULTS=$(awk '
    NR > 1 {n++}
    END {print n+0}
' "${REPORT}")


echo
echo "============================================================"
echo "CheckM2 complete"
echo "Expected bins: ${N_BINS}"
echo "Results:       ${N_RESULTS}"
echo "Report:        ${REPORT}"
echo "============================================================"
