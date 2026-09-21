#!/bin/bash

#SBATCH --job-name=sb_gtdbtk
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sb_gtdbtk_%j.out
#SBATCH --error=logs/sb_gtdbtk_%j.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc
conda activate gtdbtk_2.6.1

set -euo pipefail


## ------------------------------------------------------------
## Project paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

BINROOT="${PROJECT}/strainberry_metabat2"

OUTDIR="${PROJECT}/strainberry_gtdbtk_results"

BATCHFILE="${PROJECT}/strainberry_gtdbtk_bins.tsv"


## ------------------------------------------------------------
## GTDB-Tk database
##
## Same GTDB release used for the original MAG analysis.
## ------------------------------------------------------------

export GTDBTK_DATA_PATH=\
/home/jeppeos/methanotrophs/methanotroph_project/jeppe/reference_data/gtdbtk_r226


## ------------------------------------------------------------
## Prevent accidental overwrite
## ------------------------------------------------------------

if [[ -d "${OUTDIR}" ]] && \
   [[ -n "$(find "${OUTDIR}" -mindepth 1 -print -quit 2>/dev/null)" ]]
then

    echo "ERROR: GTDB-Tk output directory already exists and is not empty:"
    echo "${OUTDIR}"
    echo
    echo "Remove or rename it before rerunning."
    exit 1

fi

mkdir -p "${OUTDIR}"


## ------------------------------------------------------------
## Build batch file
##
## Only numbered MetaBAT2 bins are included.
##
## Excluded automatically:
##   unbinned FASTAs
##   too-short FASTAs
##   auxiliary MetaBAT files
##
## GTDB batch-file format:
##
##   FASTA_path<TAB>genome_ID
## ------------------------------------------------------------

: > "${BATCHFILE}"


find "${BINROOT}" \
    -type f \
    -regextype posix-extended \
    -regex '.*/barcode(10|11|12|13|14|15|16)_seqs_bin\.[0-9]+\.fa' \
    | sort -V \
    | while read -r fasta
do

    genome_id=$(basename "${fasta}" .fa)

    printf "%s\t%s\n" \
        "${fasta}" \
        "${genome_id}" \
        >> "${BATCHFILE}"

done


## ------------------------------------------------------------
## Validate input
## ------------------------------------------------------------

N_BINS=$(wc -l < "${BATCHFILE}")


echo "============================================================"
echo "Strainberry GTDB-Tk"
echo
echo "GTDB release:    R226"
echo "Database:        ${GTDBTK_DATA_PATH}"
echo "Bin root:        ${BINROOT}"
echo "Batch file:      ${BATCHFILE}"
echo "Output:          ${OUTDIR}"
echo "Bins:            ${N_BINS}"
echo "CPUs:            ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


if [[ "${N_BINS}" -ne 214 ]]; then

    echo
    echo "ERROR: Expected 214 numbered bins but found ${N_BINS}."
    exit 1

fi


## ------------------------------------------------------------
## Show first and last few genomes as a sanity check
## ------------------------------------------------------------

echo
echo "First five batch entries:"
head -n 5 "${BATCHFILE}"

echo
echo "Last five batch entries:"
tail -n 5 "${BATCHFILE}"


## ------------------------------------------------------------
## Run GTDB-Tk
## ------------------------------------------------------------

echo
echo "Starting GTDB-Tk classify_wf..."
echo


gtdbtk classify_wf \
    --batchfile "${BATCHFILE}" \
    --out_dir "${OUTDIR}" \
    --cpus "${SLURM_CPUS_PER_TASK}"


## ------------------------------------------------------------
## Locate summary files
## ------------------------------------------------------------

BAC="${OUTDIR}/gtdbtk.bac120.summary.tsv"
ARC="${OUTDIR}/gtdbtk.ar53.summary.tsv"


echo
echo "============================================================"
echo "GTDB-Tk finished"
echo "============================================================"


## ------------------------------------------------------------
## Count bacterial classifications
## ------------------------------------------------------------

N_BAC=0

if [[ -s "${BAC}" ]]; then

    N_BAC=$(awk '
        NR > 1 {n++}
        END {print n+0}
    ' "${BAC}")

    echo "Bacterial genomes: ${N_BAC}"
    echo "Bacterial summary: ${BAC}"

else

    echo "No bacterial summary found."

fi


## ------------------------------------------------------------
## Count archaeal classifications
## ------------------------------------------------------------

N_ARC=0

if [[ -s "${ARC}" ]]; then

    N_ARC=$(awk '
        NR > 1 {n++}
        END {print n+0}
    ' "${ARC}")

    echo "Archaeal genomes:  ${N_ARC}"
    echo "Archaeal summary:  ${ARC}"

else

    echo "No archaeal summary found."

fi


## ------------------------------------------------------------
## Total classified / processed
## ------------------------------------------------------------

N_TOTAL=$((N_BAC + N_ARC))

echo
echo "Total in bacterial + archaeal summaries: ${N_TOTAL}"
echo "Expected input bins:                     ${N_BINS}"


## ------------------------------------------------------------
## Extract Methylococcales classifications for quick inspection
## ------------------------------------------------------------

METHYLO="${OUTDIR}/Methylococcales_GTDB.tsv"

if [[ -s "${BAC}" ]]; then

    awk -F'\t' '
        NR == 1 ||
        $0 ~ /o__Methylococcales/
    ' "${BAC}" > "${METHYLO}"

    N_METHYLO=$(awk '
        NR > 1 {n++}
        END {print n+0}
    ' "${METHYLO}")

    echo
    echo "Methylococcales MAGs: ${N_METHYLO}"
    echo "Methylococcales table:"
    echo "${METHYLO}"

fi


echo
echo "============================================================"
echo "Done"
echo "============================================================"
