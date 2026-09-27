#!/bin/bash

#SBATCH --job-name=dastool_flye
#SBATCH --partition=normal
#SBATCH --array=1-7%3
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/dastool_%A_%a.out
#SBATCH --error=logs/dastool_%A_%a.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc
conda activate dastool

set -euo pipefail


## ------------------------------------------------------------
## Project
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

OUTBASE="${PROJECT}/dastool_results/flye"

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
## Input assembly
##
## All three binners below used this SAME Flye assembly.
## ------------------------------------------------------------

ASSEMBLY="${PROJECT}/flye_results/${SAMPLE}/assembly.fasta"


## ------------------------------------------------------------
## Binner directories
## ------------------------------------------------------------

METABAT_DIR="${PROJECT}/metabat2_results/${SAMPLE}/bins"

SEMIBIN_DIR="${PROJECT}/semibin2_results/flye/${SAMPLE}/output_bins"

MAXBIN_DIR="${PROJECT}/maxbin2_results/flye/${SAMPLE}"


## ------------------------------------------------------------
## Output
## ------------------------------------------------------------

OUTDIR="${OUTBASE}/${SAMPLE}"

PREFIX="${OUTDIR}/${SAMPLE}"

METABAT_C2B="${OUTDIR}/metabat2_contigs2bin.tsv"

SEMIBIN_C2B="${OUTDIR}/semibin2_contigs2bin.tsv"

MAXBIN_C2B="${OUTDIR}/maxbin2_contigs2bin.tsv"


## ------------------------------------------------------------
## Validate inputs
## ------------------------------------------------------------

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: Flye assembly missing:"
    echo "${ASSEMBLY}"
    exit 1
fi

for dir in \
    "${METABAT_DIR}" \
    "${SEMIBIN_DIR}" \
    "${MAXBIN_DIR}"
do

    if [[ ! -d "${dir}" ]]; then
        echo "ERROR: Missing bin directory:"
        echo "${dir}"
        exit 1
    fi

done


## ------------------------------------------------------------
## Prevent accidental overwrite
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
## Function: convert one bin FASTA into contig2bin entries
##
## $1 = FASTA
## $2 = bin ID
## $3 = output TSV
## ------------------------------------------------------------

add_fasta_to_c2b () {

    fasta="$1"
    bin_id="$2"
    outfile="$3"

    if [[ "${fasta}" == *.gz ]]; then

        gzip -cd "${fasta}" |

        awk \
            -v BIN="${bin_id}" \
            '
            /^>/ {
                sub(/^>/, "", $0)
                split($0, a, /[[:space:]]+/)
                print a[1] "\t" BIN
            }
            ' >> "${outfile}"

    else

        awk \
            -v BIN="${bin_id}" \
            '
            /^>/ {
                sub(/^>/, "", $0)
                split($0, a, /[[:space:]]+/)
                print a[1] "\t" BIN
            }
            ' "${fasta}" >> "${outfile}"

    fi
}


## ------------------------------------------------------------
## MetaBAT2 → contigs2bin
##
## Only numbered bins are used.
## ------------------------------------------------------------

: > "${METABAT_C2B}"

while IFS= read -r -d '' fasta
do

    base=$(basename "${fasta}" .fa)

    bin_id="metabat2.${base}"

    add_fasta_to_c2b \
        "${fasta}" \
        "${bin_id}" \
        "${METABAT_C2B}"

done < <(

    find "${METABAT_DIR}" \
        -maxdepth 1 \
        -type f \
        -regextype posix-extended \
        -regex ".*/${SAMPLE}_bin\.[0-9]+\.fa" \
        -print0
)


## ------------------------------------------------------------
## SemiBin2 → contigs2bin
## ------------------------------------------------------------

: > "${SEMIBIN_C2B}"

while IFS= read -r -d '' fasta
do

    base=$(basename "${fasta}")

    base="${base%.gz}"
    base="${base%.fasta}"
    base="${base%.fa}"

    bin_id="semibin2.${base}"

    add_fasta_to_c2b \
        "${fasta}" \
        "${bin_id}" \
        "${SEMIBIN_C2B}"

done < <(

    find "${SEMIBIN_DIR}" \
        -maxdepth 1 \
        -type f \
        \( \
            -name '*.fa' \
            -o -name '*.fasta' \
            -o -name '*.fa.gz' \
        \) \
        -print0
)


## ------------------------------------------------------------
## MaxBin2 → contigs2bin
## ------------------------------------------------------------

: > "${MAXBIN_C2B}"

while IFS= read -r -d '' fasta
do

    base=$(basename "${fasta}" .fasta)

    bin_id="maxbin2.${base}"

    add_fasta_to_c2b \
        "${fasta}" \
        "${bin_id}" \
        "${MAXBIN_C2B}"

done < <(

    find "${MAXBIN_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "${SAMPLE}_flye_maxbin.[0-9][0-9][0-9].fasta" \
        -print0
)


## ------------------------------------------------------------
## Check that tables are non-empty
## ------------------------------------------------------------

for table in \
    "${METABAT_C2B}" \
    "${SEMIBIN_C2B}" \
    "${MAXBIN_C2B}"
do

    if [[ ! -s "${table}" ]]; then
        echo "ERROR: Empty contigs2bin table:"
        echo "${table}"
        exit 1
    fi

done


## ------------------------------------------------------------
## Check each binner assigns a contig at most once
## ------------------------------------------------------------

for table in \
    "${METABAT_C2B}" \
    "${SEMIBIN_C2B}" \
    "${MAXBIN_C2B}"
do

    DUP=$(cut -f1 "${table}" \
        | sort \
        | uniq -d \
        | wc -l)

    if [[ "${DUP}" -ne 0 ]]; then

        echo "ERROR: Duplicate contig assignments found in:"
        echo "${table}"
        echo "Duplicate contigs: ${DUP}"

        exit 1

    fi

done


## ------------------------------------------------------------
## Validate contig IDs against original Flye assembly
## ------------------------------------------------------------

ASSEMBLY_IDS="${OUTDIR}/assembly_contigs.txt"

awk '
    /^>/ {
        sub(/^>/, "", $0)
        split($0, a, /[[:space:]]+/)
        print a[1]
    }
' "${ASSEMBLY}" \
    | sort -u \
    > "${ASSEMBLY_IDS}"


for table in \
    "${METABAT_C2B}" \
    "${SEMIBIN_C2B}" \
    "${MAXBIN_C2B}"
do

    TMP="${table}.contigs.tmp"

    cut -f1 "${table}" \
        | sort -u \
        > "${TMP}"

    MISSING=$(comm -23 \
        "${TMP}" \
        "${ASSEMBLY_IDS}" \
        | wc -l)

    rm -f "${TMP}"

    if [[ "${MISSING}" -ne 0 ]]; then

        echo "ERROR: ${MISSING} contigs in:"
        echo "${table}"
        echo "were not found in the Flye assembly."

        exit 1

    fi

done


## ------------------------------------------------------------
## Input statistics
## ------------------------------------------------------------

N_METABAT=$(cut -f2 "${METABAT_C2B}" \
    | sort -u \
    | wc -l)

N_SEMIBIN=$(cut -f2 "${SEMIBIN_C2B}" \
    | sort -u \
    | wc -l)

N_MAXBIN=$(cut -f2 "${MAXBIN_C2B}" \
    | sort -u \
    | wc -l)


echo "============================================================"
echo "DAS Tool"
echo
echo "Sample:       ${SAMPLE}"
echo "Assembly:     ${ASSEMBLY}"
echo
echo "MetaBAT2 bins: ${N_METABAT}"
echo "SemiBin2 bins: ${N_SEMIBIN}"
echo "MaxBin2 bins:  ${N_MAXBIN}"
echo
echo "Output:       ${OUTDIR}"
echo "Threads:      ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


## ------------------------------------------------------------
## Run DAS Tool
##
## --write_bins:
##     export the final selected bins as FASTA
##
## --write_bin_evals:
##     retain DAS Tool's evaluation of all candidate bins
##
## score threshold 0.5:
##     explicit default for reproducibility
##
## diamond:
##     recommended search engine for larger datasets
## ------------------------------------------------------------

DAS_Tool \
    -i "${METABAT_C2B},${SEMIBIN_C2B},${MAXBIN_C2B}" \
    -l "MetaBAT2,SemiBin2,MaxBin2" \
    -c "${ASSEMBLY}" \
    -o "${PREFIX}" \
    --search_engine diamond \
    --score_threshold 0.5 \
    --write_bin_evals \
    --write_bins \
    --threads "${SLURM_CPUS_PER_TASK}"


## ------------------------------------------------------------
## Final output checks
## ------------------------------------------------------------

SUMMARY="${PREFIX}_DASTool_summary.tsv"

FINAL_BIN_DIR="${PREFIX}_DASTool_bins"


if [[ ! -s "${SUMMARY}" ]]; then

    echo "ERROR: DAS Tool summary not produced:"
    echo "${SUMMARY}"
    exit 1

fi


N_SELECTED=$(awk '
    NR > 1 {n++}
    END {print n+0}
' "${SUMMARY}")


N_FASTA=0

if [[ -d "${FINAL_BIN_DIR}" ]]; then

    N_FASTA=$(find "${FINAL_BIN_DIR}" \
        -maxdepth 1 \
        -type f \
        | wc -l)

fi


echo
echo "============================================================"
echo "DAS Tool completed"
echo
echo "Sample:          ${SAMPLE}"
echo "Selected bins:   ${N_SELECTED}"
echo "FASTA bins:      ${N_FASTA}"
echo
echo "Summary:"
echo "${SUMMARY}"
echo
echo "Final bins:"
echo "${FINAL_BIN_DIR}"
echo "============================================================"
