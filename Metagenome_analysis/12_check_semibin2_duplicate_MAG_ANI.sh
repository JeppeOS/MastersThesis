#!/bin/bash

#SBATCH --job-name=sb2_ani
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sb2_ani_%j.out
#SBATCH --error=logs/sb2_ani_%j.err
#SBATCH --account=methanotrophs


## ------------------------------------------------------------
## Environment
## ------------------------------------------------------------

source ~/.bashrc
conda activate gtdbtk_2.6.1

set -euo pipefail


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

MASTER="${PROJECT}/semibin2_analysis/MAG_master_table_SemiBin2.tsv"

OUTDIR="${PROJECT}/semibin2_analysis/duplicate_MAG_ANI"

rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}"


## ------------------------------------------------------------
## Locate FASTA path from Genome_ID
## ------------------------------------------------------------

get_fasta () {

    genome="$1"

    awk -F'\t' \
        -v genome="${genome}" \
        '
        NR == 1 {
            for (i=1; i<=NF; i++) {
                if ($i == "Genome_ID") genome_col=i
                if ($i == "Original_path") path_col=i
            }
            next
        }

        $genome_col == genome {
            print $path_col
            exit
        }
        ' "${MASTER}"
}


## ------------------------------------------------------------
## Copy/decompress FASTA
##
## SemiBin2 files may be gzip compressed even if their
## extension does not clearly indicate this.
## ------------------------------------------------------------

copy_fasta () {

    source="$1"
    dest="$2"

    magic=$(od \
        -An \
        -tx1 \
        -N2 \
        "${source}" \
        | tr -d ' \n')

    if [[ "${magic}" == "1f8b" ]]; then

        gzip -cd "${source}" > "${dest}"

    else

        cp "${source}" "${dest}"

    fi
}


## ------------------------------------------------------------
## Run one pairwise ANI comparison
##
## Arguments:
##   $1 = group name
##   $2 = genome 1
##   $3 = genome 2
## ------------------------------------------------------------

run_pair () {

    group="$1"
    genome1="$2"
    genome2="$3"

    echo
    echo "============================================================"
    echo "Group: ${group}"
    echo "============================================================"


    GROUPDIR="${OUTDIR}/${group}"

    mkdir -p "${GROUPDIR}"


    ## --------------------------------------------------------
    ## Find source FASTAs
    ## --------------------------------------------------------

    source1=$(get_fasta "${genome1}")
    source2=$(get_fasta "${genome2}")


    if [[ -z "${source1}" ]]; then

        echo "ERROR: Could not find FASTA for:"
        echo "${genome1}"

        exit 1

    fi


    if [[ -z "${source2}" ]]; then

        echo "ERROR: Could not find FASTA for:"
        echo "${genome2}"

        exit 1

    fi


    if [[ ! -s "${source1}" ]]; then

        echo "ERROR: FASTA does not exist:"
        echo "${source1}"

        exit 1

    fi


    if [[ ! -s "${source2}" ]]; then

        echo "ERROR: FASTA does not exist:"
        echo "${source2}"

        exit 1

    fi


    ## --------------------------------------------------------
    ## Make plain FASTA copies
    ## --------------------------------------------------------

    fasta1="${GROUPDIR}/${genome1}.fa"
    fasta2="${GROUPDIR}/${genome2}.fa"


    copy_fasta \
        "${source1}" \
        "${fasta1}"


    copy_fasta \
        "${source2}" \
        "${fasta2}"


    ## --------------------------------------------------------
    ## Build list for skani triangle
    ## --------------------------------------------------------

    LIST="${GROUPDIR}/genomes.txt"

    printf "%s\n%s\n" \
        "${fasta1}" \
        "${fasta2}" \
        > "${LIST}"


    OUTPUT="${GROUPDIR}/skani_triangle.tsv"


    ## --------------------------------------------------------
    ## Run skani
    ## --------------------------------------------------------

    skani triangle \
        -t "${SLURM_CPUS_PER_TASK}" \
        -E \
        -l "${LIST}" \
        > "${OUTPUT}"


    ## --------------------------------------------------------
    ## Print result
    ## --------------------------------------------------------

    echo
    echo "Genome 1:"
    echo "${genome1}"

    echo
    echo "Genome 2:"
    echo "${genome2}"

    echo
    echo "skani result:"
    cat "${OUTPUT}"
}


## ------------------------------------------------------------
## barcode10 - Methylumidiphilus
## ------------------------------------------------------------

run_pair \
    "barcode10_Methylumidiphilus" \
    "semibin2__flye__barcode10_seqs__SemiBin_1" \
    "semibin2__flye__barcode10_seqs__SemiBin_8"


## ------------------------------------------------------------
## barcode12 - Methylovulum
## ------------------------------------------------------------

run_pair \
    "barcode12_Methylovulum" \
    "semibin2__flye__barcode12_seqs__SemiBin_0" \
    "semibin2__flye__barcode12_seqs__SemiBin_6"


## ------------------------------------------------------------
## barcode14 - Methylobacter_A
## ------------------------------------------------------------

run_pair \
    "barcode14_Methylobacter_A" \
    "semibin2__flye__barcode14_seqs__SemiBin_0" \
    "semibin2__flye__barcode14_seqs__SemiBin_1"


## ------------------------------------------------------------
## barcode15 - Methylobacter_A
## ------------------------------------------------------------

run_pair \
    "barcode15_Methylobacter_A" \
    "semibin2__flye__barcode15_seqs__SemiBin_0" \
    "semibin2__flye__barcode15_seqs__SemiBin_11"


echo
echo "============================================================"
echo "All ANI comparisons completed"
echo
echo "Results:"
echo "${OUTDIR}"
echo "============================================================"
