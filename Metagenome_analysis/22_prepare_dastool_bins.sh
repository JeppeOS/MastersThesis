#!/bin/bash

set -euo pipefail

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

DASROOT="${PROJECT}/dastool_results/flye"

OUTDIR="${PROJECT}/dastool_evaluation"
BINDIR="${OUTDIR}/input_bins"
MANIFEST="${OUTDIR}/bin_manifest.tsv"

rm -rf "${OUTDIR}"

mkdir -p "${BINDIR}"

printf "Genome_ID\tAssembly\tBinner\tSample\tOriginal_path\n" \
    > "${MANIFEST}"


SAMPLES=(
    barcode10_seqs
    barcode11_seqs
    barcode12_seqs
    barcode13_seqs
    barcode14_seqs
    barcode15_seqs
    barcode16_seqs
)


for sample in "${SAMPLES[@]}"
do

    bindir="${DASROOT}/${sample}/${sample}_DASTool_bins"

    if [[ ! -d "${bindir}" ]]; then
        echo "ERROR: Missing DAS Tool bin directory:"
        echo "${bindir}"
        exit 1
    fi


    while IFS= read -r -d '' fasta
    do

        base=$(basename "${fasta}")

        base="${base%.fasta}"
        base="${base%.fna}"
        base="${base%.fa}"

        base=$(echo "${base}" \
            | sed 's/[^A-Za-z0-9._-]/_/g')

        genome_id="dastool__flye__${sample}__${base}"

        dest="${BINDIR}/${genome_id}.fa"

        cp "${fasta}" "${dest}"


        printf "%s\t%s\t%s\t%s\t%s\n" \
            "${genome_id}" \
            "Flye" \
            "DAS_Tool" \
            "${sample}" \
            "${fasta}" \
            >> "${MANIFEST}"

    done < <(
        find "${bindir}" \
            -maxdepth 1 \
            -type f \
            \( \
                -name '*.fa' \
                -o -name '*.fasta' \
                -o -name '*.fna' \
            \) \
            -print0
    )

done


N=$(find "${BINDIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)


echo "============================================================"
echo "Prepared DAS Tool bins: ${N}"
echo "Input directory:        ${BINDIR}"
echo "Manifest:               ${MANIFEST}"
echo "============================================================"
