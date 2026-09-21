#!/bin/bash

set -euo pipefail

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

OUTDIR="${PROJECT}/binner_comparison"
BINDIR="${OUTDIR}/input_bins"
MANIFEST="${OUTDIR}/bin_manifest.tsv"

rm -rf "${BINDIR}"

mkdir -p "${BINDIR}"

printf "Genome_ID\tBinner\tAssembly\tSample\tOriginal_path\n" \
    > "${MANIFEST}"


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
## SemiBin2
## ------------------------------------------------------------

for mode in flye strainberry
do

    for sample in "${SAMPLES[@]}"
    do

        dir="${PROJECT}/semibin2_results/${mode}/${sample}/output_bins"

        while IFS= read -r -d '' fasta
        do

            base=$(basename "${fasta}")

            base="${base%.gz}"
            base="${base%.fasta}"
            base="${base%.fa}"

            base=$(echo "${base}" \
                | sed 's/[^A-Za-z0-9._-]/_/g')

            genome_id="semibin2__${mode}__${sample}__${base}"

            dest="${BINDIR}/${genome_id}.fa"


            if [[ "${fasta}" == *.gz ]]; then

                gzip -cd "${fasta}" > "${dest}"

            else

                cp "${fasta}" "${dest}"

            fi


            printf "%s\t%s\t%s\t%s\t%s\n" \
                "${genome_id}" \
                "SemiBin2" \
                "${mode}" \
                "${sample}" \
                "${fasta}" \
                >> "${MANIFEST}"

        done < <(
            find "${dir}" \
                -maxdepth 1 \
                -type f \
                \( \
                    -name '*.fa' \
                    -o -name '*.fasta' \
                    -o -name '*.fa.gz' \
                \) \
                -print0
        )

    done

done


## ------------------------------------------------------------
## MaxBin2
## ------------------------------------------------------------

for mode in flye strainberry
do

    for sample in "${SAMPLES[@]}"
    do

        dir="${PROJECT}/maxbin2_results/${mode}/${sample}"

        while IFS= read -r -d '' fasta
        do

            base=$(basename "${fasta}" .fasta)

            base=$(echo "${base}" \
                | sed 's/[^A-Za-z0-9._-]/_/g')

            genome_id="maxbin2__${mode}__${sample}__${base}"

            dest="${BINDIR}/${genome_id}.fa"

            cp "${fasta}" "${dest}"


            printf "%s\t%s\t%s\t%s\t%s\n" \
                "${genome_id}" \
                "MaxBin2" \
                "${mode}" \
                "${sample}" \
                "${fasta}" \
                >> "${MANIFEST}"

        done < <(
            find "${dir}" \
                -maxdepth 1 \
                -type f \
                -name '*_maxbin.[0-9][0-9][0-9].fasta' \
                -print0
        )

    done

done


## ------------------------------------------------------------
## Validate
## ------------------------------------------------------------

N=$(find "${BINDIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)


echo "Prepared bins: ${N}"
echo "Manifest:      ${MANIFEST}"
echo "Input dir:     ${BINDIR}"


if [[ "${N}" -ne 612 ]]; then

    echo "ERROR: Expected 612 bins but found ${N}."
    exit 1

fi
