#!/bin/bash

set -euo pipefail

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome
PROFILE_DIR="${PROJECT}/sylph_results/tax_profiles"
OUTDIR="${PROJECT}/abundance_comparison"
OUTFILE="${OUTDIR}/sylph_sequence_abundance_merged.tsv"

mkdir -p "${OUTDIR}"

shopt -s nullglob

FILES=(
    "${PROFILE_DIR}"/*.sylphmpa
)

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .sylphmpa files found in:"
    echo "${PROFILE_DIR}"
    exit 1
fi

echo "Found ${#FILES[@]} Sylph profiles"

sylph-tax merge \
    "${FILES[@]}" \
    --column sequence_abundance \
    -o "${OUTFILE}"

echo "Written:"
echo "${OUTFILE}"
