#!/bin/bash
set -euo pipefail


## ================================================================== ##
## STAGE 84
##
## PREPARE FINAL bac120 SPECIES-TREE INPUT
##
##   631 GlobDB Methylococcales
##    16 SemiBin2 Methylococcales MAGs
##     2 Umezawa Allocrenothrix genomes
##     1 Methylophaga outgroup
##
## Expected total = 650 genomes
## ================================================================== ##


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

OLD_WORKFLOW="$HOME/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"

OLD_PHYLO="$OLD_WORKFLOW/17_phylogeny/03_bac120_fasttree"

BASE="$PROJECT/comparative_analysis/final_species_phylogeny"

INPUT="$BASE/01_inputs"

mkdir -p "$INPUT"


OLD_MANIFEST="$OLD_PHYLO/00_inputs/genomes_632.tsv"

OUTGROUP_FASTA="$OLD_PHYLO/00_inputs/outgroup/OUT_Methylophaga_nitratireducenticrescens.fna"

UMEZAWA="$BASE/00_umezawa_genomes"

MANIFEST="$INPUT/genomes_650.tsv"


## ================================================================== ##
## Input checks
## ================================================================== ##

[[ -s "$OLD_MANIFEST" ]] || {
    echo "ERROR: old 632-genome manifest not found:"
    echo "$OLD_MANIFEST"
    exit 1
}


[[ -s "$OUTGROUP_FASTA" ]] || {
    echo "ERROR: old Methylophaga outgroup FASTA not found:"
    echo "$OUTGROUP_FASTA"
    exit 1
}


[[ -s "$UMEZAWA/UME_AF98_Allocrenothrix_methanica.fna" ]] || {
    echo "ERROR: AF98 genome missing. Run Stage 83 first."
    exit 1
}


[[ -s "$UMEZAWA/UME_MI19235_Allocrenothrix_methanica.fna" ]] || {
    echo "ERROR: MI19235 genome missing. Run Stage 83 first."
    exit 1
}


## ================================================================== ##
## Recover original 631 GlobDB genomes
##
## Preserve exactly the same FASTA paths and IDs as the previous tree.
## ================================================================== ##

awk -F '\t' '
    BEGIN {
        OFS="\t"
    }

    NF >= 2 &&
    $2 != "OUT_Methylophaga_nitratireducenticrescens" {
        print $1, $2
    }
' "$OLD_MANIFEST" \
    > "$INPUT/globdb_631.tsv"


N_GLOBDB=$(wc -l < "$INPUT/globdb_631.tsv")


if [[ "$N_GLOBDB" -ne 631 ]]
then
    echo "ERROR: expected 631 GlobDB genomes; recovered $N_GLOBDB."
    exit 1
fi


cp "$INPUT/globdb_631.tsv" "$MANIFEST"


## ================================================================== ##
## Add the 16 SemiBin2 Methylococcales MAGs
## ================================================================== ##

add_mag () {

    SAMPLE="$1"
    BIN="$2"
    ID="$3"

    FASTA="$PROJECT/semibin2_results/flye/${SAMPLE}_seqs/output_bins/SemiBin_${BIN}.fa.gz"


    if [[ ! -s "$FASTA" ]]
    then
        echo "ERROR: MAG FASTA missing:"
        echo "$FASTA"
        exit 1
    fi


    printf "%s\t%s\n" \
        "$FASTA" \
        "$ID" \
        >> "$MANIFEST"
}


add_mag barcode10 14 MAG_b10_SB14
add_mag barcode10 1  MAG_b10_SB1
add_mag barcode10 8  MAG_b10_SB8

add_mag barcode12 5  MAG_b12_SB5
add_mag barcode12 0  MAG_b12_SB0
add_mag barcode12 6  MAG_b12_SB6

add_mag barcode13 1  MAG_b13_SB1

add_mag barcode14 1  MAG_b14_SB1
add_mag barcode14 0  MAG_b14_SB0

add_mag barcode15 0  MAG_b15_SB0
add_mag barcode15 11 MAG_b15_SB11
add_mag barcode15 3  MAG_b15_SB3
add_mag barcode15 10 MAG_b15_SB10

add_mag barcode16 1  MAG_b16_SB1
add_mag barcode16 12 MAG_b16_SB12
add_mag barcode16 11 MAG_b16_SB11


## ================================================================== ##
## Add Umezawa et al. isolates
## ================================================================== ##

printf "%s\t%s\n" \
    "$UMEZAWA/UME_AF98_Allocrenothrix_methanica.fna" \
    "UME_AF98_Allocrenothrix_methanica" \
    >> "$MANIFEST"


printf "%s\t%s\n" \
    "$UMEZAWA/UME_MI19235_Allocrenothrix_methanica.fna" \
    "UME_MI19235_Allocrenothrix_methanica" \
    >> "$MANIFEST"


## ================================================================== ##
## Add the SAME outgroup as the previous tree
## ================================================================== ##

printf "%s\t%s\n" \
    "$OUTGROUP_FASTA" \
    "OUT_Methylophaga_nitratireducenticrescens" \
    >> "$MANIFEST"


## ================================================================== ##
## QC
## ================================================================== ##

TOTAL=$(wc -l < "$MANIFEST")

UNIQUE=$(cut -f2 "$MANIFEST" | sort -u | wc -l)


echo
echo "============================================================"
echo "FINAL PHYLOGENY INPUT"
echo "============================================================"
echo
echo "GlobDB genomes:       631"
echo "SemiBin2 MAGs:         16"
echo "Umezawa genomes:        2"
echo "Outgroup:               1"
echo "--------------------------"
echo "Expected total:       650"
echo
echo "Observed total:       $TOTAL"
echo "Unique IDs:           $UNIQUE"
echo


if [[ "$TOTAL" -ne 650 ]]
then
    echo "ERROR: expected 650 manifest entries."
    exit 1
fi


if [[ "$UNIQUE" -ne 650 ]]
then
    echo "ERROR: duplicate genome IDs detected."
    exit 1
fi


echo "Checking FASTA paths..."

while IFS=$'\t' read -r FASTA ID
do
    if [[ ! -s "$FASTA" ]]
    then
        echo "ERROR: missing FASTA for $ID"
        echo "$FASTA"
        exit 1
    fi
done < "$MANIFEST"


echo
echo "All 650 FASTA files found."
echo
echo "Manifest:"
echo "  $MANIFEST"
