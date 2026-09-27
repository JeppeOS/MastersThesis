#!/bin/bash

#SBATCH --job-name=MCA0421_umezawa
#SBATCH --account=methanotrophs
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=comparative_analysis/final_species_phylogeny/logs/90_MCA0421_%j.out
#SBATCH --error=comparative_analysis/final_species_phylogeny/logs/90_MCA0421_%j.err


set -euo pipefail


## ================================================================== ##
## STAGE 90
##
## TARGETED MCA0421 HOMOLOG SEARCH
##
## Query genomes:
##
##   AF98
##   MI19235
##
## Reference:
##
##   Methylococcus capsulatus Bath MCA0421
##   WP_041360666.1 / MCA_RS02090
##
## Previous workflow established that MCA0421 maps to:
##
##   Cluster_00063
##   Module_06
##
## Search threshold:
##
##   >= 40% sequence identity
##   >= 80% coverage
##
## No broad annotation is performed.
## ================================================================== ##


ROOT="$HOME/methanotrophs/methanotroph_project/jeppe"

PROJECT="$ROOT/metagenome"

OLDWF="$ROOT/genome_analysis_workflow"

BASE="$PROJECT/comparative_analysis/final_species_phylogeny"

GENOMES="$BASE/00_umezawa_genomes"

OUT="$BASE/07_MCA0421_umezawa"

PROTEINS="$OUT/00_proteins"

SEARCH="$OUT/01_search"

TMP="$OUT/tmp"


REFERENCE="$OLDWF/19_Methylococcus_OMC/reference/MCA0421.faa"


mkdir -p \
    "$PROTEINS" \
    "$SEARCH" \
    "$TMP"


## ================================================================== ##
## QC reference
## ================================================================== ##

if [[ ! -s "$REFERENCE" ]]
then
    echo "ERROR: MCA0421 reference FASTA not found:"
    echo "$REFERENCE"
    exit 1
fi


echo
echo "MCA0421 reference:"
echo "$REFERENCE"
grep '^>' "$REFERENCE"
echo


## ================================================================== ##
## Predict proteins with Prodigal
##
## These are complete isolate genomes rather than MAGs, so -p single
## is appropriate.
## ================================================================== ##

source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda activate checkm2_cpu


for SAMPLE in \
    UME_AF98_Allocrenothrix_methanica \
    UME_MI19235_Allocrenothrix_methanica

do

    GENOME="$GENOMES/${SAMPLE}.fna"

    FAA="$PROTEINS/${SAMPLE}.faa"


    if [[ ! -s "$GENOME" ]]
    then
        echo "ERROR: missing genome:"
        echo "$GENOME"
        exit 1
    fi


    echo "Predicting proteins: $SAMPLE"


    prodigal \
        -i "$GENOME" \
        -a "$FAA" \
        -p single \
        -q


    N=$(grep -c '^>' "$FAA")


    echo "  proteins: $N"

done


## ================================================================== ##
## Combine proteins
## ================================================================== ##

ALLFAA="$PROTEINS/UMEZAWA_all_proteins.faa"

cat \
    "$PROTEINS/UME_AF98_Allocrenothrix_methanica.faa" \
    "$PROTEINS/UME_MI19235_Allocrenothrix_methanica.faa" \
    > "$ALLFAA"


## ================================================================== ##
## MMseqs2 search
## ================================================================== ##

conda activate MMseqs2


RESULTS="$SEARCH/MCA0421_hits.tsv"


mmseqs easy-search \
    "$ALLFAA" \
    "$REFERENCE" \
    "$RESULTS" \
    "$TMP" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --min-seq-id 0.40 \
    -c 0.80 \
    --cov-mode 0 \
    --format-output \
    "query,target,fident,alnlen,qcov,tcov,evalue,bits"


## ================================================================== ##
## Genome-level presence / absence summary
## ================================================================== ##

SUMMARY="$SEARCH/MCA0421_presence_absence.tsv"


python - \
    "$RESULTS" \
    "$SUMMARY" \
<<'PY'
import csv
import sys


hits_file = sys.argv[1]
out_file = sys.argv[2]


genomes = {
    "UME_AF98_Allocrenothrix_methanica": [],
    "UME_MI19235_Allocrenothrix_methanica": [],
}


with open(hits_file) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t",
        fieldnames=[
            "query",
            "target",
            "fident",
            "alnlen",
            "qcov",
            "tcov",
            "evalue",
            "bits",
        ],
    )


    for row in reader:

        query = row["query"]


        if query.startswith("AF98_"):

            genome = "UME_AF98_Allocrenothrix_methanica"

        elif query.startswith("MI19235_"):

            genome = "UME_MI19235_Allocrenothrix_methanica"

        else:

            raise RuntimeError(
                f"Cannot assign protein to genome: {query}"
            )


        genomes[genome].append(row)


with open(out_file, "w", newline="") as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
    )


    writer.writerow([
        "genome",
        "MCA0421_present",
        "n_passing_hits",
        "best_protein",
        "best_fident",
        "best_qcov",
        "best_tcov",
        "best_evalue",
        "best_bits",
    ])


    for genome, hits in genomes.items():

        if hits:

            best = max(
                hits,
                key=lambda x: float(x["bits"])
            )


            writer.writerow([
                genome,
                1,
                len(hits),
                best["query"],
                best["fident"],
                best["qcov"],
                best["tcov"],
                best["evalue"],
                best["bits"],
            ])

        else:

            writer.writerow([
                genome,
                0,
                0,
                "",
                "",
                "",
                "",
                "",
                "",
            ])
PY


## ================================================================== ##
## Report
## ================================================================== ##

echo
echo "============================================================"
echo "STAGE 90 - MCA0421 SEARCH"
echo "============================================================"
echo

echo "Passing protein hits:"
echo

if [[ -s "$RESULTS" ]]
then
    column -t -s $'\t' "$RESULTS"
else
    echo "None."
fi


echo
echo "Presence / absence:"
echo

column -t -s $'\t' "$SUMMARY"

echo
echo "Outputs:"
echo "  $RESULTS"
echo "  $SUMMARY"
echo
echo "Stage 90 complete."
