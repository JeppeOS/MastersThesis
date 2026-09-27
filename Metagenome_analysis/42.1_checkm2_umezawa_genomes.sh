#!/bin/bash

#SBATCH --job-name=checkm2_umezawa
#SBATCH --account=methanotrophs
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=comparative_analysis/final_species_phylogeny/logs/89_checkm2_%j.out
#SBATCH --error=comparative_analysis/final_species_phylogeny/logs/89_checkm2_%j.err


set -euo pipefail


## ================================================================== ##
## STAGE 89
##
## CheckM2 quality for the two Umezawa Allocrenothrix genomes.
##
## Uses the SAME CheckM2 v1.1.0 database as the previous 631-genome
## quality analysis.
## ================================================================== ##


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

BASE="$PROJECT/comparative_analysis/final_species_phylogeny"

SOURCE="$BASE/00_umezawa_genomes"

CHECKM_BASE="$BASE/06_checkm2_umezawa"

INPUT="$CHECKM_BASE/00_input"

RESULTS="$CHECKM_BASE/01_results"


DB="$HOME/methanotrophs/methanotroph_project/jeppe/reference_data/checkm2_v1.1.0_db/CheckM2_database/uniref100.KO.1.dmnd"


mkdir -p "$INPUT"

rm -rf "$RESULTS"


## ================================================================== ##
## Copy only the two genome FASTAs
## ================================================================== ##

cp \
    "$SOURCE/UME_AF98_Allocrenothrix_methanica.fna" \
    "$INPUT/"


cp \
    "$SOURCE/UME_MI19235_Allocrenothrix_methanica.fna" \
    "$INPUT/"


## ================================================================== ##
## Environment
## ================================================================== ##

source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda activate checkm2_cpu


echo
echo "CheckM2 version:"
checkm2 --version

echo
echo "Database:"
echo "$DB"

echo
echo "Input genomes:"
ls -lh "$INPUT"
echo


## ================================================================== ##
## Run CheckM2
## ================================================================== ##

checkm2 predict \
    --input "$INPUT" \
    --output-directory "$RESULTS" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --database_path "$DB" \
    --extension fna


## ================================================================== ##
## QC
## ================================================================== ##

REPORT="$RESULTS/quality_report.tsv"


if [[ ! -s "$REPORT" ]]
then
    echo "ERROR: CheckM2 quality_report.tsv was not produced."
    exit 1
fi


N=$(($(wc -l < "$REPORT") - 1))


if [[ "$N" -ne 2 ]]
then
    echo "ERROR: expected 2 CheckM2 genome rows; found $N."
    exit 1
fi


echo
echo "============================================================"
echo "STAGE 89 - CHECKM2 RESULTS"
echo "============================================================"
echo

column -t -s $'\t' "$REPORT"

echo
echo "Quality report:"
echo "  $REPORT"
echo
echo "Stage 89 complete."
