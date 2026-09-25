#!/bin/bash

#SBATCH --job-name=final_tree_identify
#SBATCH --account=methanotrophs
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=comparative_analysis/final_species_phylogeny/logs/85_identify_%j.out
#SBATCH --error=comparative_analysis/final_species_phylogeny/logs/85_identify_%j.err


set -euo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

BASE="$PROJECT/comparative_analysis/final_species_phylogeny"

MANIFEST="$BASE/01_inputs/genomes_650.tsv"

OUTDIR="$BASE/02_gtdbtk_identify"


source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda activate gtdbtk_2.6.1


export GTDBTK_DATA_PATH="$HOME/methanotrophs/methanotroph_project/jeppe/reference_data/gtdbtk_r226"


echo "GTDB-Tk version:"
gtdbtk --version

echo
echo "GTDB data:"
echo "$GTDBTK_DATA_PATH"

echo
echo "Input genomes:"
wc -l "$MANIFEST"


rm -rf "$OUTDIR"


gtdbtk identify \
    --batchfile "$MANIFEST" \
    --out_dir "$OUTDIR" \
    --cpus "$SLURM_CPUS_PER_TASK"


echo
echo "Stage 85 complete."
