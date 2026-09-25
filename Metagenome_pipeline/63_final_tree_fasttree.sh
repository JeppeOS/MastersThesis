#!/bin/bash

#SBATCH --job-name=final_bac120_fasttree
#SBATCH --account=methanotrophs
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=comparative_analysis/final_species_phylogeny/logs/87_fasttree_%j.out
#SBATCH --error=comparative_analysis/final_species_phylogeny/logs/87_fasttree_%j.err


set -euo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

BASE="$PROJECT/comparative_analysis/final_species_phylogeny"


MSA="$BASE/03_gtdbtk_align/align/gtdbtk.bac120.user_msa.fasta.gz"

OUTDIR="$BASE/04_fasttree"


mkdir -p "$OUTDIR"


TREE="$OUTDIR/final_methylococcales_GlobDB_MAGs_Umezawa_bac120.tree"

FTLOG="$OUTDIR/final_methylococcales_GlobDB_MAGs_Umezawa_bac120.fasttree.log"


source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda activate gtdbtk_2.6.1


export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"


echo "FastTree:"
which FastTreeMP

echo
echo "Threads:"
echo "$OMP_NUM_THREADS"

echo
echo "Input MSA:"
echo "$MSA"

echo


zcat "$MSA" \
    | FastTreeMP \
        -wag \
        -gamma \
        -boot 1000 \
        -log "$FTLOG" \
    > "$TREE"


echo
echo "============================================================"
echo "FASTTREE COMPLETE"
echo "============================================================"
echo
echo "Tree:"
echo "$TREE"
echo
echo "Log:"
echo "$FTLOG"
