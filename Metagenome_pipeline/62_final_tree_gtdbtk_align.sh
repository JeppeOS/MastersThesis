#!/bin/bash

#SBATCH --job-name=final_tree_align
#SBATCH --account=methanotrophs
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=comparative_analysis/final_species_phylogeny/logs/86_align_%j.out
#SBATCH --error=comparative_analysis/final_species_phylogeny/logs/86_align_%j.err


set -euo pipefail


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

BASE="$PROJECT/comparative_analysis/final_species_phylogeny"

IDENTIFY="$BASE/02_gtdbtk_identify"

OUTDIR="$BASE/03_gtdbtk_align"


source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda activate gtdbtk_2.6.1


export GTDBTK_DATA_PATH="$HOME/methanotrophs/methanotroph_project/jeppe/reference_data/gtdbtk_r226"


rm -rf "$OUTDIR"


gtdbtk align \
    --identify_dir "$IDENTIFY" \
    --out_dir "$OUTDIR" \
    --skip_gtdb_refs \
    --cpus "$SLURM_CPUS_PER_TASK"


MSA="$OUTDIR/align/gtdbtk.bac120.user_msa.fasta.gz"


echo
echo "============================================================"
echo "bac120 ALIGNMENT QC"
echo "============================================================"
echo


python - "$MSA" <<'PY'
import gzip
import sys

path = sys.argv[1]

seqs = {}
name = None
parts = []

with gzip.open(path, "rt") as fh:

    for line in fh:

        line = line.strip()

        if line.startswith(">"):

            if name is not None:
                seqs[name] = "".join(parts)

            name = line[1:].split()[0]
            parts = []

        else:
            parts.append(line)

    if name is not None:
        seqs[name] = "".join(parts)


lengths = sorted(
    set(
        map(
            len,
            seqs.values()
        )
    )
)


print("Taxa:", len(seqs))
print("Alignment lengths:", lengths)

if len(lengths) != 1:
    raise SystemExit(
        "ERROR: sequences do not have a common alignment length."
    )


wanted = {
    "UME_AF98_Allocrenothrix_methanica",
    "UME_MI19235_Allocrenothrix_methanica",
    "MAG_b16_SB11",
    "OUT_Methylophaga_nitratireducenticrescens",
}

print()
print("Key taxa:")

for taxon in sorted(wanted):
    print(
        taxon,
        "PRESENT" if taxon in seqs else "MISSING"
    )
PY


echo
echo "MSA:"
echo "$MSA"
echo
echo "Stage 86 complete."
