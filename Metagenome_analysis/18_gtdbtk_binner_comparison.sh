#!/bin/bash

#SBATCH --job-name=bin_gtdbtk
#SBATCH --partition=normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bin_gtdbtk_%j.out
#SBATCH --error=logs/bin_gtdbtk_%j.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate gtdbtk_2.6.1

set -euo pipefail


PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

INPUT="${PROJECT}/binner_comparison/input_bins"

OUTDIR="${PROJECT}/binner_comparison/gtdbtk"

BATCH="${PROJECT}/binner_comparison/gtdbtk_batch.tsv"


export GTDBTK_DATA_PATH=\
/home/jeppeos/methanotrophs/methanotroph_project/jeppe/reference_data/gtdbtk_r226


## ------------------------------------------------------------
## Build batch file
## ------------------------------------------------------------

: > "${BATCH}"


find "${INPUT}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | sort \
    | while read -r fasta
do

    genome_id=$(basename "${fasta}" .fa)

    printf "%s\t%s\n" \
        "${fasta}" \
        "${genome_id}" \
        >> "${BATCH}"

done


N=$(wc -l < "${BATCH}")

echo "GTDB-Tk input bins: ${N}"

if [[ "${N}" -ne 612 ]]; then
    echo "ERROR: Expected 612 bins."
    exit 1
fi


rm -rf "${OUTDIR}"


gtdbtk classify_wf \
    --batchfile "${BATCH}" \
    --out_dir "${OUTDIR}" \
    --cpus "${SLURM_CPUS_PER_TASK}"


echo
echo "GTDB-Tk complete."


## ------------------------------------------------------------
## Pull out Methylococcales immediately
## ------------------------------------------------------------

BAC="${OUTDIR}/gtdbtk.bac120.summary.tsv"

METH="${PROJECT}/binner_comparison/Methylococcales_GTDB.tsv"


if [[ -s "${BAC}" ]]; then

    awk -F'\t' '
        NR == 1 ||
        $0 ~ /o__Methylococcales/
    ' "${BAC}" > "${METH}"

    echo
    echo "Methylococcales:"
    awk 'NR > 1 {n++} END {print n+0}' "${METH}"

fi
