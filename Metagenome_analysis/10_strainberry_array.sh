#!/bin/bash

#SBATCH --job-name=strainberry
#SBATCH --partition=normal
#SBATCH --array=1-7%3
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/strainberry_%A_%a.out
#SBATCH --error=logs/strainberry_%A_%a.err
#SBATCH --account=methanotrophs

source ~/.bashrc
conda activate sberry

set -euo pipefail

## Strainberry executable ##
export PATH="$HOME/software/strainberry:$PATH"


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

FLYEDIR="${PROJECT}/flye_results"
MAPDIR="${PROJECT}/mapping"
OUTBASE="${PROJECT}/strainberry_results"

mkdir -p "${OUTBASE}"


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

SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"


## ------------------------------------------------------------
## Inputs
## ------------------------------------------------------------

ASSEMBLY="${FLYEDIR}/${SAMPLE}/assembly.fasta"

BAM="${MAPDIR}/${SAMPLE}/${SAMPLE}.sorted.bam"

OUTDIR="${OUTBASE}/${SAMPLE}"


## ------------------------------------------------------------
## Validate inputs
## ------------------------------------------------------------

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found:"
    echo "${ASSEMBLY}"
    exit 1
fi

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found:"
    echo "${BAM}"
    exit 1
fi


## ------------------------------------------------------------
## Check BAM
## ------------------------------------------------------------

samtools quickcheck -v "${BAM}"


## ------------------------------------------------------------
## Index assembly if necessary
## ------------------------------------------------------------

if [[ ! -f "${ASSEMBLY}.fai" ]]; then
    echo "Creating FASTA index..."
    samtools faidx "${ASSEMBLY}"
fi


## ------------------------------------------------------------
## Index BAM if necessary
## ------------------------------------------------------------

if [[ ! -f "${BAM}.bai" ]]; then
    echo "Creating BAM index..."
    samtools index \
        -@ "${SLURM_CPUS_PER_TASK}" \
        "${BAM}"
fi


## ------------------------------------------------------------
## Prevent overwrite
## ------------------------------------------------------------

if [[ -e "${OUTDIR}" ]]; then
    echo "ERROR: Output directory already exists:"
    echo "${OUTDIR}"
    echo "Remove or rename it before rerunning."
    exit 1
fi


## ------------------------------------------------------------
## Summary
## ------------------------------------------------------------

echo "============================================================"
echo "Job:        ${SLURM_JOB_ID}"
echo "Task:       ${SLURM_ARRAY_TASK_ID}"
echo "Sample:     ${SAMPLE}"
echo "Assembly:   ${ASSEMBLY}"
echo "BAM:        ${BAM}"
echo "Output:     ${OUTDIR}"
echo "CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


## ------------------------------------------------------------
## Run Strainberry
## ------------------------------------------------------------

strainberry \
    --nanopore \
    -r "${ASSEMBLY}" \
    -b "${BAM}" \
    -o "${OUTDIR}" \
    -c "${SLURM_CPUS_PER_TASK}" \
    -v


## ------------------------------------------------------------
## Check output
## ------------------------------------------------------------

FINAL="${OUTDIR}/assembly.scaffolds.fa"

if [[ ! -s "${FINAL}" ]]; then
    echo "ERROR: Final Strainberry assembly not found:"
    echo "${FINAL}"
    exit 1
fi


N_SCAFFOLDS=$(grep -c '^>' "${FINAL}")

TOTAL_BP=$(awk '
    /^>/ {next}
    {
        gsub(/[[:space:]]/, "")
        n += length($0)
    }
    END {
        print n
    }
' "${FINAL}")


echo "============================================================"
echo "Strainberry completed successfully"
echo "Sample:       ${SAMPLE}"
echo "Scaffolds:    ${N_SCAFFOLDS}"
echo "Total bp:     ${TOTAL_BP}"
echo "Final FASTA:  ${FINAL}"
echo "============================================================"
