#!/bin/bash

#SBATCH --job-name=metabat2
#SBATCH --partition=normal
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/metabat_%A_%a.out
#SBATCH --error=logs/metabat_%A_%a.err
#SBATCH --account=methanotrophs


## Activate environment ##
source ~/.bashrc
conda activate ont_metagenome

set -euo pipefail

## Limit OpenMP programs to allocated CPUs ##
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT=~/methanotrophs/methanotroph_project/jeppe/metagenome

MANIFEST="${PROJECT}/trimmed_fastqs.txt"
FLYEDIR="${PROJECT}/flye_results"
MAPDIR="${PROJECT}/mapping"
OUTDIR="${PROJECT}/metabat2_results"

mkdir -p "${OUTDIR}"


## ------------------------------------------------------------
## Get sample for this array task
## ------------------------------------------------------------

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")

if [[ -z "${FASTQ}" ]]; then
    echo "ERROR: No FASTQ found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


## barcode15_seqs_trimmed.fastq.gz -> barcode15_seqs ##
BASENAME=$(basename "${FASTQ}")
SAMPLE=${BASENAME%_trimmed.fastq.gz}


## ------------------------------------------------------------
## Input files
## ------------------------------------------------------------

ASSEMBLY="${FLYEDIR}/${SAMPLE}/assembly.fasta"
BAM="${MAPDIR}/${SAMPLE}/${SAMPLE}.sorted.bam"


## Check assembly ##
if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "ERROR: Assembly not found:"
    echo "${ASSEMBLY}"
    exit 1
fi


## Check BAM ##
if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found:"
    echo "${BAM}"
    exit 1
fi


## ------------------------------------------------------------
## Output directories/files
## ------------------------------------------------------------

SAMPLEDIR="${OUTDIR}/${SAMPLE}"
BINDIR="${SAMPLEDIR}/bins"

DEPTH="${SAMPLEDIR}/${SAMPLE}_depth.txt"

mkdir -p "${BINDIR}"


echo "============================================================"
echo "SLURM job:    ${SLURM_JOB_ID}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Sample:       ${SAMPLE}"
echo "Assembly:     ${ASSEMBLY}"
echo "BAM:          ${BAM}"
echo "Depth file:   ${DEPTH}"
echo "Bin dir:      ${BINDIR}"
echo "Threads:      ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


## ------------------------------------------------------------
## Calculate contig depth
##
## Use 97% minimum alignment identity --> Based on 
## This filters lower-identity read mappings before calculating
## coverage used by MetaBAT2.
## ------------------------------------------------------------

jgi_summarize_bam_contig_depths \
    --outputDepth "${DEPTH}" \
    --percentIdentity 97 \
    --referenceFasta "${ASSEMBLY}" \
    "${BAM}"


## Check depth file ##
if [[ ! -s "${DEPTH}" ]]; then
    echo "ERROR: Depth file was not created or is empty."
    exit 1
fi


## ------------------------------------------------------------
## Run MetaBAT2
##
## -m 2500:
## Only contigs >= 2500 bp are considered for binning.
##
## --seed 42:
## Reproducible binning.
##
## --saveCls:
## Save contig-to-bin classification information.
##
## --unbinned:
## Also save contigs that were not assigned to a bin.
## ------------------------------------------------------------

metabat2 \
    -i "${ASSEMBLY}" \
    -a "${DEPTH}" \
    -o "${BINDIR}/${SAMPLE}_bin" \
    -m 2500 \
    -t "${SLURM_CPUS_PER_TASK}" \
    --seed 42 \
    --saveCls \
    --unbinned \
    --verbose \
    > "${SAMPLEDIR}/${SAMPLE}_metabat2.log" 2>&1


## ------------------------------------------------------------
## Summary
## ------------------------------------------------------------

NBINS=$(find "${BINDIR}" \
    -maxdepth 1 \
    -type f \
    -name "${SAMPLE}_bin.*.fa" \
    ! -name "*unbinned*" \
    | wc -l)


echo "============================================================"
echo "MetaBAT2 completed successfully"
echo "Sample:         ${SAMPLE}"
echo "Number of bins: ${NBINS}"
echo "============================================================"