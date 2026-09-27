#!/usr/bin/env bash

#SBATCH --job-name=fegenie_semibin2
#SBATCH --partition=normal
#SBATCH --account=methanotrophs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fegenie_semibin2_%j.out
#SBATCH --error=logs/fegenie_semibin2_%j.err


set -Eeo pipefail
shopt -s nullglob
umask 0027


###############################################################
## User settings ##
###############################################################

BASE_DIR="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

INPUT_ROOT="$BASE_DIR/functional_analysis/input_bins"


## ------------------------------------------------------------
## This is a completely independent FeGenie run.
##
## Change the name if we intentionally rerun it later.
## ------------------------------------------------------------

RUN_ROOT="$BASE_DIR/functional_analysis/FeGenie_SemiBin2_full_run_20260918"

FEGENIE_OUT="$RUN_ROOT/fegenie_results"

FAA_OUT="$RUN_ROOT/predicted_proteins_faa"

METADATA_OUT="$RUN_ROOT/run_metadata"


CONDA_ENV="fegenie_conda"

EXPECTED_INPUT_COUNT=187


###############################################################
## Small helper functions ##
###############################################################

log() {

    printf \
        '[%s] %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$*"
}


fail() {

    printf \
        '[%s] ERROR: %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$*" \
        >&2

    exit 1
}


###############################################################
## Activate Conda ##
###############################################################

CONDA_SH=""


for candidate in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/mambaforge/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh"
do

    if [[ -f "$candidate" ]]; then

        CONDA_SH="$candidate"

        break
    fi

done


if [[ -n "$CONDA_SH" ]]; then

    # shellcheck source=/dev/null
    source "$CONDA_SH"

elif command -v conda >/dev/null 2>&1; then

    # shellcheck source=/dev/null
    source "$(conda info --base)/etc/profile.d/conda.sh"

else

    fail \
        "Conda was not found. Edit the Conda activation section."

fi


conda activate "$CONDA_ENV"


###############################################################
## Preflight checks ##
###############################################################

for program in \
    FeGenie.py \
    prodigal \
    hmmsearch \
    diamond \
    blastp
do

    command -v "$program" >/dev/null 2>&1 \
        || fail \
        "Required program is not available in the active environment: $program"

done


[[ -d "$INPUT_ROOT" ]] \
    || fail \
    "Input directory does not exist: $INPUT_ROOT"


## ------------------------------------------------------------
## SemiBin2 inputs are stored as:
##
## input_bins/
## ├── barcode10_seqs/*.fa
## ├── barcode11_seqs/*.fa
## ...
## └── barcode16_seqs/*.fa
## ------------------------------------------------------------

inputs=(
    "$INPUT_ROOT"/*/*.fa
)


input_count=${#inputs[@]}


(( input_count > 0 )) \
    || fail \
    "No input FASTAs matching $INPUT_ROOT/*/*.fa were found."


if (( input_count != EXPECTED_INPUT_COUNT )); then

    log \
        "WARNING: found $input_count input files; expected $EXPECTED_INPUT_COUNT."

    log \
        "The run will continue, but check for missing or stray MAGs."

fi


## ------------------------------------------------------------
## FeGenie HMM environment
## ------------------------------------------------------------

[[ -n "${iron_hmms:-}" ]] \
    || fail \
    "The FeGenie environment variable iron_hmms is not set after Conda activation."


[[ -d "$iron_hmms" ]] \
    || fail \
    "The FeGenie HMM directory does not exist: $iron_hmms"


[[ -f "$iron_hmms/HMM-bitcutoffs.txt" ]] \
    || fail \
    "HMM-bitcutoffs.txt is missing from: $iron_hmms"


## ------------------------------------------------------------
## Refuse to overwrite an earlier run
## ------------------------------------------------------------

[[ ! -e "$RUN_ROOT" ]] \
    || fail \
    "Run directory already exists: $RUN_ROOT. Rename RUN_ROOT for a new run."


mkdir -p \
    "$RUN_ROOT" \
    "$FAA_OUT" \
    "$METADATA_OUT"


## ------------------------------------------------------------
## IMPORTANT:
##
## FeGenie writes helper files such as HMMlib.txt into the
## current working directory.
##
## This mirrors the previous successful 631-genome workflow.
## ------------------------------------------------------------

cd "$RUN_ROOT"


###############################################################
## Record software and input provenance ##
###############################################################

log \
    "Recording environment and input information."


{
    echo "run_started=$(date --iso-8601=seconds)"
    echo "hostname=$(hostname)"
    echo "slurm_job_id=${SLURM_JOB_ID:-not_slurm}"
    echo "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-not_set}"

    echo "input_root=$INPUT_ROOT"
    echo "input_count=$input_count"

    echo "fegenie_output=$FEGENIE_OUT"
    echo "protein_output=$FAA_OUT"

    echo "conda_environment=$CONDA_ENV"
    echo "conda_prefix=${CONDA_PREFIX:-not_set}"

    echo "iron_hmms=$iron_hmms"

} > "$METADATA_OUT/run_info.txt"


printf \
    '%s\n' \
    "${inputs[@]}" \
    > "$METADATA_OUT/input_files.txt"


conda list \
    > "$METADATA_OUT/conda_list.txt"


conda list --explicit \
    > "$METADATA_OUT/conda_environment_explicit.txt"


FeGenie.py -h \
    > "$METADATA_OUT/FeGenie_help.txt" \
    2>&1


FEGENIE_EXE="$(
    command -v FeGenie.py
)"


printf \
    '%s\n' \
    "$FEGENIE_EXE" \
    > "$METADATA_OUT/FeGenie_executable_path.txt"


sha256sum \
    "$FEGENIE_EXE" \
    > "$METADATA_OUT/FeGenie_executable.sha256"


find \
    "$iron_hmms" \
    -type f \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$METADATA_OUT/FeGenie_HMM_files.sha256"


###############################################################
## Stage the 187 MAG FASTAs into one flat directory ##
###############################################################

STAGE_PARENT="${SLURM_TMPDIR:-$RUN_ROOT}"

STAGE_DIR="$STAGE_PARENT/fegenie_inputs_${SLURM_JOB_ID:-$$}"


mkdir -p "$STAGE_DIR"


cleanup() {

    rm -rf "$STAGE_DIR"
}


trap cleanup EXIT


log \
    "Staging $input_count SemiBin2 MAG FASTAs into: $STAGE_DIR"


## ------------------------------------------------------------
## Keep original source paths so we can reconstruct provenance
## after FeGenie creates its protein FASTAs.
## ------------------------------------------------------------

declare -A ORIGINAL_BY_STAGED


for src in "${inputs[@]}"
do

    staged_name="$(
        basename "$src"
    )"


    dest="$STAGE_DIR/$staged_name"


    ## --------------------------------------------------------
    ## Every Genome_ID should be unique. Refuse collisions.
    ## --------------------------------------------------------

    if [[ -e "$dest" ]]; then

        fail \
            "Duplicate staged FASTA basename detected: $staged_name"

    fi


    cp -p \
        -- "$src" \
        "$dest.tmp"


    mv \
        -- "$dest.tmp" \
        "$dest"


    ORIGINAL_BY_STAGED[
        "$staged_name"
    ]="$src"


    [[ -s "$dest" ]] \
        || fail \
        "Staged FASTA is empty: $src"


    first_character="$(
        head -c 1 "$dest"
    )"


    [[ "$first_character" == ">" ]] \
        || fail \
        "File does not begin with a FASTA header: $src"


    ## --------------------------------------------------------
    ## FeGenie does not tolerate "|" in FASTA identifiers.
    ## --------------------------------------------------------

    bad_header="$(
        awk \
            '/^>/ && index($0, "|") {print; exit}' \
            "$dest"
    )"


    if [[ -n "$bad_header" ]]; then

        fail \
            "FeGenie cannot use FASTA headers containing '|'. File: $src; header: $bad_header"

    fi

done


staged=(
    "$STAGE_DIR"/*.fa
)


staged_count=${#staged[@]}


(( staged_count == input_count )) \
    || fail \
    "Staging count mismatch: $staged_count staged files versus $input_count inputs."


###############################################################
## Run the complete FeGenie workflow ##
###############################################################

threads=${SLURM_CPUS_PER_TASK:-16}


if (( threads > 16 )); then

    threads=16

fi


log \
    "Starting FeGenie with $threads threads."


## ------------------------------------------------------------
## These are reconstructed MAGs.
##
## Therefore:
##
##   --meta        intentionally omitted
##   --all_results intentionally omitted
##
## This mirrors the authoritative 631-genome workflow.
## ------------------------------------------------------------

FeGenie.py \
    -bin_dir "$STAGE_DIR" \
    -bin_ext fa \
    -out "$FEGENIE_OUT" \
    -t "$threads"


log \
    "FeGenie command finished."


###############################################################
## Validate the run ##
###############################################################

[[ -d "$FEGENIE_OUT/ORF_calls" ]] \
    || fail \
    "FeGenie did not create the expected ORF_calls directory."


[[ -s "$FEGENIE_OUT/FinalSummary.csv" ]] \
    || fail \
    "FeGenie did not create a non-empty FinalSummary.csv."


faa_files=(
    "$FEGENIE_OUT"/ORF_calls/*-proteins.faa
)


faa_count=${#faa_files[@]}


(( faa_count > 0 )) \
    || fail \
    "No FeGenie-predicted protein FASTA files were found."


###############################################################
## Collect FeGenie-predicted proteins ##
###############################################################

manifest="$METADATA_OUT/predicted_protein_manifest.tsv"


printf \
    'genome\toriginal_input\tfegenie_orf_file\tcollected_faa\n' \
    > "$manifest"


for faa in "${faa_files[@]}"
do

    base="$(
        basename "$faa"
    )"


    ## --------------------------------------------------------
    ## Example:
    ##
    ## semibin2__flye__barcode10_seqs__SemiBin_14.fa-proteins.faa
    ##
    ## becomes:
    ##
    ## semibin2__flye__barcode10_seqs__SemiBin_14
    ## --------------------------------------------------------

    staged_name="${base%-proteins.faa}"

    genome="${staged_name%.fa}"


    collected="$FAA_OUT/${genome}.faa"


    original="${ORIGINAL_BY_STAGED[$staged_name]:-UNKNOWN}"


    cp -p \
        -- "$faa" \
        "$collected"


    printf \
        '%s\t%s\t%s\t%s\n' \
        "$genome" \
        "$original" \
        "$faa" \
        "$collected" \
        >> "$manifest"

done


###############################################################
## Validate proteome count ##
###############################################################

if (( faa_count != input_count )); then

    fail \
        "Protein FASTA count mismatch: $faa_count predicted proteomes versus $input_count input MAGs. FeGenie outputs were retained for inspection."

fi


###############################################################
## Completion marker ##
###############################################################

{
    echo "run_finished=$(date --iso-8601=seconds)"

    echo "input_MAGs=$input_count"

    echo "predicted_proteomes=$faa_count"

    echo "fegenie_results=$FEGENIE_OUT"

    echo "collected_proteins=$FAA_OUT"

    echo "manifest=$manifest"

} > "$RUN_ROOT/RUN_COMPLETED.txt"


###############################################################
## Finish ##
###############################################################

log \
    "Run completed successfully."


log \
    "FeGenie results: $FEGENIE_OUT"


log \
    "Proteins for SignalP/FindMeHemes: $FAA_OUT"


log \
    "Run metadata: $METADATA_OUT"
