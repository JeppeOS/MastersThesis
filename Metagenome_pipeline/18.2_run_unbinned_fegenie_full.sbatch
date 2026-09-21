#!/usr/bin/env bash

#SBATCH --job-name=fegenie_unbinned
#SBATCH --partition=normal
#SBATCH --account=methanotrophs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fegenie_unbinned_%j.out
#SBATCH --error=logs/fegenie_unbinned_%j.err


set -Eeo pipefail
shopt -s nullglob
umask 0027


###############################################################
## User settings ##
###############################################################

BASE_DIR="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

INPUT_ROOT="$BASE_DIR/functional_analysis/unbinned_contigs/input"


## ------------------------------------------------------------
## Independent FeGenie run for unbinned metagenomic contigs.
##
## Change this name for a later intentional rerun.
## ------------------------------------------------------------

RUN_ROOT="$BASE_DIR/functional_analysis/unbinned_contigs/FeGenie_unbinned_full_run_20260918"

FEGENIE_OUT="$RUN_ROOT/fegenie_results"

FAA_OUT="$RUN_ROOT/predicted_proteins_faa"

METADATA_OUT="$RUN_ROOT/run_metadata"


CONDA_ENV="fegenie_conda"

EXPECTED_INPUT_COUNT=7


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
## Expected structure:
##
## input/
## ├── barcode10_seqs/
## │   └── barcode10_seqs__unbinned_ge1000.fa
## ├── barcode11_seqs/
## │   └── barcode11_seqs__unbinned_ge1000.fa
## ...
## └── barcode16_seqs/
##     └── barcode16_seqs__unbinned_ge1000.fa
## ------------------------------------------------------------

inputs=(
    "$INPUT_ROOT"/*/*__unbinned_ge1000.fa
)


input_count=${#inputs[@]}


(( input_count > 0 )) \
    || fail \
    "No unbinned FASTAs matching $INPUT_ROOT/*/*__unbinned_ge1000.fa were found."


if (( input_count != EXPECTED_INPUT_COUNT )); then

    log \
        "WARNING: found $input_count input files; expected $EXPECTED_INPUT_COUNT."

    log \
        "The run will continue, but check for missing or stray samples."

fi


###############################################################
## FeGenie HMM environment ##
###############################################################

[[ -n "${iron_hmms:-}" ]] \
    || fail \
    "The FeGenie environment variable iron_hmms is not set after Conda activation."


[[ -d "$iron_hmms" ]] \
    || fail \
    "The FeGenie HMM directory does not exist: $iron_hmms"


[[ -f "$iron_hmms/HMM-bitcutoffs.txt" ]] \
    || fail \
    "HMM-bitcutoffs.txt is missing from: $iron_hmms"


###############################################################
## Refuse to overwrite an earlier run ##
###############################################################

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
## We intentionally reproduce the working-directory behaviour
## of the previously successful 631-genome workflow.
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

    echo "fegenie_meta_mode=yes"

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
## Stage the seven unbinned FASTAs into one flat directory ##
###############################################################

STAGE_PARENT="${SLURM_TMPDIR:-$RUN_ROOT}"

STAGE_DIR="$STAGE_PARENT/fegenie_unbinned_inputs_${SLURM_JOB_ID:-$$}"


mkdir -p "$STAGE_DIR"


cleanup() {

    rm -rf "$STAGE_DIR"

}


trap cleanup EXIT


log \
    "Staging $input_count unbinned metagenomic FASTAs into: $STAGE_DIR"


## ------------------------------------------------------------
## Store source path for each staged FASTA.
## ------------------------------------------------------------

declare -A ORIGINAL_BY_STAGED


for src in "${inputs[@]}"
do

    staged_name="$(
        basename "$src"
    )"


    dest="$STAGE_DIR/$staged_name"


    ## --------------------------------------------------------
    ## Refuse basename collisions.
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


    ## --------------------------------------------------------
    ## Basic FASTA validation.
    ## --------------------------------------------------------

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
## Record staged input statistics ##
###############################################################

input_stats="$METADATA_OUT/unbinned_input_statistics.tsv"


printf \
    'dataset\toriginal_input\tcontigs\tbp\n' \
    > "$input_stats"


for staged_fasta in "${staged[@]}"
do

    staged_name="$(
        basename "$staged_fasta"
    )"


    dataset="${staged_name%.fa}"


    original="${ORIGINAL_BY_STAGED[$staged_name]}"


    n_contigs="$(
        grep -c '^>' "$staged_fasta" \
        || true
    )"


    total_bp="$(
        awk '
            /^>/ {
                next
            }

            {
                n += length($0)
            }

            END {
                print n + 0
            }
        ' "$staged_fasta"
    )"


    printf \
        '%s\t%s\t%s\t%s\n' \
        "$dataset" \
        "$original" \
        "$n_contigs" \
        "$total_bp" \
        >> "$input_stats"

done


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
## IMPORTANT DIFFERENCE FROM THE MAG RUN:
##
## These are NOT reconstructed genomes.
##
## Each FASTA contains the unbinned metagenomic contigs from
## one barcode.
##
## Therefore:
##
##     --meta        INCLUDED
##     --all_results OMITTED
##
## The latter keeps the primary reporting rules consistent
## with the MAG FeGenie analysis.
## ------------------------------------------------------------

FeGenie.py \
    -bin_dir "$STAGE_DIR" \
    -bin_ext fa \
    -out "$FEGENIE_OUT" \
    -t "$threads" \
    --meta


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
    'dataset\toriginal_input\tfegenie_orf_file\tcollected_faa\n' \
    > "$manifest"


for faa in "${faa_files[@]}"
do

    base="$(
        basename "$faa"
    )"


    ## --------------------------------------------------------
    ## Example:
    ##
    ## barcode16_seqs__unbinned_ge1000.fa-proteins.faa
    ##
    ## becomes:
    ##
    ## barcode16_seqs__unbinned_ge1000
    ## --------------------------------------------------------

    staged_name="${base%-proteins.faa}"

    dataset="${staged_name%.fa}"


    collected="$FAA_OUT/${dataset}.faa"


    original="${ORIGINAL_BY_STAGED[$staged_name]:-UNKNOWN}"


    cp -p \
        -- "$faa" \
        "$collected"


    printf \
        '%s\t%s\t%s\t%s\n' \
        "$dataset" \
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
        "Protein FASTA count mismatch: $faa_count predicted proteomes versus $input_count input metagenomes. FeGenie outputs were retained for inspection."

fi


###############################################################
## Completion marker ##
###############################################################

{
    echo "run_finished=$(date --iso-8601=seconds)"

    echo "input_metagenomes=$input_count"

    echo "predicted_proteomes=$faa_count"

    echo "meta_mode=yes"

    echo "fegenie_results=$FEGENIE_OUT"

    echo "collected_proteins=$FAA_OUT"

    echo "manifest=$manifest"

    echo "input_statistics=$input_stats"

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
