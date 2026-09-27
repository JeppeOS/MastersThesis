#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob


###############################################################################
## 01_run_sylph.sh
##
## Sylph abundance profiling of the Nanopore metagenomes.
##
## Part 1:
##   - Sketch Porechop-trimmed Nanopore reads
##   - Profile read sketches against the GTDB r226 Sylph database
##
## Part 2:
##   - Build a Sylph database from the recovered MAGs
##   - Query the MAG database against the read sketches
##
## Main outputs:
##
##   sylph_results/read_sketches/*.sylsp
##   sylph_results/gtdb_r226_profile.tsv
##   sylph_results/MAG_database.syldb
##   sylph_results/MAG_vs_reads.tsv
###############################################################################


###############################################################################
## User settings
###############################################################################

PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

REFERENCE_DIR="$HOME/methanotrophs/methanotroph_project/jeppe/reference_data/sylph"

SYLPHDB="$REFERENCE_DIR/gtdb-r226-c200-dbv2.syl2db"

READ_DIR="$PROJECT/porechop_results"

OUTDIR="$PROJECT/sylph_results"

READ_SKETCH_DIR="$OUTDIR/read_sketches"

MAG_LIST="$OUTDIR/MAG_list.txt"

MAG_DB_PREFIX="$OUTDIR/MAG_database"

THREADS=32


###############################################################################
## Helper functions
###############################################################################

log() {

    printf '[%s] %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$*"
}


fail() {

    printf '[%s] ERROR: %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$*" >&2

    exit 1
}


###############################################################################
## Activate Sylph Conda environment
###############################################################################

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

    fail "Conda could not be found."

fi


conda activate sylph


command -v sylph >/dev/null 2>&1 \
    || fail "Sylph is not available in the active Conda environment."


log "Sylph version:"
sylph --version || true


###############################################################################
## Create output directories
###############################################################################

mkdir -p \
    "$OUTDIR" \
    "$READ_SKETCH_DIR"


###############################################################################
## Locate trimmed Nanopore reads
###############################################################################

READS=(
    "$READ_DIR"/*_trimmed.fastq.gz
)


if (( ${#READS[@]} == 0 )); then

    fail \
        "No Porechop-trimmed FASTQ files found in:
         $READ_DIR"
fi


log "Trimmed read files found: ${#READS[@]}"

printf '  %s\n' "${READS[@]}"


###############################################################################
## Check GTDB r226 Sylph database
###############################################################################

[[ -f "$SYLPHDB" ]] || \
    fail "GTDB r226 Sylph database not found: $SYLPHDB"


log "GTDB r226 Sylph database:"
log "  $SYLPHDB"


###############################################################################
## STEP 1 — Sketch Nanopore reads
##
## This is the original command used in the analysis:
##
## sylph sketch \
##     -r porechop_results/*_trimmed.fastq.gz \
##     -d sylph_results/read_sketches \
##     -t 32
###############################################################################

log "============================================================"
log "STEP 1 — Sketching trimmed Nanopore reads"
log "============================================================"


sylph sketch \
    -r "${READS[@]}" \
    -d "$READ_SKETCH_DIR" \
    -t "$THREADS"


###############################################################################
## Check read sketches
###############################################################################

READ_SKETCHES=(
    "$READ_SKETCH_DIR"/*.sylsp
)


if (( ${#READ_SKETCHES[@]} == 0 )); then

    fail "Sylph did not generate any .sylsp read sketches."
fi


log "Read sketches generated: ${#READ_SKETCHES[@]}"


###############################################################################
## STEP 2 — Profile reads against GTDB r226
##
## Original command:
##
## sylph profile \
##     "$SYLPHDB" \
##     sylph_results/read_sketches/*.sylsp \
##     -t 32 \
##     -o sylph_results/gtdb_r226_profile.tsv
###############################################################################

log "============================================================"
log "STEP 2 — Profiling reads against GTDB r226"
log "============================================================"


sylph profile \
    "$SYLPHDB" \
    "${READ_SKETCHES[@]}" \
    -t "$THREADS" \
    -o "$OUTDIR/gtdb_r226_profile.tsv"


[[ -s "$OUTDIR/gtdb_r226_profile.tsv" ]] || \
    fail "GTDB r226 Sylph profile was not generated."


log "GTDB profile written to:"
log "  $OUTDIR/gtdb_r226_profile.tsv"


###############################################################################
## STEP 3 — Build Sylph database from recovered MAGs
##
## MAG_list.txt contains one MAG FASTA path per line.
##
## Original command:
##
## sylph sketch \
##     -l sylph_results/MAG_list.txt \
##     -o sylph_results/MAG_database \
##     -t 32
###############################################################################

log "============================================================"
log "STEP 3 — Building MAG Sylph database"
log "============================================================"


[[ -s "$MAG_LIST" ]] || \
    fail "MAG list not found or empty: $MAG_LIST"


N_MAGS=$(
    grep -cvE '^[[:space:]]*$' "$MAG_LIST"
)


log "MAGs listed: $N_MAGS"


## Check that every FASTA listed in MAG_list.txt exists.
while IFS= read -r mag; do

    [[ -n "$mag" ]] || continue

    [[ -f "$mag" ]] || \
        fail "MAG FASTA listed in MAG_list.txt does not exist: $mag"

done < "$MAG_LIST"


sylph sketch \
    -l "$MAG_LIST" \
    -o "$MAG_DB_PREFIX" \
    -t "$THREADS"


[[ -s "${MAG_DB_PREFIX}.syldb" ]] || \
    fail "MAG Sylph database was not generated."


log "MAG Sylph database written to:"
log "  ${MAG_DB_PREFIX}.syldb"


###############################################################################
## STEP 4 — Query MAGs against the metagenome read sketches
##
## Original command:
##
## sylph query \
##     MAG_database.syldb \
##     read_sketches/*.sylsp \
##     -t 32 \
##     -o MAG_vs_reads.tsv
###############################################################################

log "============================================================"
log "STEP 4 — Querying MAGs against read sketches"
log "============================================================"


sylph query \
    "${MAG_DB_PREFIX}.syldb" \
    "${READ_SKETCHES[@]}" \
    -t "$THREADS" \
    -o "$OUTDIR/MAG_vs_reads.tsv"


[[ -s "$OUTDIR/MAG_vs_reads.tsv" ]] || \
    fail "MAG_vs_reads.tsv was not generated."


###############################################################################
## Final report
###############################################################################

log "============================================================"
log "Sylph analysis completed successfully"
log "============================================================"

log "Read sketches:"
log "  $READ_SKETCH_DIR"

log "GTDB r226 profile:"
log "  $OUTDIR/gtdb_r226_profile.tsv"

log "MAG database:"
log "  ${MAG_DB_PREFIX}.syldb"

log "MAG vs reads:"
log "  $OUTDIR/MAG_vs_reads.tsv"
