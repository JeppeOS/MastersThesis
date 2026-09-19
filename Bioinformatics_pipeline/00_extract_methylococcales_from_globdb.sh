#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob
umask 0027


###############################################################################
## Stage 00
##
## Extract the complete Methylococcales genome set from GLOBDB r226.
##
## Selection criterion:
##     GTDB/GLOBDB taxonomy order == o__Methylococcales
##
## Expected result:
##     631 unique Methylococcales genomes
##
## GLOBDB genome structure:
##     globdb_r226_genome_fasta/
##         R226CHUNK000/
##         R226CHUNK001/
##         ...
##         R226CHUNK030/
##
## Genome files:
##     <genome_id>.fa.gz
##
## Final output:
##     methylococcales_fastas/
###############################################################################


###############################################################################
## User settings
###############################################################################

BASE_DIR="$HOME/methanotrophs/methanotroph_project/jeppe"

GLOBDB_DIR="$BASE_DIR/globdb"

## Extracted GLOBDB r226 genome FASTA collection.
GENOME_SOURCE_DIR="$GLOBDB_DIR/globdb_r226_genome_fasta"

## The taxonomy was originally distributed compressed.
TAXONOMY_GZ="$GLOBDB_DIR/globdb_r226_taxonomy.tsv.gz"

## Also allow an already-uncompressed copy.
TAXONOMY_TSV="$GLOBDB_DIR/globdb_r226_taxonomy.tsv"

## Authoritative genome collection used by Stage 01 onward.
OUT_DIR="$BASE_DIR/methylococcales_fastas"

## Keep Stage-00 manifests/QC with the genome-analysis workflow.
WORKFLOW_DIR="$BASE_DIR/genome_analysis_workflow"
STAGE_DIR="$WORKFLOW_DIR/00_methylococcales_input"

EXPECTED_GENOMES=631


###############################################################################
## Output files
###############################################################################

IDS_FILE="$STAGE_DIR/methylococcales_ids.txt"

TAXONOMY_OUT="$STAGE_DIR/methylococcales_taxonomy.tsv"

SOURCE_MANIFEST="$STAGE_DIR/methylococcales_source_manifest.tsv"

OBSERVED_IDS="$STAGE_DIR/methylococcales_extracted_ids.txt"

MISSING_IDS="$STAGE_DIR/missing_genomes.txt"

EXTRA_IDS="$STAGE_DIR/unexpected_genomes.txt"

QC_FILE="$STAGE_DIR/extraction_qc.tsv"


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
## Check input files
###############################################################################

[[ -d "$GENOME_SOURCE_DIR" ]] || \
    fail "GLOBDB genome directory not found: $GENOME_SOURCE_DIR"


if [[ -f "$TAXONOMY_GZ" ]]; then

    TAXONOMY_FILE="$TAXONOMY_GZ"
    TAXONOMY_COMMAND=(gzip -cd "$TAXONOMY_FILE")

elif [[ -f "$TAXONOMY_TSV" ]]; then

    TAXONOMY_FILE="$TAXONOMY_TSV"
    TAXONOMY_COMMAND=(cat "$TAXONOMY_FILE")

else

    fail \
        "Could not find either:
         $TAXONOMY_GZ
         $TAXONOMY_TSV"
fi


mkdir -p \
    "$STAGE_DIR" \
    "$OUT_DIR"


log "GLOBDB genome directory:"
log "  $GENOME_SOURCE_DIR"

log "Taxonomy:"
log "  $TAXONOMY_FILE"

log "Output genome directory:"
log "  $OUT_DIR"


###############################################################################
## Extract Methylococcales taxonomy rows
##
## GLOBDB taxonomy format:
##
## genome_id<TAB>d__...;p__...;c__...;o__...;f__...;g__...;s__...
##
## We explicitly split the taxonomy string and require the exact rank:
##
##     o__Methylococcales
##
## rather than using a loose substring search.
###############################################################################

log "Selecting genomes assigned to o__Methylococcales..."


"${TAXONOMY_COMMAND[@]}" \
    | awk -F '\t' '
        {
            n = split($2, ranks, ";")

            for (i = 1; i <= n; i++) {

                if (ranks[i] == "o__Methylococcales") {

                    print $0
                    break
                }
            }
        }
    ' \
    | LC_ALL=C sort -t $'\t' -k1,1 \
    > "$TAXONOMY_OUT"


awk -F '\t' '
    NF >= 1 && $1 != "" {
        print $1
    }
' "$TAXONOMY_OUT" \
    | LC_ALL=C sort -u \
    > "$IDS_FILE"


###############################################################################
## Validate taxonomy selection
###############################################################################

N_TAXONOMY_ROWS=$(
    wc -l < "$TAXONOMY_OUT"
)

N_IDS=$(
    wc -l < "$IDS_FILE"
)


log "Methylococcales taxonomy rows: $N_TAXONOMY_ROWS"
log "Unique Methylococcales genome IDs: $N_IDS"


if (( N_IDS != EXPECTED_GENOMES )); then

    fail \
        "Expected $EXPECTED_GENOMES unique Methylococcales genomes, \
but taxonomy selection returned $N_IDS."
fi


###############################################################################
## Load requested genome IDs
###############################################################################

declare -A WANTED
declare -A FOUND


while IFS= read -r genome; do

    [[ -n "$genome" ]] || continue

    WANTED["$genome"]=1

done < "$IDS_FILE"


###############################################################################
## Search GLOBDB chunks once and copy selected genome FASTAs
###############################################################################

log "Scanning GLOBDB r226 genome chunks..."

printf 'genome\tsource_file\toutput_file\n' \
    > "$SOURCE_MANIFEST"


N_FOUND=0


while IFS= read -r -d '' source_file; do

    filename="$(
        basename "$source_file"
    )"

    genome="${
        filename%.fa.gz
    }"

    ## Ignore genomes outside Methylococcales.
    if [[ -z "${WANTED[$genome]+x}" ]]; then
        continue
    fi


    ## A genome ID should occur exactly once in the source database.
    if [[ -n "${FOUND[$genome]+x}" ]]; then

        fail \
            "Genome $genome occurs more than once in the GLOBDB FASTA collection:
             first:  ${FOUND[$genome]}
             second: $source_file"
    fi


    destination="$OUT_DIR/$filename"


    cp -p \
        "$source_file" \
        "$destination"


    FOUND["$genome"]="$source_file"

    N_FOUND=$((N_FOUND + 1))


    printf '%s\t%s\t%s\n' \
        "$genome" \
        "$source_file" \
        "$destination" \
        >> "$SOURCE_MANIFEST"

done < <(
    find "$GENOME_SOURCE_DIR" \
        -type f \
        -name '*.fa.gz' \
        -print0
)


log "Matching genome FASTAs copied: $N_FOUND"


###############################################################################
## Build list of extracted genome IDs
###############################################################################

find "$OUT_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.fa.gz' \
    -printf '%f\n' \
    | sed 's/\.fa\.gz$//' \
    | LC_ALL=C sort -u \
    > "$OBSERVED_IDS"


N_OUTPUT=$(
    wc -l < "$OBSERVED_IDS"
)


###############################################################################
## Exact expected-vs-observed comparison
###############################################################################

comm -23 \
    "$IDS_FILE" \
    "$OBSERVED_IDS" \
    > "$MISSING_IDS"


comm -13 \
    "$IDS_FILE" \
    "$OBSERVED_IDS" \
    > "$EXTRA_IDS"


N_MISSING=$(
    wc -l < "$MISSING_IDS"
)

N_EXTRA=$(
    wc -l < "$EXTRA_IDS"
)


###############################################################################
## Check for unrelated files in the output directory
##
## Historically this caught screenlog.0, which made a simple `ls | wc -l`
## report 632 entries even though there were correctly 631 genome FASTAs.
###############################################################################

NON_FASTA_COUNT=$(
    find "$OUT_DIR" \
        -maxdepth 1 \
        -type f \
        ! -name '*.fa.gz' \
        | wc -l
)


###############################################################################
## Write QC summary
###############################################################################

{
    printf 'metric\tvalue\n'

    printf 'target_taxon\t%s\n' \
        'o__Methylococcales'

    printf 'expected_genomes\t%d\n' \
        "$EXPECTED_GENOMES"

    printf 'taxonomy_rows\t%d\n' \
        "$N_TAXONOMY_ROWS"

    printf 'unique_taxonomy_ids\t%d\n' \
        "$N_IDS"

    printf 'matched_source_fastas\t%d\n' \
        "$N_FOUND"

    printf 'output_fastas\t%d\n' \
        "$N_OUTPUT"

    printf 'missing_genomes\t%d\n' \
        "$N_MISSING"

    printf 'unexpected_genomes\t%d\n' \
        "$N_EXTRA"

    printf 'non_fasta_files_in_output_directory\t%d\n' \
        "$NON_FASTA_COUNT"

} > "$QC_FILE"


###############################################################################
## Final validation
###############################################################################

if (( N_MISSING > 0 )); then

    printf '\nMissing genomes:\n' >&2
    cat "$MISSING_IDS" >&2

    fail \
        "$N_MISSING expected Methylococcales genomes were not found."
fi


if (( N_EXTRA > 0 )); then

    printf '\nUnexpected genome FASTAs:\n' >&2
    cat "$EXTRA_IDS" >&2

    fail \
        "$N_EXTRA unexpected .fa.gz genome files are present in $OUT_DIR."
fi


if (( N_OUTPUT != EXPECTED_GENOMES )); then

    fail \
        "Final output contains $N_OUTPUT genome FASTAs; \
expected $EXPECTED_GENOMES."
fi


###############################################################################
## Report unrelated files without counting them as genomes
###############################################################################

if (( NON_FASTA_COUNT > 0 )); then

    log \
        "WARNING: $NON_FASTA_COUNT non-.fa.gz file(s) are also present in \
$OUT_DIR."

    find "$OUT_DIR" \
        -maxdepth 1 \
        -type f \
        ! -name '*.fa.gz' \
        -printf '  %f\n'
fi


###############################################################################
## Success
###############################################################################

log "Stage 00 completed successfully."

log \
    "Recovered $N_OUTPUT / $EXPECTED_GENOMES \
Methylococcales genome FASTAs."

log "Genome IDs:"
log "  $IDS_FILE"

log "Selected taxonomy:"
log "  $TAXONOMY_OUT"

log "Source manifest:"
log "  $SOURCE_MANIFEST"

log "QC:"
log "  $QC_FILE"

log "Final genome FASTAs:"
log "  $OUT_DIR"
