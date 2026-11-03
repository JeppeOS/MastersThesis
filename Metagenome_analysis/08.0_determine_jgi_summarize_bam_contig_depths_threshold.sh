#!/usr/bin/env bash

set -Eeuo pipefail


###############################################################################
## QC_barcode15_mapping_identity.sh
##
## Purpose
## -------
## Validate the alignment-identity threshold used by
## jgi_summarize_bam_contig_depths for long-read depth estimation.
##
## Barcode15 was used as a representative sample.
##
## The script:
##
##   1. Uses the barcode15 minimap2 map-ont BAM.
##   2. Generates JGI read statistics for all mappings.
##   3. Removes secondary and supplementary alignments.
##   4. Generates JGI read statistics for primary alignments only.
##   5. Summarizes the empirical alignment-identity distribution.
##   6. Compares retention at 90, 95, 97, 98 and 99% identity.
##
## The final pipeline retained 97% identity for contig-depth
## calculation.
###############################################################################


###############################################################################
## Settings
###############################################################################

PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

THREADS=8

ASSEMBLY="$PROJECT/flye_results/barcode15_seqs/assembly.fasta"

OUTDIR="$PROJECT/mapping_identity_qc"

mkdir -p "$OUTDIR"


###############################################################################
## Locate barcode15 BAM
##
## Different stages of the project used slightly different locations.
## Check the final mapping directory first, then the old root-level BAM.
###############################################################################

BAM_CANDIDATES=(

    "$PROJECT/mapping/barcode15_seqs/barcode15_seqs.mapont.sorted.bam"

    "$PROJECT/mapping/barcode15_mapont.sorted.bam"

    "$PROJECT/barcode15_mapont.sorted.bam"
)


BAM=""

for candidate in "${BAM_CANDIDATES[@]}"; do

    if [[ -f "$candidate" ]]; then

        BAM="$candidate"
        break

    fi

done


if [[ -z "$BAM" ]]; then

    echo "ERROR: Could not locate the barcode15 map-ont BAM." >&2

    printf 'Checked:\n' >&2

    printf '  %s\n' "${BAM_CANDIDATES[@]}" >&2

    exit 1
fi


###############################################################################
## Output files
###############################################################################

PRIMARY_BAM="$OUTDIR/barcode15_primary.bam"


FULL_READSTATS="$OUTDIR/barcode15_readstats.tsv"

FULL_DEPTH="$OUTDIR/barcode15_test_depth.tsv"

FULL_JGI_LOG="$OUTDIR/barcode15_jgi_depth.log"


PRIMARY_READSTATS="$OUTDIR/barcode15_primary_readstats.tsv"

PRIMARY_DEPTH="$OUTDIR/barcode15_primary_test_depth.tsv"

PRIMARY_JGI_LOG="$OUTDIR/barcode15_primary_jgi_depth.log"


SUMMARY="$OUTDIR/barcode15_mapping_identity_summary.tsv"

FLAGSTAT_FULL="$OUTDIR/barcode15_full.flagstat.txt"

FLAGSTAT_PRIMARY="$OUTDIR/barcode15_primary.flagstat.txt"


###############################################################################
## Helper
###############################################################################

log() {

    printf '[%s] %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$*"
}


###############################################################################
## Check software
###############################################################################

command -v samtools >/dev/null 2>&1 || {

    echo "ERROR: samtools not found." >&2
    exit 1
}


command -v jgi_summarize_bam_contig_depths >/dev/null 2>&1 || {

    echo "ERROR: jgi_summarize_bam_contig_depths not found." >&2
    exit 1
}


[[ -f "$ASSEMBLY" ]] || {

    echo "ERROR: Assembly not found: $ASSEMBLY" >&2
    exit 1
}


###############################################################################
## Report inputs
###############################################################################

log "Assembly:"
log "  $ASSEMBLY"

log "BAM:"
log "  $BAM"

log "Output directory:"
log "  $OUTDIR"


###############################################################################
## 1. Full BAM statistics
###############################################################################

log "============================================================"
log "1. Full barcode15 BAM"
log "============================================================"


samtools flagstat \
    -@ "$THREADS" \
    "$BAM" \
    > "$FLAGSTAT_FULL"


cat "$FLAGSTAT_FULL"


###############################################################################
## 2. JGI statistics for the full BAM
##
## --percentIdentity 97 is written explicitly here so the method does
## not depend on the program default.
###############################################################################

log "============================================================"
log "2. JGI statistics — all alignments"
log "============================================================"


jgi_summarize_bam_contig_depths \
    --referenceFasta "$ASSEMBLY" \
    --percentIdentity 97 \
    --outputReadStats "$FULL_READSTATS" \
    --outputDepth "$FULL_DEPTH" \
    "$BAM" \
    2> >(tee "$FULL_JGI_LOG" >&2)


###############################################################################
## 3. Construct primary-only BAM
##
## SAM flags:
##
##   0x100 = secondary alignment = 256
##   0x800 = supplementary alignment = 2048
##
## 256 + 2048 = 2304
##
## Therefore:
##
##   -F 2304
##
## excludes both secondary and supplementary alignments.
###############################################################################

log "============================================================"
log "3. Creating primary-only BAM"
log "============================================================"


samtools view \
    -@ "$THREADS" \
    -bh \
    -F 2304 \
    "$BAM" \
    > "$PRIMARY_BAM"


samtools index \
    -@ "$THREADS" \
    "$PRIMARY_BAM"


samtools flagstat \
    -@ "$THREADS" \
    "$PRIMARY_BAM" \
    > "$FLAGSTAT_PRIMARY"


cat "$FLAGSTAT_PRIMARY"


###############################################################################
## 4. JGI statistics for primary alignments
###############################################################################

log "============================================================"
log "4. JGI statistics — primary alignments only"
log "============================================================"


jgi_summarize_bam_contig_depths \
    --referenceFasta "$ASSEMBLY" \
    --percentIdentity 97 \
    --outputReadStats "$PRIMARY_READSTATS" \
    --outputDepth "$PRIMARY_DEPTH" \
    "$PRIMARY_BAM" \
    2> >(tee "$PRIMARY_JGI_LOG" >&2)


###############################################################################
## 5. Function for summarizing PctId
##
## We locate the PctId column from the header rather than assuming
## that it will always be column 4.
###############################################################################

summarize_identity() {

    local label="$1"
    local infile="$2"


    awk \
        -v dataset="$label" \
        'BEGIN {

            FS  = "\t"
            OFS = "\t"
        }

        NR == 1 {

            pct_col = 0

            for (i = 1; i <= NF; i++) {

                if ($i == "PctId") {

                    pct_col = i
                    break
                }
            }

            if (pct_col == 0) {

                print \
                    "ERROR: PctId column not found in " FILENAME \
                    > "/dev/stderr"

                exit 1
            }

            next
        }

        {
            id = $pct_col + 0

            n++
            sum += id

            if (id >= 0.90) n90++
            if (id >= 0.95) n95++
            if (id >= 0.97) n97++
            if (id >= 0.98) n98++
            if (id >= 0.99) n99++
        }

        END {

            if (n == 0) {

                print \
                    "ERROR: No alignments found in " FILENAME \
                    > "/dev/stderr"

                exit 1
            }

            print dataset, \
                  n, \
                  sum / n, \
                  n90 + 0, 100 * (n90 + 0) / n, \
                  n95 + 0, 100 * (n95 + 0) / n, \
                  n97 + 0, 100 * (n97 + 0) / n, \
                  n98 + 0, 100 * (n98 + 0) / n, \
                  n99 + 0, 100 * (n99 + 0) / n
        }' \
        "$infile"
}


###############################################################################
## 6. Build comparison table
###############################################################################

log "============================================================"
log "5. Calculating alignment-identity distributions"
log "============================================================"


printf \
'dataset\talignments\tmean_identity\tN_ge90\tpct_ge90\tN_ge95\tpct_ge95\tN_ge97\tpct_ge97\tN_ge98\tpct_ge98\tN_ge99\tpct_ge99\n' \
    > "$SUMMARY"


summarize_identity \
    "all_alignments" \
    "$FULL_READSTATS" \
    >> "$SUMMARY"


summarize_identity \
    "primary_only" \
    "$PRIMARY_READSTATS" \
    >> "$SUMMARY"


###############################################################################
## 7. Pretty-print results
###############################################################################

log "============================================================"
log "Identity summary"
log "============================================================"


awk '
BEGIN {
    FS = "\t"
}

NR == 1 {
    next
}

{
    printf "\n%s\n", $1
    printf "---------------------------------------------\n"

    printf "Alignments:       %d\n", $2
    printf "Mean identity:    %.2f%%\n", $3 * 100

    printf ">=90%% identity:   %d  (%.2f%%)\n", $4,  $5
    printf ">=95%% identity:   %d  (%.2f%%)\n", $6,  $7
    printf ">=97%% identity:   %d  (%.2f%%)\n", $8,  $9
    printf ">=98%% identity:   %d  (%.2f%%)\n", $10, $11
    printf ">=99%% identity:   %d  (%.2f%%)\n", $12, $13
}
' "$SUMMARY"


###############################################################################
## 8. Specific 95% vs 97% comparison for primary alignments
###############################################################################

log "============================================================"
log "Primary-alignment threshold comparison"
log "============================================================"


awk '
BEGIN {
    FS = "\t"
}

NR > 1 && $1 == "primary_only" {

    printf "Primary alignments: %d\n\n", $2

    printf "95%% threshold:\n"
    printf "  retained = %d\n", $6
    printf "  fraction = %.2f%%\n\n", $7

    printf "97%% threshold:\n"
    printf "  retained = %d\n", $8
    printf "  fraction = %.2f%%\n\n", $9

    printf "Additional alignments admitted at 95%% instead of 97%%:\n"
    printf "  %d\n", $6 - $8
}
' "$SUMMARY"


###############################################################################
## 9. Reference values from the original barcode15 check
##
## These are printed only as a comparison. They are NOT used to alter
## or filter the results.
###############################################################################

cat <<'EOF'

Expected values from the original barcode15 diagnostic:

Primary alignments:
    123,815

Mean PctId:
    92.07%

>=95%:
    62.06%

>=97%:
    53.37%

JGI readsWellMapped at 97%:
    66,082

The >=97% primary-alignment count and the JGI readsWellMapped count
should therefore be approximately identical.

EOF


###############################################################################
## Complete
###############################################################################

log "============================================================"
log "Barcode15 mapping-identity QC complete"
log "============================================================"

log "Summary:"
log "  $SUMMARY"

log "Full read statistics:"
log "  $FULL_READSTATS"

log "Primary read statistics:"
log "  $PRIMARY_READSTATS"

log "Full JGI log:"
log "  $FULL_JGI_LOG"

log "Primary JGI log:"
log "  $PRIMARY_JGI_LOG"

log "Primary BAM:"
log "  $PRIMARY_BAM"
