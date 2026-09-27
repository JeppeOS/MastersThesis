#!/bin/bash

set -euo pipefail


## ================================================================== ##
## STAGE 93E
##
## FIX MMO CLASSIFICATION
##
## BUG:
##
## Stage 93 used MMseqs2 "fident", which is fractional:
##
##     0.35 = 35% identity
##     0.70 = 70% identity
##     1.00 = 100% identity
##
## but the downstream Python filter incorrectly used:
##
##     MIN_IDENTITY = 35.0
##
## This rejected every hit.
##
## This script:
##
##   1. permanently fixes Stage 93 to MIN_IDENTITY = 0.35
##   2. reuses the EXISTING raw MMseqs output
##   3. rebuilds the best-hit table
##   4. rebuilds MMO status for all 650 genomes
##   5. rebuilds the final extra tree-track table
##
## No MMseqs search is rerun.
## ================================================================== ##


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

PHYLO="$PROJECT/comparative_analysis/final_species_phylogeny"

OUTDIR="$PHYLO/08_visualization/extra_tracks"


STAGE93="scripts/93_build_cluster00121_mmo_tree_tracks.sh"


META="$PHYLO/08_visualization/final_650_tip_metadata.tsv"


RAW_HITS="$OUTDIR/mmseqs/MMO_all650_raw_hits.tsv"

BEST_HITS="$OUTDIR/mmseqs/MMO_all650_best_component_hits.tsv"

MMO_SUMMARY="$OUTDIR/all_650_MMO_status.tsv"

FINAL_TRACKS="$OUTDIR/final_extra_tree_tracks.tsv"


C121_FOCAL="$PROJECT/comparative_analysis/Cluster_00121/Cluster_00121_focal_regions.tsv"

C121_COMBINED="$PROJECT/comparative_analysis/Cluster_00121/Cluster_00121_combined_neighborhood_genes_local_families.tsv"


## ================================================================== ##
## QC
## ================================================================== ##

for F in \
    "$STAGE93" \
    "$META" \
    "$RAW_HITS"

do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: required file missing or empty:"
        echo "$F"
        exit 1
    fi

done


echo
echo "Raw MMseqs hits:"
wc -l "$RAW_HITS"

echo
echo "First five raw hits:"
head -n 5 "$RAW_HITS"


## ================================================================== ##
## PART A
##
## Permanently repair Stage 93.
## ================================================================== ##

python - \
    "$STAGE93" \
<<'PY'
import sys
from pathlib import Path


path = Path(
    sys.argv[1]
)

text = path.read_text()


old = "MIN_IDENTITY = 35.0"

new = "MIN_IDENTITY = 0.35"


if old in text:

    backup = Path(
        str(path)
        +
        ".pre_fident_fix"
    )


    if not backup.exists():

        backup.write_text(
            text
        )


    text = text.replace(
        old,
        new,
        1
    )


    path.write_text(
        text
    )


    print()
    print(
        "Fixed Stage 93:"
    )

    print(
        "  MIN_IDENTITY = 35.0"
    )

    print(
        "  ->"
    )

    print(
        "  MIN_IDENTITY = 0.35"
    )


elif new in text:

    print()
    print(
        "Stage 93 already contains "
        "MIN_IDENTITY = 0.35"
    )


else:

    raise SystemExit(
        "ERROR: could not locate MIN_IDENTITY setting "
        "in Stage 93."
    )
PY


## ================================================================== ##
## PART B
##
## Reclassify EXISTING raw MMseqs hits.
## ================================================================== ##

python - \
    "$META" \
    "$RAW_HITS" \
    "$C121_FOCAL" \
    "$C121_COMBINED" \
    "$BEST_HITS" \
    "$MMO_SUMMARY" \
    "$FINAL_TRACKS" \
<<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path


meta_file = Path(
    sys.argv[1]
)

hits_file = Path(
    sys.argv[2]
)

c121_focal = Path(
    sys.argv[3]
)

c121_combined = Path(
    sys.argv[4]
)

best_out = Path(
    sys.argv[5]
)

mmo_out = Path(
    sys.argv[6]
)

tracks_out = Path(
    sys.argv[7]
)


## ================================================================== ##
## Correct MMseqs thresholds
##
## fident, qcov and tcov are fractions.
## ================================================================== ##

MIN_IDENTITY = 0.35

MIN_QCOV = 0.70

MIN_TCOV = 0.70


VALID_COMPONENTS = {
    "pmoA",
    "pmoB",
    "pmoC",
    "mmoX",
    "mmoY",
    "mmoZ",
}


## ================================================================== ##
## MAG ID crosswalk for Cluster_00121
## ================================================================== ##

crosswalk = {

    "semibin2__flye__barcode10_seqs__SemiBin_14":
        "MAG_b10_SB14",

    "semibin2__flye__barcode10_seqs__SemiBin_1":
        "MAG_b10_SB1",

    "semibin2__flye__barcode10_seqs__SemiBin_8":
        "MAG_b10_SB8",

    "semibin2__flye__barcode12_seqs__SemiBin_5":
        "MAG_b12_SB5",

    "semibin2__flye__barcode12_seqs__SemiBin_0":
        "MAG_b12_SB0",

    "semibin2__flye__barcode12_seqs__SemiBin_6":
        "MAG_b12_SB6",

    "semibin2__flye__barcode13_seqs__SemiBin_1":
        "MAG_b13_SB1",

    "semibin2__flye__barcode14_seqs__SemiBin_1":
        "MAG_b14_SB1",

    "semibin2__flye__barcode14_seqs__SemiBin_0":
        "MAG_b14_SB0",

    "semibin2__flye__barcode15_seqs__SemiBin_0":
        "MAG_b15_SB0",

    "semibin2__flye__barcode15_seqs__SemiBin_11":
        "MAG_b15_SB11",

    "semibin2__flye__barcode15_seqs__SemiBin_3":
        "MAG_b15_SB3",

    "semibin2__flye__barcode15_seqs__SemiBin_10":
        "MAG_b15_SB10",

    "semibin2__flye__barcode16_seqs__SemiBin_1":
        "MAG_b16_SB1",

    "semibin2__flye__barcode16_seqs__SemiBin_12":
        "MAG_b16_SB12",

    "semibin2__flye__barcode16_seqs__SemiBin_11":
        "MAG_b16_SB11",
}


## ================================================================== ##
## Read tree metadata
## ================================================================== ##

with meta_file.open() as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )


    meta = list(
        reader
    )


tree_genomes = [
    row["genome"]
    for row in meta
]


if len(
    tree_genomes
) != 650:

    raise SystemExit(
        f"ERROR: expected 650 tree genomes; "
        f"found {len(tree_genomes)}."
    )


if len(
    set(
        tree_genomes
    )
) != 650:

    raise SystemExit(
        "ERROR: duplicate genome IDs in tree metadata."
    )


tree_genome_set = set(
    tree_genomes
)


## ================================================================== ##
## Parse raw MMseqs results
##
## Stage 93 wrote:
##
## query
## target
## fident
## alnlen
## qcov
## tcov
## evalue
## bits
## ================================================================== ##

best_by_query = {}


raw_count = 0

threshold_pass_count = 0


with hits_file.open() as handle:

    reader = csv.reader(
        handle,
        delimiter="\t"
    )


    for row in reader:

        if len(row) != 8:
            continue


        raw_count += 1


        (
            query,
            target,
            fident,
            alnlen,
            qcov,
            tcov,
            evalue,
            bits,
        ) = row


        try:

            fident_f = float(
                fident
            )

            qcov_f = float(
                qcov
            )

            tcov_f = float(
                tcov
            )

            bits_f = float(
                bits
            )

        except ValueError:

            continue


        ## ---------------------------------------------------------- ##
        ## Correct fractional thresholds
        ## ---------------------------------------------------------- ##

        if fident_f < MIN_IDENTITY:
            continue


        if qcov_f < MIN_QCOV:
            continue


        if tcov_f < MIN_TCOV:
            continue


        threshold_pass_count += 1


        if "|" not in query:

            continue


        if "|" not in target:

            continue


        genome, protein = query.split(
            "|",
            1
        )


        component, accession = target.split(
            "|",
            1
        )


        if genome not in tree_genome_set:

            continue


        if component not in VALID_COMPONENTS:

            continue


        ## ---------------------------------------------------------- ##
        ## One protein can hit several references.
        ##
        ## Keep its best bitscore assignment.
        ## ---------------------------------------------------------- ##

        current = best_by_query.get(
            query
        )


        if (
            current is None
            or
            bits_f >
            current["bits"]
        ):

            best_by_query[
                query
            ] = {

                "genome":
                    genome,

                "protein":
                    protein,

                "component":
                    component,

                "accession":
                    accession,

                "fident":
                    fident_f,

                "qcov":
                    qcov_f,

                "tcov":
                    tcov_f,

                "evalue":
                    evalue,

                "bits":
                    bits_f,
            }


print()
print(
    "Raw MMseqs alignments:",
    raw_count
)

print(
    "Alignments passing 35% identity + 70/70% coverage:",
    threshold_pass_count
)

print(
    "Best retained query proteins:",
    len(
        best_by_query
    )
)


if len(
    best_by_query
) == 0:

    raise SystemExit(
        "\nERROR: still no passing MMO hits after correcting "
        "fident to fractional units.\n"
        "Do not plot the MMO track; inspect raw MMseqs output."
    )


## ================================================================== ##
## Write best hit table
## ================================================================== ##

best_fields = [

    "genome",
    "protein",
    "component",
    "accession",
    "fident",
    "qcov",
    "tcov",
    "evalue",
    "bits",
]


with best_out.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=best_fields
    )


    writer.writeheader()


    for query in sorted(
        best_by_query
    ):

        writer.writerow(
            best_by_query[
                query
            ]
        )


## ================================================================== ##
## Components per genome
## ================================================================== ##

components = defaultdict(
    set
)


for hit in best_by_query.values():

    components[
        hit["genome"]
    ].add(
        hit["component"]
    )


## ================================================================== ##
## Genome-level MMO classification
## ================================================================== ##

def classify(comp):

    pmmo = (
        "pmoA" in comp
        and
        (
            "pmoB" in comp
            or
            "pmoC" in comp
        )
    )


    smmo = (
        "mmoX" in comp
        and
        (
            "mmoY" in comp
            or
            "mmoZ" in comp
        )
    )


    if pmmo and smmo:
        return "Both"


    if pmmo:
        return "pMMO"


    if smmo:
        return "sMMO"


    return "None"


mmo_rows = []


for genome in tree_genomes:

    comp = components[
        genome
    ]


    status = classify(
        comp
    )


    partial = (
        bool(
            comp
        )
        and
        status ==
            "None"
    )


    mmo_rows.append(
        {

            "genome":
                genome,

            "pmoA":
                int(
                    "pmoA" in comp
                ),

            "pmoB":
                int(
                    "pmoB" in comp
                ),

            "pmoC":
                int(
                    "pmoC" in comp
                ),

            "mmoX":
                int(
                    "mmoX" in comp
                ),

            "mmoY":
                int(
                    "mmoY" in comp
                ),

            "mmoZ":
                int(
                    "mmoZ" in comp
                ),

            "MMO_status":
                status,

            "partial_component_hits":
                int(
                    partial
                ),

            "detected_components":
                ",".join(
                    sorted(
                        comp
                    )
                ),
        }
    )


mmo_fields = [

    "genome",
    "pmoA",
    "pmoB",
    "pmoC",
    "mmoX",
    "mmoY",
    "mmoZ",
    "MMO_status",
    "partial_component_hits",
    "detected_components",
]


with mmo_out.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=mmo_fields
    )


    writer.writeheader()

    writer.writerows(
        mmo_rows
    )


## ================================================================== ##
## Cluster_00121 positives
## ================================================================== ##

c121_positive = set()


source = (
    c121_focal
    if c121_focal.exists()
    else c121_combined
)


with source.open() as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )


    if "genome" not in (
        reader.fieldnames
        or []
    ):

        raise SystemExit(
            "ERROR: Cluster_00121 table lacks genome column."
        )


    for row in reader:

        if source == c121_combined:

            focal = str(
                row.get(
                    "is_focal",
                    ""
                )
            ).lower()


            if focal not in {
                "1",
                "true",
                "t",
                "yes",
            }:

                continue


        genome = (
            row.get(
                "genome",
                ""
            )
            .strip()
        )


        if not genome:

            continue


        genome = crosswalk.get(
            genome,
            genome
        )


        if genome in tree_genome_set:

            c121_positive.add(
                genome
            )


## ================================================================== ##
## Final plotting table
## ================================================================== ##

mmo_by_genome = {

    row["genome"]:
        row["MMO_status"]

    for row in mmo_rows
}


with tracks_out.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=[
            "genome",
            "Cluster00121_plot",
            "MMO_status",
        ]
    )


    writer.writeheader()


    for genome in tree_genomes:

        writer.writerow(
            {

                "genome":
                    genome,

                "Cluster00121_plot":
                    (
                        "Present"
                        if genome in
                        c121_positive
                        else ""
                    ),

                "MMO_status":
                    mmo_by_genome[
                        genome
                    ],
            }
        )


## ================================================================== ##
## Final QC
## ================================================================== ##

counts = defaultdict(
    int
)


for row in mmo_rows:

    counts[
        row["MMO_status"]
    ] += 1


partials = [

    row

    for row in mmo_rows

    if row[
        "partial_component_hits"
    ] ==
    1
]


print()
print("=" * 80)
print("CORRECTED MMO STATUS")
print("=" * 80)
print()


for status in (
    "pMMO",
    "sMMO",
    "Both",
    "None",
):

    print(
        f"{status:<5}: "
        f"{counts[status]}"
    )


print()
print(
    "Cluster_00121-positive tips:",
    len(
        c121_positive
    )
)


print()
print(
    "Partial component patterns:",
    len(
        partials
    )
)


for row in partials[:50]:

    print(
        f"  {row['genome']:<45} "
        f"{row['detected_components']}"
    )


## -------------------------------------------------------------- ##
## Refuse to silently accept another impossible all-None result.
## -------------------------------------------------------------- ##

if counts["None"] == 650:

    raise SystemExit(
        "\nERROR: all 650 genomes still classify as None.\n"
        "Do not plot this track."
    )


print()
print(
    "Corrected MMO table:"
)

print(
    " ",
    mmo_out
)


print()
print(
    "Corrected final track table:"
)

print(
    " ",
    tracks_out
)
PY


echo
echo "============================================================"
echo "STAGE 93E COMPLETE"
echo "============================================================"
echo

echo "MMO status counts:"
cut -f8 "$MMO_SUMMARY" \
    | tail -n +2 \
    | sort \
    | uniq -c

echo
echo "sMMO / Both genomes:"
awk -F '\t' \
    'NR == 1 || $8 == "sMMO" || $8 == "Both"' \
    "$MMO_SUMMARY"

echo
echo "Files:"
echo "  $BEST_HITS"
echo "  $MMO_SUMMARY"
echo "  $FINAL_TRACKS"
