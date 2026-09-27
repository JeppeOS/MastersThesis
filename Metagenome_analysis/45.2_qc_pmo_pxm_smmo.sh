#!/bin/bash

#SBATCH --job-name=pmo_pxm_smmo_qc
#SBATCH --account=methanotrophs
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=comparative_analysis/final_species_phylogeny/logs/93F_pmo_pxm_smmo_%j.out
#SBATCH --error=comparative_analysis/final_species_phylogeny/logs/93F_pmo_pxm_smmo_%j.err


set -euo pipefail


## ================================================================== ##
## STAGE 93F
##
## QC / RECLASSIFICATION OF METHANE MONOOXYGENASE SYSTEMS
##
## Across all 650 final-tree genomes:
##
##   canonical pMMO:
##       pmoA / pmoB / pmoC
##
##   divergent pXMO:
##       pxmA / pxmB / pxmC
##
##   soluble MMO:
##       mmoX / mmoY / mmoZ
##
##
## WHY:
##
## Pmo and Pxm are homologous CuMMO proteins.
##
## A search containing only canonical Pmo references can therefore
## incorrectly classify Pxm proteins as Pmo.
##
## We instead search every protein simultaneously against:
##
##   pmoA/B/C
##   pxmA/B/C
##   mmoX/Y/Z
##
## and assign each query protein to its BEST reference-family hit.
##
##
## FINAL SYSTEM RULES:
##
## canonical pMMO =
##     pmoA AND (pmoB OR pmoC)
##
## pXMO =
##     pxmA AND (pxmB OR pxmC)
##
## sMMO =
##     mmoX AND (mmoY OR mmoZ)
##
##
## This is a QC stage.
##
## It DOES NOT overwrite the final plotting track yet.
## ================================================================== ##


ROOT="$HOME/methanotrophs/methanotroph_project/jeppe"

PROJECT="$ROOT/metagenome"

PHYLO="$PROJECT/comparative_analysis/final_species_phylogeny"

EXTRA="$PHYLO/08_visualization/extra_tracks"

REFDIR="$EXTRA/references"

MMSEQSDIR="$EXTRA/mmseqs"

TMPDIR="$EXTRA/tmp/93F_pmo_pxm_smmo"

LOGDIR="$PHYLO/logs"


mkdir -p \
    "$REFDIR" \
    "$MMSEQSDIR" \
    "$TMPDIR" \
    "$LOGDIR"


## ================================================================== ##
## Existing Stage 93 inputs
## ================================================================== ##

ALL_PROTEINS="$EXTRA/proteins/all_650_tree_genomes_proteins.faa"

CANONICAL_REFS="$REFDIR/MMO_reference_proteins.faa"

OLD_STATUS="$EXTRA/all_650_MMO_status.tsv"

META="$PHYLO/08_visualization/final_650_tip_metadata.tsv"


## ================================================================== ##
## New outputs
## ================================================================== ##

BG8_GENBANK="$REFDIR/Methylomicrobium_album_BG8_NZ_CM001475.gb"

PXM_REFS="$REFDIR/BG8_pxmABC_reference_proteins.faa"

ALL_REFS="$REFDIR/pmo_pxm_smmo_reference_proteins.faa"


RAW_HITS="$MMSEQSDIR/pmo_pxm_smmo_all650_raw_hits.tsv"

BEST_HITS="$MMSEQSDIR/pmo_pxm_smmo_all650_best_hits.tsv"


SYSTEM_TABLE="$EXTRA/all_650_pmo_pxm_smmo_QC.tsv"

OLD_NONE_QC="$EXTRA/original_209_None_QC.tsv"


## ================================================================== ##
## Input QC
## ================================================================== ##

for F in \
    "$ALL_PROTEINS" \
    "$CANONICAL_REFS" \
    "$OLD_STATUS" \
    "$META"

do

    if [[ ! -s "$F" ]]
    then

        echo "ERROR: required file missing:"
        echo "$F"

        exit 1

    fi

done


echo
echo "============================================================"
echo "STAGE 93F"
echo "============================================================"
echo

echo "Existing canonical/sMMO references:"

grep '^>' "$CANONICAL_REFS"

echo


## ================================================================== ##
## PART A
##
## Download the annotated Methylomicrobium album BG8 chromosome.
##
## BG8 contains:
##
##   pmoCAB:
##       METAL_RS17430
##       METAL_RS17425
##       METAL_RS17420
##
##   pxmABC:
##       METAL_RS06980
##       METAL_RS06975
##       METAL_RS06970
##
## We only extract pxmABC here because canonical Pmo references
## already exist from Stage 93.
## ================================================================== ##

curl \
    -fL \
    --retry 5 \
    --retry-delay 3 \
    --connect-timeout 30 \
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CM001475&rettype=gbwithparts&retmode=text' \
    -o "$BG8_GENBANK"


if [[ ! -s "$BG8_GENBANK" ]]
then

    echo "ERROR: BG8 GenBank download failed."

    exit 1

fi


echo "Downloaded BG8 GenBank record:"
echo "  $BG8_GENBANK"
echo


## ================================================================== ##
## PART B
##
## Extract pxmA / pxmB / pxmC translations from the GenBank record.
##
## Literature mapping:
##
##   pxmA = METAL_RS06980
##   pxmB = METAL_RS06975
##   pxmC = METAL_RS06970
##
## Standard-library Python only.
## ================================================================== ##

python - \
    "$BG8_GENBANK" \
    "$PXM_REFS" \
<<'PY'
import re
import sys
from pathlib import Path


gb_file = Path(sys.argv[1])
out_file = Path(sys.argv[2])


wanted = {

    "METAL_RS06980":
        "pxmA",

    "METAL_RS06975":
        "pxmB",

    "METAL_RS06970":
        "pxmC",
}


text = gb_file.read_text(
    errors="ignore"
)


## ------------------------------------------------------------------ ##
## Split GenBank FEATURES into CDS-sized chunks.
## ------------------------------------------------------------------ ##

chunks = re.split(
    r"\n\s{5}(?=CDS\s)",
    text
)


found = {}


for chunk in chunks:

    if not chunk.startswith(
        "CDS"
    ):

        continue


    locus_match = re.search(
        r'/locus_tag="([^"]+)"',
        chunk
    )


    if not locus_match:

        continue


    locus = locus_match.group(
        1
    )


    if locus not in wanted:

        continue


    translation_match = re.search(
        r'/translation="([^"]+)"',
        chunk,
        flags=re.S
    )


    if not translation_match:

        raise SystemExit(
            f"ERROR: no translation found for {locus}"
        )


    sequence = re.sub(
        r"\s+",
        "",
        translation_match.group(
            1
        )
    )


    found[
        locus
    ] = sequence


missing = [
    locus
    for locus in wanted
    if locus not in found
]


if missing:

    raise SystemExit(
        "ERROR: failed to recover BG8 pXMO proteins: "
        +
        ", ".join(
            missing
        )
    )


with out_file.open(
    "w"
) as handle:

    for locus, component in wanted.items():

        sequence = found[
            locus
        ]


        handle.write(
            f">{component}|BG8_{locus}\n"
        )


        for i in range(
            0,
            len(sequence),
            80
        ):

            handle.write(
                sequence[i:i+80]
                +
                "\n"
            )


print()
print("Recovered BG8 pXMO references:")

for locus, component in wanted.items():

    print(
        f"  {component:<4} "
        f"{locus:<15} "
        f"{len(found[locus])} aa"
    )
PY


echo
echo "pXMO reference FASTA:"
grep '^>' "$PXM_REFS"


N_PXM=$(grep -c '^>' "$PXM_REFS")


if [[ "$N_PXM" -ne 3 ]]
then

    echo "ERROR: expected exactly 3 pXMO proteins."

    exit 1

fi


## ================================================================== ##
## PART C
##
## Build one nine-protein reference FASTA:
##
##   pmoA
##   pmoB
##   pmoC
##
##   pxmA
##   pxmB
##   pxmC
##
##   mmoX
##   mmoY
##   mmoZ
## ================================================================== ##

cat \
    "$CANONICAL_REFS" \
    "$PXM_REFS" \
    > "$ALL_REFS"


echo
echo "Combined reference set:"

grep '^>' "$ALL_REFS"


N_REF=$(grep -c '^>' "$ALL_REFS")


if [[ "$N_REF" -ne 9 ]]
then

    echo
    echo "ERROR: expected exactly 9 reference proteins;"
    echo "found $N_REF."

    exit 1

fi


## ================================================================== ##
## PART D
##
## Search ALL 650 genome proteins against all nine references.
##
## Initial MMseqs filter is deliberately slightly permissive.
##
## Final thresholding is performed below:
##
##   identity >= 35%
##   qcov     >= 70%
##   tcov     >= 70%
## ================================================================== ##

source "$HOME/miniforge3/etc/profile.d/conda.sh"

conda activate MMseqs2


rm -rf "$TMPDIR/mmseqs"


mmseqs easy-search \
    "$ALL_PROTEINS" \
    "$ALL_REFS" \
    "$RAW_HITS" \
    "$TMPDIR/mmseqs" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --min-seq-id 0.30 \
    -c 0.60 \
    --cov-mode 0 \
    --max-seqs 20 \
    --format-output \
    "query,target,fident,alnlen,qcov,tcov,evalue,bits"


echo
echo "Raw search hits:"

wc -l "$RAW_HITS"


## ================================================================== ##
## PART E
##
## Best-family assignment + genome-level classification.
##
## IMPORTANT:
##
## A query protein can hit both Pmo and Pxm references.
##
## We therefore:
##
##   1. retain all hits passing 35% / 70% / 70%
##   2. compare them by query protein
##   3. retain the highest-bitscore reference
##
## This explicitly separates canonical Pmo-like and Pxm-like proteins.
## ================================================================== ##

python - \
    "$META" \
    "$RAW_HITS" \
    "$OLD_STATUS" \
    "$BEST_HITS" \
    "$SYSTEM_TABLE" \
    "$OLD_NONE_QC" \
<<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path


meta_file = Path(
    sys.argv[1]
)

raw_hits_file = Path(
    sys.argv[2]
)

old_status_file = Path(
    sys.argv[3]
)

best_hits_file = Path(
    sys.argv[4]
)

system_table_file = Path(
    sys.argv[5]
)

old_none_file = Path(
    sys.argv[6]
)


## ================================================================== ##
## Thresholds
##
## MMseqs fident / qcov / tcov are fractions.
## ================================================================== ##

MIN_IDENT = 0.35

MIN_QCOV = 0.70

MIN_TCOV = 0.70


VALID = {

    "pmoA",
    "pmoB",
    "pmoC",

    "pxmA",
    "pxmB",
    "pxmC",

    "mmoX",
    "mmoY",
    "mmoZ",
}


## ================================================================== ##
## Tree genomes
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
        f"ERROR: expected 650 genomes; "
        f"found {len(tree_genomes)}."
    )


tree_genome_set = set(
    tree_genomes
)


## ================================================================== ##
## Existing status
##
## We retain it ONLY to identify the original 209 "None" genomes.
## ================================================================== ##

old_status = {}


with old_status_file.open() as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )


    for row in reader:

        old_status[
            row["genome"]
        ] = row["MMO_status"]


original_none = {

    genome

    for genome, status
    in old_status.items()

    if status ==
        "None"
}


print()
print(
    "Original None genomes:",
    len(
        original_none
    )
)


## ================================================================== ##
## Parse hits
##
## Store all PASSING hits for each query protein.
## ================================================================== ##

passing_by_query = defaultdict(
    list
)


raw_n = 0

pass_n = 0


with raw_hits_file.open() as handle:

    reader = csv.reader(
        handle,
        delimiter="\t"
    )


    for row in reader:

        if len(row) != 8:

            continue


        raw_n += 1


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

            fident = float(
                fident
            )

            qcov = float(
                qcov
            )

            tcov = float(
                tcov
            )

            bits = float(
                bits
            )

        except ValueError:

            continue


        if fident < MIN_IDENT:

            continue


        if qcov < MIN_QCOV:

            continue


        if tcov < MIN_TCOV:

            continue


        if "|" not in query:

            continue


        if "|" not in target:

            continue


        genome, protein = query.split(
            "|",
            1
        )


        component, reference = target.split(
            "|",
            1
        )


        if genome not in tree_genome_set:

            continue


        if component not in VALID:

            continue


        pass_n += 1


        passing_by_query[
            query
        ].append(
            {

                "genome":
                    genome,

                "protein":
                    protein,

                "component":
                    component,

                "reference":
                    reference,

                "fident":
                    fident,

                "qcov":
                    qcov,

                "tcov":
                    tcov,

                "evalue":
                    evalue,

                "bits":
                    bits,
            }
        )


print(
    "Raw alignments:",
    raw_n
)

print(
    "Passing alignments:",
    pass_n
)


## ================================================================== ##
## Best assignment per query protein.
##
## Also retain the second-best hit so we can quantify whether a
## Pmo/Pxm classification is marginal.
## ================================================================== ##

best_hits = []


for query, hits in passing_by_query.items():

    hits = sorted(
        hits,
        key=lambda x: x["bits"],
        reverse=True
    )


    best = hits[0]


    second = (
        hits[1]
        if len(hits) > 1
        else None
    )


    second_component = (
        second["component"]
        if second
        else ""
    )


    second_bits = (
        second["bits"]
        if second
        else 0.0
    )


    bit_margin = (
        best["bits"]
        -
        second_bits
    )


    relative_margin = (
        bit_margin
        /
        best["bits"]
        if best["bits"] > 0
        else 0.0
    )


    best_hits.append(
        {

            **best,

            "second_component":
                second_component,

            "second_bits":
                second_bits,

            "bit_margin":
                bit_margin,

            "relative_bit_margin":
                relative_margin,
        }
    )


## ================================================================== ##
## Write detailed protein classification
## ================================================================== ##

best_fields = [

    "genome",
    "protein",
    "component",
    "reference",

    "fident",
    "qcov",
    "tcov",
    "evalue",
    "bits",

    "second_component",
    "second_bits",
    "bit_margin",
    "relative_bit_margin",
]


with best_hits_file.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=best_fields
    )


    writer.writeheader()

    writer.writerows(
        best_hits
    )


## ================================================================== ##
## Components by genome
## ================================================================== ##

components = defaultdict(
    set
)


for hit in best_hits:

    components[
        hit["genome"]
    ].add(
        hit["component"]
    )


## ================================================================== ##
## System calls
## ================================================================== ##

def canonical_pmmo(comp):

    return (
        "pmoA" in comp
        and
        (
            "pmoB" in comp
            or
            "pmoC" in comp
        )
    )


def pxmo(comp):

    return (
        "pxmA" in comp
        and
        (
            "pxmB" in comp
            or
            "pxmC" in comp
        )
    )


def smmo(comp):

    return (
        "mmoX" in comp
        and
        (
            "mmoY" in comp
            or
            "mmoZ" in comp
        )
    )


rows = []


for genome in tree_genomes:

    comp = components[
        genome
    ]


    has_pmo = canonical_pmmo(
        comp
    )

    has_pxm = pxmo(
        comp
    )

    has_smmo = smmo(
        comp
    )


    canonical_components = sorted(
        comp.intersection(
            {
                "pmoA",
                "pmoB",
                "pmoC",
            }
        )
    )


    pxm_components = sorted(
        comp.intersection(
            {
                "pxmA",
                "pxmB",
                "pxmC",
            }
        )
    )


    smmo_components = sorted(
        comp.intersection(
            {
                "mmoX",
                "mmoY",
                "mmoZ",
            }
        )
    )


    partial_pmo = (
        bool(
            canonical_components
        )
        and
        not has_pmo
    )


    partial_pxm = (
        bool(
            pxm_components
        )
        and
        not has_pxm
    )


    partial_smmo = (
        bool(
            smmo_components
        )
        and
        not has_smmo
    )


    if has_pmo and has_smmo:

        canonical_mmo_status = (
            "Both"
        )

    elif has_pmo:

        canonical_mmo_status = (
            "pMMO"
        )

    elif has_smmo:

        canonical_mmo_status = (
            "sMMO"
        )

    else:

        canonical_mmo_status = (
            "None"
        )


    if has_pmo:

        particulate_status = (
            "canonical_pMMO"
        )

    elif has_pxm:

        particulate_status = (
            "pXMO_only"
        )

    elif partial_pmo:

        particulate_status = (
            "partial_pMMO"
        )

    elif partial_pxm:

        particulate_status = (
            "partial_pXMO"
        )

    else:

        particulate_status = (
            "none_detected"
        )


    rows.append(
        {

            "genome":
                genome,

            "original_MMO_status":
                old_status.get(
                    genome,
                    ""
                ),

            "canonical_pMMO":
                int(
                    has_pmo
                ),

            "pXMO":
                int(
                    has_pxm
                ),

            "sMMO":
                int(
                    has_smmo
                ),

            "partial_pMMO":
                int(
                    partial_pmo
                ),

            "partial_pXMO":
                int(
                    partial_pxm
                ),

            "partial_sMMO":
                int(
                    partial_smmo
                ),

            "canonical_MMO_status":
                canonical_mmo_status,

            "particulate_status":
                particulate_status,

            "pmo_components":
                ",".join(
                    canonical_components
                ),

            "pxm_components":
                ",".join(
                    pxm_components
                ),

            "smmo_components":
                ",".join(
                    smmo_components
                ),

            "all_components":
                ",".join(
                    sorted(
                        comp
                    )
                ),
        }
    )


fields = [

    "genome",

    "original_MMO_status",

    "canonical_pMMO",
    "pXMO",
    "sMMO",

    "partial_pMMO",
    "partial_pXMO",
    "partial_sMMO",

    "canonical_MMO_status",

    "particulate_status",

    "pmo_components",
    "pxm_components",
    "smmo_components",

    "all_components",
]


with system_table_file.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fields
    )


    writer.writeheader()

    writer.writerows(
        rows
    )


## ================================================================== ##
## Focus specifically on original 209 "None" genomes
## ================================================================== ##

none_rows = [

    row

    for row in rows

    if row[
        "genome"
    ] in original_none
]


with old_none_file.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fields
    )


    writer.writeheader()

    writer.writerows(
        none_rows
    )


## ================================================================== ##
## Summaries
## ================================================================== ##

def count(
    rows,
    column,
    value=1
):

    return sum(
        row[column] ==
        value

        for row in rows
    )


print()
print("=" * 80)
print("ALL 650 GENOMES")
print("=" * 80)

print()
print(
    "Canonical pMMO:",
    count(
        rows,
        "canonical_pMMO"
    )
)

print(
    "pXMO:",
    count(
        rows,
        "pXMO"
    )
)

print(
    "sMMO:",
    count(
        rows,
        "sMMO"
    )
)


print()
print(
    "pMMO + pXMO:",
    sum(
        row["canonical_pMMO"] == 1
        and
        row["pXMO"] == 1
        for row in rows
    )
)

print(
    "pXMO only:",
    sum(
        row["canonical_pMMO"] == 0
        and
        row["pXMO"] == 1
        for row in rows
    )
)


print()
print("=" * 80)
print("ORIGINAL 209 NONE GENOMES")
print("=" * 80)

print()
print(
    "Original None count:",
    len(
        none_rows
    )
)


print(
    "Now canonical pMMO:",
    count(
        none_rows,
        "canonical_pMMO"
    )
)

print(
    "Now pXMO:",
    count(
        none_rows,
        "pXMO"
    )
)

print(
    "Now sMMO:",
    count(
        none_rows,
        "sMMO"
    )
)


print()
print(
    "pXMO-only:",
    sum(
        row["canonical_pMMO"] == 0
        and
        row["pXMO"] == 1
        for row in none_rows
    )
)


print(
    "Partial canonical pMMO:",
    count(
        none_rows,
        "partial_pMMO"
    )
)


print(
    "Partial pXMO:",
    count(
        none_rows,
        "partial_pXMO"
    )
)


print(
    "Still no pMMO, pXMO or sMMO components:",
    sum(
        not row[
            "all_components"
        ]
        for row in none_rows
    )
)


print()
print("Original-None particulate-status breakdown:")


status_counts = defaultdict(
    int
)


for row in none_rows:

    status_counts[
        row[
            "particulate_status"
        ]
    ] += 1


for status in sorted(
    status_counts
):

    print(
        f"  {status:<20} "
        f"{status_counts[status]}"
    )


print()
print(
    "Detailed QC table:"
)

print(
    " ",
    system_table_file
)


print()
print(
    "Original 209 None QC:"
)

print(
    " ",
    old_none_file
)
PY


echo
echo "============================================================"
echo "STAGE 93F COMPLETE"
echo "============================================================"
echo

echo "Main QC table:"
echo "  $SYSTEM_TABLE"

echo

echo "Original-None QC:"
echo "  $OLD_NONE_QC"

echo

echo "First 20 original None genomes after QC:"
echo

column -t -s $'\t' "$OLD_NONE_QC" \
    | head -n 21
