#!/usr/bin/env python3

from pathlib import Path
import gzip
import pandas as pd


## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

MASTER = (
    PROJECT
    / "semibin2_analysis"
    / "MAG_master_table_SemiBin2.tsv"
)

DEPTH_DIR = (
    PROJECT
    / "semibin2_analysis"
    / "depth"
)

OUTDIR = (
    PROJECT
    / "semibin2_analysis"
)

PER_MAG_OUT = (
    OUTDIR
    / "SemiBin2_MAG_coverage_per_bin.tsv"
)

CATEGORY_OUT = (
    OUTDIR
    / "SemiBin2_MAG_coverage_by_category.tsv"
)


## ------------------------------------------------------------
## Read master table
## ------------------------------------------------------------

master = pd.read_csv(
    MASTER,
    sep="\t"
)


## ------------------------------------------------------------
## Open FASTA
##
## SemiBin2 may store bins as gzip-compressed FASTA files.
##
## Detect gzip from the file magic bytes rather than trusting
## the filename extension.
##
## gzip magic:
##
##   1f 8b
## ------------------------------------------------------------

def open_fasta(path):

    path = Path(path)

    with open(
        path,
        "rb"
    ) as handle:

        magic = handle.read(2)


    if magic == b"\x1f\x8b":

        return gzip.open(
            path,
            "rt"
        )

    return open(
        path,
        "rt"
    )


## ------------------------------------------------------------
## FASTA parser
##
## Returns:
##
##   [(contig_name, contig_length), ...]
## ------------------------------------------------------------

def fasta_contigs(path):

    records = []

    name = None
    length = 0


    with open_fasta(
        path
    ) as handle:

        for line in handle:

            if line.startswith(">"):

                if name is not None:

                    records.append(
                        (
                            name,
                            length
                        )
                    )

                name = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

                length = 0

            else:

                length += len(
                    line.strip()
                )


        if name is not None:

            records.append(
                (
                    name,
                    length
                )
            )


    return records


## ------------------------------------------------------------
## Cache depth tables
## ------------------------------------------------------------

depth_cache = {}


def get_depth(sample):

    if sample in depth_cache:

        return depth_cache[
            sample
        ]


    path = (
        DEPTH_DIR
        / f"{sample}_depth.txt"
    )


    if not path.exists():

        raise RuntimeError(
            f"Missing depth table: {path}"
        )


    depth = pd.read_csv(
        path,
        sep="\t"
    )


    required = {
        "contigName",
        "contigLen",
        "totalAvgDepth"
    }


    missing = (
        required
        - set(
            depth.columns
        )
    )


    if missing:

        raise RuntimeError(
            f"Missing columns in {path}: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )


    depth = (
        depth[
            [
                "contigName",
                "contigLen",
                "totalAvgDepth"
            ]
        ]
        .set_index(
            "contigName"
        )
    )


    depth_cache[
        sample
    ] = depth


    return depth


## ------------------------------------------------------------
## Calculate coverage
##
## MAG mean coverage:
##
##   sum(contig length × contig mean depth)
##   ---------------------------------------
##              sum(contig length)
##
## This is length-weighted mean coverage.
##
## We do NOT sum contig coverage directly because highly
## fragmented MAGs would otherwise appear artificially more
## abundant.
## ------------------------------------------------------------

rows = []


for _, mag in master.iterrows():

    genome_id = str(
        mag[
            "Genome_ID"
        ]
    )

    sample = str(
        mag[
            "Sample"
        ]
    )

    fasta = Path(
        str(
            mag[
                "Original_path"
            ]
        )
    )


    if not fasta.exists():

        raise RuntimeError(
            f"{genome_id}: FASTA does not exist:\n"
            f"{fasta}"
        )


    depth = get_depth(
        sample
    )

    contigs = fasta_contigs(
        fasta
    )


    if len(contigs) == 0:

        raise RuntimeError(
            f"{genome_id}: no contigs found in FASTA:\n"
            f"{fasta}"
        )


    total_length = 0.0

    weighted_depth = 0.0

    missing_contigs = []


    for contig, fasta_length in contigs:

        if contig not in depth.index:

            missing_contigs.append(
                contig
            )

            continue


        contig_length = float(
            depth.loc[
                contig,
                "contigLen"
            ]
        )

        contig_depth = float(
            depth.loc[
                contig,
                "totalAvgDepth"
            ]
        )


        total_length += (
            contig_length
        )

        weighted_depth += (
            contig_length
            * contig_depth
        )


    ## --------------------------------------------------------
    ## Fail rather than silently calculate partial coverage
    ## --------------------------------------------------------

    if missing_contigs:

        examples = ", ".join(
            missing_contigs[:5]
        )

        raise RuntimeError(
            f"{genome_id}: "
            f"{len(missing_contigs)} contigs "
            f"were missing from the depth table.\n"
            f"Examples: {examples}"
        )


    if total_length <= 0:

        raise RuntimeError(
            f"{genome_id}: total genome length is zero."
        )


    mean_coverage = (
        weighted_depth
        / total_length
    )


    ## --------------------------------------------------------
    ## Taxonomic plotting category
    ##
    ## Keep Methylococcales at genus level.
    ## Collapse everything else into Other.
    ## --------------------------------------------------------

    order = mag.get(
        "Order"
    )


    if (
        not pd.isna(order)
        and str(order)
        == "Methylococcales"
    ):

        genus = mag.get(
            "Genus"
        )


        if pd.isna(
            genus
        ):

            category = (
                "Methylococcales_unclassified"
            )

        else:

            category = str(
                genus
            )

    else:

        category = "Other"


    rows.append({

        "Genome_ID":
            genome_id,

        "Sample":
            sample,

        "Category":
            category,

        "Genus":
            mag.get(
                "Genus"
            ),

        "Species":
            mag.get(
                "Species"
            ),

        "Completeness":
            mag[
                "Completeness"
            ],

        "Contamination":
            mag[
                "Contamination"
            ],

        "Quality":
            mag[
                "Quality"
            ],

        "Genome_length_bp":
            int(
                total_length
            ),

        "N_contigs":
            len(
                contigs
            ),

        "Mean_coverage_x":
            mean_coverage
    })


## ------------------------------------------------------------
## Per-MAG table
## ------------------------------------------------------------

per_mag = pd.DataFrame(
    rows
)


per_mag = per_mag.sort_values(
    [
        "Sample",
        "Category",
        "Mean_coverage_x"
    ],
    ascending=[
        True,
        True,
        False
    ]
)


per_mag.to_csv(
    PER_MAG_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Sum MAG coverage by category
##
## IMPORTANT:
##
## This is preliminary.
##
## Multiple same-genus MAGs are currently treated as separate
## populations and their mean coverages are summed.
##
## We will subsequently check same-genus MAGs using ANI,
## completeness, genome size and coverage before deciding
## whether any represent fragments of the same population.
## ------------------------------------------------------------

category = (
    per_mag
    .groupby(
        [
            "Sample",
            "Category"
        ],
        as_index=False
    )
    ["Mean_coverage_x"]
    .sum()
    .rename(
        columns={
            "Mean_coverage_x":
                "Summed_MAG_mean_coverage_x"
        }
    )
)


category.to_csv(
    CATEGORY_OUT,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Print Methylococcales only
## ------------------------------------------------------------

meth = per_mag[
    per_mag[
        "Category"
    ] != "Other"
].copy()


print()

print(
    "SemiBin2 Methylococcales MAG coverage"
)

print(
    "=" * 120
)


if len(meth) == 0:

    print(
        "No Methylococcales MAGs found."
    )

else:

    print(

        meth[
            [
                "Sample",
                "Genome_ID",
                "Category",
                "Species",
                "Completeness",
                "Contamination",
                "Genome_length_bp",
                "N_contigs",
                "Mean_coverage_x"
            ]
        ].to_string(
            index=False
        )
    )


print()

print(
    f"Per-MAG table:\n"
    f"{PER_MAG_OUT}"
)

print()

print(
    f"Category table:\n"
    f"{CATEGORY_OUT}"
)
