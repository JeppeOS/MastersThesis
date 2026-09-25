#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd


## ================================================================== ##
## STAGE 50
##
## BUILD METHYLOCOCCALES ACTUAL-GENE CATALOGUE
##
## Purpose
## -------
##
## Build a complete coordinate catalogue for every Prodigal-predicted
## protein in the 16 SemiBin2 Methylococcales MAGs.
##
## Existing functional evidence is attached when available:
##
##   - FeGenie
##   - FindMeHemes
##   - SignalP
##   - DeepTMHMM
##   - integrated localization
##   - fixed GlobDB MMseq cluster mapping
##   - old GlobDB MCL module membership of the TOP cluster
##
## Ordinary genes are retained even when they have no functional
## annotation.
##
## The 21 current priority focal proteins are also extracted:
##
##   - 6 Cyc2_repCluster2 proteins
##   - 15 proteins with >=5 predicted hemes
##
## No COG annotation is added here.
## No neighborhood extraction is performed here.
## ================================================================== ##


## ================================================================== ##
## Paths
## ================================================================== ##

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)


OLD_WORKFLOW = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "genome_analysis_workflow"
)


ORF_DIR = (
    PROJECT
    / "functional_analysis"
    / "FeGenie_SemiBin2_full_run_20260918"
    / "fegenie_results"
    / "ORF_calls"
)


METHYLOCOCCALES_MAGS = (
    PROJECT
    / "semibin2_analysis"
    / "Methylococcales_MAGs_SemiBin2.tsv"
)


INTEGRATED_MASTER = (
    PROJECT
    / "functional_analysis"
    / "integrated_annotations"
    / "Methylococcales"
    / "Methylococcales_MASTER_PROTEIN_ANNOTATIONS.tsv"
)


GLOBDB_ASSIGNMENTS = (
    PROJECT
    / "functional_analysis"
    / "globdb_cluster_mapping"
    / "Methylococcales_GlobDB_cluster_assignments.tsv"
)


OLD_MODULE_MEMBERSHIP = (
    OLD_WORKFLOW
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)


OUT_DIR = (
    PROJECT
    / "functional_analysis"
    / "methylococcales_neighborhoods"
)


OUT_COORDINATES = (
    OUT_DIR
    / "all_methylococcales_gene_coordinates.tsv"
)


OUT_CATALOG = (
    OUT_DIR
    / "methylococcales_gene_catalog.tsv"
)


OUT_FOCALS = (
    OUT_DIR
    / "priority_focal_genes.tsv"
)


OUT_MAG_METADATA = (
    OUT_DIR
    / "methylococcales_mag_metadata.tsv"
)


OUT_GENOME_QC = (
    OUT_DIR
    / "genome_coordinate_qc.tsv"
)


OUT_QC = (
    OUT_DIR
    / "gene_catalog_qc.tsv"
)


OUT_UNMATCHED_CANDIDATES = (
    OUT_DIR
    / "unmatched_integrated_candidates.tsv"
)


OUT_UNMATCHED_FOCALS = (
    OUT_DIR
    / "unmatched_priority_focal_genes.tsv"
)


## ================================================================== ##
## Locked expectations from the current analysis
## ================================================================== ##

EXPECTED_MAGS = 16

EXPECTED_CYC2 = 6

EXPECTED_GE5_HEME = 15

EXPECTED_PRIORITY_FOCALS = 21


## ================================================================== ##
## Prodigal FASTA header parser
##
## Example:
##
## >contig_25_317 # 1234 # 2763 # 1 # ID=25_317;partial=00;...
##
## ================================================================== ##

HEADER_RE = re.compile(
    r"^>(\S+)"
    r"\s*#\s*(\d+)"
    r"\s*#\s*(\d+)"
    r"\s*#\s*(-?1)"
    r"\s*#\s*(.*)$"
)


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr,
    )

    sys.exit(1)


def clean(value):

    if pd.isna(value):
        return ""

    return str(value).strip()


def require_file(path):

    if not path.is_file():

        fail(
            f"Required file does not exist:\n"
            f"{path}"
        )


def parse_attributes(text):

    attributes = {}

    for field in text.split(";"):

        field = field.strip()

        if not field:
            continue

        if "=" not in field:
            continue

        key, value = field.split(
            "=",
            1,
        )

        attributes[
            key.strip()
        ] = value.strip()

    return attributes


def normalise_genome_filename(path):

    name = path.name

    suffixes = (
        "-proteins.faa",
        ".fasta",
        ".fna",
        ".faa",
        ".fa",
    )

    changed = True

    while changed:

        changed = False

        for suffix in suffixes:

            if name.endswith(suffix):

                name = name[
                    :-len(suffix)
                ]

                changed = True

                break

    return name


def parse_header(
    header,
    genome,
    source_fasta,
):

    match = HEADER_RE.match(
        header.rstrip()
    )

    if match is None:

        fail(
            "Could not parse Prodigal FASTA header:\n"
            f"{header}"
        )


    protein_id = match.group(1)

    start = int(
        match.group(2)
    )

    end = int(
        match.group(3)
    )

    strand_numeric = int(
        match.group(4)
    )

    attribute_string = (
        match.group(5)
    )


    if start < 1:

        fail(
            f"Invalid start coordinate for "
            f"{genome} / {protein_id}: "
            f"{start}"
        )


    if end < start:

        fail(
            f"End < start for "
            f"{genome} / {protein_id}: "
            f"{start}-{end}"
        )


    attributes = parse_attributes(
        attribute_string
    )


    prodigal_id = attributes.get(
        "ID",
        ""
    )


    if not prodigal_id:

        fail(
            f"No Prodigal ID= field for "
            f"{genome} / {protein_id}"
        )


    if "_" not in prodigal_id:

        fail(
            f"Unexpected Prodigal ID for "
            f"{genome} / {protein_id}: "
            f"{prodigal_id}"
        )


    seqnum, gene_ordinal = (
        prodigal_id.rsplit(
            "_",
            1,
        )
    )


    if (
        not seqnum.isdigit()
        or
        not gene_ordinal.isdigit()
    ):

        fail(
            f"Unexpected Prodigal ID for "
            f"{genome} / {protein_id}: "
            f"{prodigal_id}"
        )


    ## -------------------------------------------------------------- ##
    ## Recover original contig ID.
    ##
    ## Prodigal protein ID:
    ##
    ##     contig_25_317
    ##
    ## Prodigal ID:
    ##
    ##     25_317
    ##
    ## Gene ordinal:
    ##
    ##     317
    ##
    ## Therefore contig:
    ##
    ##     contig_25
    ## -------------------------------------------------------------- ##

    expected_suffix = (
        "_"
        + gene_ordinal
    )


    if not protein_id.endswith(
        expected_suffix
    ):

        fail(
            "Protein ID does not end in its "
            "Prodigal gene ordinal:\n"
            f"genome      = {genome}\n"
            f"protein_id  = {protein_id}\n"
            f"Prodigal ID = {prodigal_id}"
        )


    contig = protein_id[
        :-len(expected_suffix)
    ]


    if strand_numeric == 1:

        strand = "+"

    elif strand_numeric == -1:

        strand = "-"

    else:

        fail(
            f"Unexpected strand for "
            f"{genome} / {protein_id}: "
            f"{strand_numeric}"
        )


    return {

        "genome":
            genome,

        "protein_id":
            protein_id,

        "contig":
            contig,

        "start":
            start,

        "end":
            end,

        "gene_length_nt":
            end - start + 1,

        "strand":
            strand,

        "prodigal_strand":
            strand_numeric,

        "prodigal_id":
            prodigal_id,

        "prodigal_seqnum":
            int(
                seqnum
            ),

        "prodigal_gene_ordinal":
            int(
                gene_ordinal
            ),

        "partial":
            attributes.get(
                "partial",
                ""
            ),

        "start_type":
            attributes.get(
                "start_type",
                ""
            ),

        "rbs_motif":
            attributes.get(
                "rbs_motif",
                ""
            ),

        "rbs_spacer":
            attributes.get(
                "rbs_spacer",
                ""
            ),

        "gc_cont":
            attributes.get(
                "gc_cont",
                ""
            ),

        "source_protein_fasta":
            str(
                source_fasta
            ),
    }


def numeric_flag(series):

    return (
        pd.to_numeric(
            series,
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )


## ================================================================== ##
## Start
## ================================================================== ##

print("=" * 80)
print("STAGE 50 - BUILD METHYLOCOCCALES GENE CATALOGUE")
print("=" * 80)


for path in (
    METHYLOCOCCALES_MAGS,
    INTEGRATED_MASTER,
    GLOBDB_ASSIGNMENTS,
    OLD_MODULE_MEMBERSHIP,
):

    require_file(
        path
    )


if not ORF_DIR.is_dir():

    fail(
        f"FeGenie ORF directory does not exist:\n"
        f"{ORF_DIR}"
    )


OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


## ================================================================== ##
## 1. Read authoritative 16-MAG table
## ================================================================== ##

print()
print("Reading authoritative Methylococcales MAG table...")


mag_meta = pd.read_csv(
    METHYLOCOCCALES_MAGS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


if "Genome_ID" not in mag_meta.columns:

    fail(
        "Methylococcales MAG table does not contain "
        "'Genome_ID'.\n"
        f"Available columns:\n"
        + "\n".join(
            mag_meta.columns
        )
    )


if mag_meta[
    "Genome_ID"
].duplicated().any():

    fail(
        "Duplicate Genome_ID values in "
        "Methylococcales MAG table."
    )


if len(
    mag_meta
) != EXPECTED_MAGS:

    fail(
        f"Expected {EXPECTED_MAGS} Methylococcales MAGs "
        f"but found {len(mag_meta)}."
    )


genomes = sorted(
    mag_meta[
        "Genome_ID"
    ]
    .astype(str)
    .str.strip()
)


genome_set = set(
    genomes
)


print(
    f"  Methylococcales MAGs: {len(genomes)}"
)


mag_meta.to_csv(
    OUT_MAG_METADATA,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 2. Read integrated candidate annotations
## ================================================================== ##

print()
print("Reading integrated candidate annotations...")


integrated = pd.read_csv(
    INTEGRATED_MASTER,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_integrated = {
    "genome",
    "protein_id",
    "candidate_source",
    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "export_evidence",
    "localization_class",
}


missing = (
    required_integrated
    - set(
        integrated.columns
    )
)


if missing:

    fail(
        "Integrated Methylococcales table is missing:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


if integrated[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys in "
        "integrated annotation table."
    )


unexpected_genomes = sorted(
    set(
        integrated[
            "genome"
        ]
    )
    - genome_set
)


if unexpected_genomes:

    fail(
        "Integrated table contains genomes outside the "
        "authoritative 16-MAG set:\n"
        + "\n".join(
            unexpected_genomes
        )
    )


print(
    f"  Integrated candidate proteins: "
    f"{len(integrated):,}"
)


## ================================================================== ##
## 3. Define current priority focal proteins
## ================================================================== ##

print()
print("Defining priority focal proteins...")


integrated[
    "number_of_hemes_numeric"
] = pd.to_numeric(
    integrated[
        "number_of_hemes"
    ],
    errors="coerce",
).fillna(0).astype(int)


integrated[
    "priority_cyc2"
] = (
    integrated[
        "fegenie_HMMs"
    ]
    .fillna("")
    .str.split(";")
    .apply(
        lambda values:
            int(
                "Cyc2_repCluster2"
                in {
                    clean(x)
                    for x
                    in values
                    if clean(x)
                }
            )
    )
)


integrated[
    "priority_ge5_hemes"
] = (
    integrated[
        "number_of_hemes_numeric"
    ]
    >= 5
).astype(int)


def priority_reason(row):

    reasons = []

    if int(
        row[
            "priority_cyc2"
        ]
    ) == 1:

        reasons.append(
            "Cyc2"
        )


    if int(
        row[
            "priority_ge5_hemes"
        ]
    ) == 1:

        reasons.append(
            ">=5_hemes"
        )


    return ";".join(
        reasons
    )


integrated[
    "priority_reason"
] = integrated.apply(
    priority_reason,
    axis=1,
)


priority = integrated[
    integrated[
        "priority_reason"
    ]
    != ""
].copy()


n_cyc2 = int(
    priority[
        "priority_cyc2"
    ].sum()
)


n_ge5 = int(
    priority[
        "priority_ge5_hemes"
    ].sum()
)


n_priority = len(
    priority
)


print(
    f"  Cyc2 proteins:       {n_cyc2}"
)

print(
    f"  >=5-heme proteins:   {n_ge5}"
)

print(
    f"  Unique focal genes:  {n_priority}"
)


if n_cyc2 != EXPECTED_CYC2:

    fail(
        f"Expected {EXPECTED_CYC2} Cyc2 proteins "
        f"but found {n_cyc2}."
    )


if n_ge5 != EXPECTED_GE5_HEME:

    fail(
        f"Expected {EXPECTED_GE5_HEME} >=5-heme proteins "
        f"but found {n_ge5}."
    )


if n_priority != EXPECTED_PRIORITY_FOCALS:

    fail(
        f"Expected {EXPECTED_PRIORITY_FOCALS} unique "
        f"priority proteins but found {n_priority}."
    )


priority_keys = set(
    zip(
        priority[
            "genome"
        ],
        priority[
            "protein_id"
        ],
    )
)


## ================================================================== ##
## 4. Read fixed GlobDB cluster assignments
## ================================================================== ##

print()
print("Reading fixed GlobDB cluster assignments...")


globdb = pd.read_csv(
    GLOBDB_ASSIGNMENTS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_globdb = {
    "genome",
    "protein_id",
    "assignment_status",
    "n_clusters_hit",
    "top_scoring_cluster",
    "top_target_protein",
    "top_pident",
    "top_query_coverage",
    "top_target_coverage",
    "top_bits",
    "all_clusters_hit",
}


missing = (
    required_globdb
    - set(
        globdb.columns
    )
)


if missing:

    fail(
        "GlobDB assignment table is missing:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


if globdb[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys in "
        "GlobDB assignment table."
    )


globdb_keep = [
    "genome",
    "protein_id",
    "assignment_status",
    "n_passing_hits",
    "n_clusters_hit",
    "top_scoring_cluster",
    "top_target_protein",
    "top_cluster_representative",
    "top_pident",
    "top_query_coverage",
    "top_target_coverage",
    "top_evalue",
    "top_bits",
    "all_clusters_hit",
]


globdb_keep = [
    column

    for column
    in globdb_keep

    if column
    in globdb.columns
]


globdb = globdb[
    globdb_keep
].copy()


rename_globdb = {

    column:
        f"globdb_{column}"

    for column
    in globdb.columns

    if column
    not in {
        "genome",
        "protein_id",
    }
}


globdb = globdb.rename(
    columns=rename_globdb
)


## ================================================================== ##
## 5. Recover old MCL module of each top-scoring cluster
## ================================================================== ##

print()
print("Reading old GlobDB cluster -> MCL module mapping...")


old_pm = pd.read_csv(
    OLD_MODULE_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


for column in (
    "cluster",
    "module",
):

    if column not in old_pm.columns:

        fail(
            f"Old module membership table lacks "
            f"'{column}'."
        )


cluster_module = (
    old_pm[
        [
            "cluster",
            "module",
        ]
    ]
    .drop_duplicates()
)


module_counts = (
    cluster_module
    .groupby(
        "cluster"
    )[
        "module"
    ]
    .nunique()
)


bad_clusters = (
    module_counts[
        module_counts > 1
    ]
)


if len(
    bad_clusters
) > 0:

    fail(
        "One or more old MMseq clusters belong to "
        "multiple MCL modules."
    )


cluster_to_module = dict(
    zip(
        cluster_module[
            "cluster"
        ],
        cluster_module[
            "module"
        ],
    )
)


## ================================================================== ##
## 6. Locate the authoritative FeGenie protein FASTA for each MAG
## ================================================================== ##

print()
print("Locating FeGenie Prodigal FASTAs...")


all_faa = sorted(
    path

    for path
    in ORF_DIR.glob(
        "*-proteins.faa"
    )

    if (
        path.is_file()
        and
        path.stat().st_size > 0
    )
)


if not all_faa:

    fail(
        f"No non-empty *-proteins.faa files found in:\n"
        f"{ORF_DIR}"
    )


genome_to_faa = {}


for path in all_faa:

    genome = normalise_genome_filename(
        path
    )


    if genome not in genome_set:

        continue


    if genome in genome_to_faa:

        fail(
            f"More than one FeGenie protein FASTA "
            f"matches MAG:\n{genome}"
        )


    genome_to_faa[
        genome
    ] = path


missing_fastas = sorted(
    genome_set
    - set(
        genome_to_faa
    )
)


if missing_fastas:

    fail(
        "Methylococcales MAGs lacking an authoritative "
        "FeGenie protein FASTA:\n"
        + "\n".join(
            missing_fastas
        )
    )


if len(
    genome_to_faa
) != EXPECTED_MAGS:

    fail(
        f"Expected {EXPECTED_MAGS} Methylococcales "
        f"protein FASTAs but found "
        f"{len(genome_to_faa)}."
    )


print(
    f"  Protein FASTAs found: "
    f"{len(genome_to_faa)}"
)


## ================================================================== ##
## 7. Parse every predicted ORF in the 16 MAGs
## ================================================================== ##

print()
print("Parsing all Prodigal coordinates...")


coordinate_rows = []

genome_qc_rows = []


for number, genome in enumerate(
    genomes,
    start=1,
):

    faa = genome_to_faa[
        genome
    ]


    seen_ids = set()

    genome_rows = []


    with faa.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:

            if not line.startswith(
                ">"
            ):

                continue


            record = parse_header(
                line,
                genome,
                faa,
            )


            protein_id = record[
                "protein_id"
            ]


            if protein_id in seen_ids:

                fail(
                    f"Duplicate protein ID in "
                    f"{faa.name}:\n"
                    f"{protein_id}"
                )


            seen_ids.add(
                protein_id
            )


            genome_rows.append(
                record
            )


    if not genome_rows:

        fail(
            f"No Prodigal headers parsed from:\n"
            f"{faa}"
        )


    coordinate_rows.extend(
        genome_rows
    )


    genome_qc_rows.append(
        {
            "genome":
                genome,

            "n_predicted_proteins":
                len(
                    genome_rows
                ),

            "n_contigs_with_predicted_proteins":
                len(
                    {
                        row[
                            "contig"
                        ]
                        for row
                        in genome_rows
                    }
                ),

            "source_protein_fasta":
                str(
                    faa
                ),
        }
    )


    print(
        f"  Parsed {number:>2}/"
        f"{len(genomes)}: "
        f"{genome} "
        f"({len(genome_rows):,} proteins)"
    )


coords = pd.DataFrame(
    coordinate_rows
)


## ================================================================== ##
## 8. Coordinate QC and contig-relative gene ranks
## ================================================================== ##

print()
print("Validating coordinate catalogue...")


if coords[
    [
        "genome",
        "protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + protein_id keys in "
        "coordinate catalogue."
    )


if (
    pd.to_numeric(
        coords[
            "start"
        ],
        errors="raise",
    )
    < 1
).any():

    fail(
        "Invalid start coordinate detected."
    )


if (
    pd.to_numeric(
        coords[
            "end"
        ],
        errors="raise",
    )
    <
    pd.to_numeric(
        coords[
            "start"
        ],
        errors="raise",
    )
).any():

    fail(
        "Coordinate row with end < start detected."
    )


coords = (
    coords
    .sort_values(
        [
            "genome",
            "contig",
            "start",
            "end",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


coords[
    "contig_gene_rank"
] = (
    coords
    .groupby(
        [
            "genome",
            "contig",
        ],
        sort=False,
    )
    .cumcount()
    + 1
)


coords[
    "contig_gene_count"
] = (
    coords
    .groupby(
        [
            "genome",
            "contig",
        ],
        sort=False,
    )[
        "protein_id"
    ]
    .transform(
        "size"
    )
)


## Prodigal gene ordinal and coordinate-derived rank should normally agree.
## We record disagreements instead of silently assuming equivalence.

coords[
    "rank_matches_prodigal_ordinal"
] = (
    pd.to_numeric(
        coords[
            "contig_gene_rank"
        ],
        errors="raise",
    )
    ==
    pd.to_numeric(
        coords[
            "prodigal_gene_ordinal"
        ],
        errors="raise",
    )
).astype(int)


n_rank_disagreements = int(
    (
        coords[
            "rank_matches_prodigal_ordinal"
        ]
        == 0
    ).sum()
)


print(
    f"  Total predicted proteins: "
    f"{len(coords):,}"
)

print(
    f"  Gene-containing contigs:  "
    f"{coords[['genome','contig']].drop_duplicates().shape[0]:,}"
)

print(
    f"  Rank/ordinal disagreements:"
    f" {n_rank_disagreements:,}"
)


## ================================================================== ##
## 9. Validate every integrated candidate against coordinates
## ================================================================== ##

coord_keys = set(
    zip(
        coords[
            "genome"
        ],
        coords[
            "protein_id"
        ],
    )
)


integrated_keys = set(
    zip(
        integrated[
            "genome"
        ],
        integrated[
            "protein_id"
        ],
    )
)


unmatched_candidates = sorted(
    integrated_keys
    - coord_keys
)


if unmatched_candidates:

    unmatched_candidate_df = pd.DataFrame(
        unmatched_candidates,
        columns=[
            "genome",
            "protein_id",
        ],
    )

else:

    unmatched_candidate_df = pd.DataFrame(
        columns=[
            "genome",
            "protein_id",
        ],
    )


unmatched_candidate_df.to_csv(
    OUT_UNMATCHED_CANDIDATES,
    sep="\t",
    index=False,
)


if unmatched_candidates:

    fail(
        f"{len(unmatched_candidates)} integrated candidates "
        f"could not be mapped to the authoritative "
        f"Prodigal coordinate catalogue.\n"
        f"See:\n{OUT_UNMATCHED_CANDIDATES}"
    )


## ================================================================== ##
## 10. Build lean integrated gene catalogue
## ================================================================== ##

print()
print("Attaching functional evidence...")


## Keep evidence useful for neighborhood interpretation.
## Do NOT copy full DeepTMHMM sequence/topology strings into every table.

candidate_columns = [
    "genome",
    "protein_id",
    "length",
    "candidate_source",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",
    "fegenie_max_bitscore",
    "fegenie_max_cutoff",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "signalp_cs_position",

    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "deeptmhmm_sp_length",
    "deeptmhmm_side_after_sp",
    "deeptmhmm_c_terminal_side",

    "signalp_export_positive",
    "deeptmhmm_sp_positive",
    "export_evidence",
    "localization_class",

    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",

    "priority_cyc2",
    "priority_ge5_hemes",
    "priority_reason",
]


candidate_columns = [
    column

    for column
    in candidate_columns

    if column
    in integrated.columns
]


candidate_meta = integrated[
    candidate_columns
].copy()


catalog = coords.merge(
    candidate_meta,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


catalog = catalog.merge(
    globdb,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one",
)


## ================================================================== ##
## 11. Fill explicit evidence flags for ordinary genes
## ================================================================== ##

catalog[
    "is_integrated_candidate"
] = (
    catalog[
        "candidate_source"
    ]
    .fillna("")
    .astype(str)
    .str.strip()
    != ""
).astype(int)


flag_columns = [
    "fegenie_positive",
    "findmehemes_positive",
    "signalp_export_positive",
    "deeptmhmm_sp_positive",
    "export_evidence",
    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",
    "priority_cyc2",
    "priority_ge5_hemes",
]


for column in flag_columns:

    if column not in catalog.columns:

        catalog[
            column
        ] = 0

    catalog[
        column
    ] = numeric_flag(
        catalog[
            column
        ]
    )


catalog[
    "number_of_hemes"
] = (
    pd.to_numeric(
        catalog[
            "number_of_hemes"
        ],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
)


## Ordinary genes have no priority reason.
catalog[
    "priority_reason"
] = (
    catalog[
        "priority_reason"
    ]
    .fillna("")
)


## ================================================================== ##
## 12. Add top old MCL module assignment
## ================================================================== ##

catalog[
    "globdb_top_module"
] = (
    catalog[
        "globdb_top_scoring_cluster"
    ]
    .fillna("")
    .map(
        cluster_to_module
    )
    .fillna("")
)


## ================================================================== ##
## 13. Final ordering
## ================================================================== ##

first_columns = [
    "genome",
    "protein_id",
    "contig",
    "contig_gene_rank",
    "contig_gene_count",

    "start",
    "end",
    "gene_length_nt",
    "strand",

    "prodigal_id",
    "prodigal_seqnum",
    "prodigal_gene_ordinal",
    "rank_matches_prodigal_ordinal",

    "partial",
    "start_type",

    "is_integrated_candidate",

    "fegenie_positive",
    "fegenie_HMMs",
    "fegenie_categories",

    "findmehemes_positive",
    "number_of_hemes",

    "signalp_prediction",
    "deeptmhmm_class",
    "deeptmhmm_n_tm_helices",
    "localization_class",

    "globdb_assignment_status",
    "globdb_top_scoring_cluster",
    "globdb_top_module",
    "globdb_all_clusters_hit",
    "globdb_top_pident",
    "globdb_top_query_coverage",
    "globdb_top_target_coverage",

    "priority_cyc2",
    "priority_ge5_hemes",
    "priority_reason",
]


first_columns = [
    column

    for column
    in first_columns

    if column
    in catalog.columns
]


remaining_columns = [
    column

    for column
    in catalog.columns

    if column
    not in first_columns
]


catalog = catalog[
    first_columns
    + remaining_columns
]


catalog = (
    catalog
    .sort_values(
        [
            "genome",
            "contig",
            "start",
            "end",
            "protein_id",
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 14. Write structural coordinate table
## ================================================================== ##

coordinate_columns = [
    "genome",
    "protein_id",
    "contig",
    "contig_gene_rank",
    "contig_gene_count",

    "start",
    "end",
    "gene_length_nt",
    "strand",

    "prodigal_strand",
    "prodigal_id",
    "prodigal_seqnum",
    "prodigal_gene_ordinal",
    "rank_matches_prodigal_ordinal",

    "partial",
    "start_type",
    "rbs_motif",
    "rbs_spacer",
    "gc_cont",

    "source_protein_fasta",
]


catalog[
    coordinate_columns
].to_csv(
    OUT_COORDINATES,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 15. Write complete gene catalogue
## ================================================================== ##

catalog.to_csv(
    OUT_CATALOG,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 16. Extract and validate the 21 priority focal genes
## ================================================================== ##

focals = catalog[
    catalog[
        "priority_reason"
    ]
    != ""
].copy()


focal_keys = set(
    zip(
        focals[
            "genome"
        ],
        focals[
            "protein_id"
        ],
    )
)


unmatched_focals = sorted(
    priority_keys
    - focal_keys
)


if unmatched_focals:

    unmatched_focal_df = pd.DataFrame(
        unmatched_focals,
        columns=[
            "genome",
            "protein_id",
        ],
    )

else:

    unmatched_focal_df = pd.DataFrame(
        columns=[
            "genome",
            "protein_id",
        ],
    )


unmatched_focal_df.to_csv(
    OUT_UNMATCHED_FOCALS,
    sep="\t",
    index=False,
)


if unmatched_focals:

    fail(
        f"{len(unmatched_focals)} priority focal genes "
        f"could not be recovered from the complete catalogue.\n"
        f"See:\n{OUT_UNMATCHED_FOCALS}"
    )


if len(
    focals
) != EXPECTED_PRIORITY_FOCALS:

    fail(
        f"Final focal catalogue contains "
        f"{len(focals)} rows; expected "
        f"{EXPECTED_PRIORITY_FOCALS}."
    )


## Add MAG-level metadata to the small focal table.
## Prefix fields so they cannot collide with gene-level fields.

mag_for_join = mag_meta.copy()


mag_for_join = mag_for_join.rename(
    columns={
        column:
            (
                "genome"
                if column == "Genome_ID"
                else f"mag_{column}"
            )

        for column
        in mag_for_join.columns
    }
)


focals = focals.merge(
    mag_for_join,
    on="genome",
    how="left",
    validate="many_to_one",
)


focals = (
    focals
    .sort_values(
        [
            "priority_reason",
            "number_of_hemes",
            "genome",
            "protein_id",
        ],
        ascending=[
            True,
            False,
            True,
            True,
        ],
        kind="stable",
    )
    .reset_index(
        drop=True
    )
)


focals.to_csv(
    OUT_FOCALS,
    sep="\t",
    index=False,
    na_rep="",
)


## ================================================================== ##
## 17. Genome-level QC
## ================================================================== ##

genome_qc = pd.DataFrame(
    genome_qc_rows
)


candidate_counts = (
    catalog
    .groupby(
        "genome"
    )
    .agg(
        integrated_candidates=(
            "is_integrated_candidate",
            "sum",
        ),

        fegenie_positive=(
            "fegenie_positive",
            "sum",
        ),

        findmehemes_positive=(
            "findmehemes_positive",
            "sum",
        ),

        priority_focals=(
            "priority_reason",
            lambda values:
                int(
                    (
                        values
                        .astype(str)
                        .str.strip()
                        != ""
                    ).sum()
                ),
        ),
    )
    .reset_index()
)


genome_qc = genome_qc.merge(
    candidate_counts,
    on="genome",
    how="left",
    validate="one_to_one",
)


genome_qc = (
    genome_qc
    .sort_values(
        "genome"
    )
    .reset_index(
        drop=True
    )
)


genome_qc.to_csv(
    OUT_GENOME_QC,
    sep="\t",
    index=False,
)


## ================================================================== ##
## 18. Global QC summary
## ================================================================== ##

qc_rows = [
    (
        "methylococcales_MAGs",
        len(
            genomes
        ),
    ),

    (
        "protein_FASTAs",
        len(
            genome_to_faa
        ),
    ),

    (
        "predicted_proteins",
        len(
            catalog
        ),
    ),

    (
        "gene_contigs",
        catalog[
            [
                "genome",
                "contig",
            ]
        ]
        .drop_duplicates()
        .shape[0],
    ),

    (
        "integrated_candidates",
        int(
            catalog[
                "is_integrated_candidate"
            ].sum()
        ),
    ),

    (
        "fegenie_positive",
        int(
            catalog[
                "fegenie_positive"
            ].sum()
        ),
    ),

    (
        "findmehemes_positive",
        int(
            catalog[
                "findmehemes_positive"
            ].sum()
        ),
    ),

    (
        "priority_Cyc2",
        int(
            catalog[
                "priority_cyc2"
            ].sum()
        ),
    ),

    (
        "priority_ge5_hemes",
        int(
            catalog[
                "priority_ge5_hemes"
            ].sum()
        ),
    ),

    (
        "priority_focals_total",
        len(
            focals
        ),
    ),

    (
        "unmatched_integrated_candidates",
        len(
            unmatched_candidates
        ),
    ),

    (
        "unmatched_priority_focals",
        len(
            unmatched_focals
        ),
    ),

    (
        "coordinate_rank_ordinal_disagreements",
        n_rank_disagreements,
    ),
]


qc = pd.DataFrame(
    qc_rows,
    columns=[
        "metric",
        "value",
    ],
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False,
)


## ================================================================== ##
## Final report
## ================================================================== ##

print()
print("=" * 80)
print("STAGE 50 COMPLETE")
print("=" * 80)

print(
    f"Methylococcales MAGs:        "
    f"{len(genomes):,}"
)

print(
    f"Predicted proteins:          "
    f"{len(catalog):,}"
)

print(
    f"Integrated candidates:       "
    f"{int(catalog['is_integrated_candidate'].sum()):,}"
)

print(
    f"FeGenie-positive proteins:   "
    f"{int(catalog['fegenie_positive'].sum()):,}"
)

print(
    f"FindMeHemes-positive:        "
    f"{int(catalog['findmehemes_positive'].sum()):,}"
)

print(
    f"Cyc2 focal proteins:         "
    f"{int(catalog['priority_cyc2'].sum()):,}"
)

print(
    f">=5-heme focal proteins:     "
    f"{int(catalog['priority_ge5_hemes'].sum()):,}"
)

print(
    f"Priority focal genes:        "
    f"{len(focals):,}"
)

print(
    f"Rank/ordinal disagreements:  "
    f"{n_rank_disagreements:,}"
)

print()

print(
    f"Coordinates:\n  "
    f"{OUT_COORDINATES}"
)

print(
    f"Gene catalogue:\n  "
    f"{OUT_CATALOG}"
)

print(
    f"Priority focals:\n  "
    f"{OUT_FOCALS}"
)

print(
    f"Genome QC:\n  "
    f"{OUT_GENOME_QC}"
)

print(
    f"Global QC:\n  "
    f"{OUT_QC}"
)
