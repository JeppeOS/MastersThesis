#!/usr/bin/env python3

## ================================================================== ##
## GENERIC LOCAL NEIGHBORHOOD FAMILY + ARCHITECTURE PIPELINE
##
## Purpose
## -------
## Generalized replacement for module-specific local-family and
## architecture-profiling scripts.
##
## Input:
##   A focal-gene-aligned neighborhood table containing actual genes.
##
## Primary observational unit:
##   focal region, NOT genome.
##
## This is important because one genome may contain multiple focal
## paralogues and therefore multiple biologically distinct regions.
##
## Workflow:
##
##   1. Define focal regions.
##   2. Extract every unique protein occurring in those regions.
##   3. Cluster proteins locally with MMseqs2.
##   4. Assign module-specific local-family IDs.
##   5. Map family identity back to every neighborhood observation.
##   6. Profile family prevalence across regions and genomes.
##   7. Group exact presence patterns.
##   8. Calculate exhaustive pairwise Jaccard.
##   9. Calculate pairwise order / spacing / orientation conservation.
##  10. Generate ordered architecture signatures per focal region.
##
## Local-family IDs are MODULE-SPECIFIC.
##
## Example:
##   M35_local_008
##
## must NOT automatically be compared with:
##   M20_local_008
##
## Cross-module homology requires a separate sequence comparison or
## common cross-module family catalogue.
## ================================================================== ##


from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import itertools
import json
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd


## ================================================================== ##
## Helpers
## ================================================================== ##

def fail(message):
    print(
        f"ERROR: {message}",
        file=sys.stderr
    )
    sys.exit(1)


def clean(value):

    if value is None:
        return ""

    if pd.isna(value):
        return ""

    return str(value).strip()


def expand_path(value):

    return Path(
        str(value)
    ).expanduser().resolve()


def safe_sd(values):

    values = list(values)

    if len(values) <= 1:
        return 0.0

    return float(
        np.std(
            values,
            ddof=1
        )
    )


def dominant_nonempty(values):

    values = [
        clean(x)
        for x in values
        if clean(x) != ""
    ]

    if not values:
        return "", 0

    value, count = Counter(
        values
    ).most_common(1)[0]

    return value, count


def fasta_iter(path):

    header = None
    sequence = []

    with path.open("r") as handle:

        for raw in handle:

            line = raw.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    yield header, "".join(sequence)

                header = line[1:]
                sequence = []

            else:

                sequence.append(
                    line.strip()
                )

    if header is not None:
        yield header, "".join(sequence)


def prodigal_protein_id(header):

    return header.split(
        " # ",
        1
    )[0].strip()


def find_orf_fasta(
    orf_dir,
    genome
):

    candidates = [

        orf_dir
        /
        f"{genome}.fa-proteins.faa",

        orf_dir
        /
        f"{genome}-proteins.faa",

    ]

    existing = [
        path
        for path in candidates
        if path.exists()
    ]

    if len(existing) == 1:
        return existing[0]

    if len(existing) > 1:
        fail(
            f"Multiple exact ORF FASTA matches for {genome}: "
            +
            ", ".join(
                str(x)
                for x in existing
            )
        )


    globbed = list(
        orf_dir.glob(
            f"{genome}*-proteins.faa"
        )
    )

    if len(globbed) == 1:
        return globbed[0]

    if len(globbed) == 0:
        fail(
            f"No ORF FASTA found for genome {genome}"
        )

    fail(
        f"Multiple ORF FASTA matches for genome {genome}: "
        +
        ", ".join(
            str(x)
            for x in globbed
        )
    )


def detect_column(
    df,
    configured,
    candidates,
    description,
    required=True
):

    if configured:

        if configured not in df.columns:
            fail(
                f"Configured {description} column "
                f"'{configured}' was not found."
            )

        return configured


    for candidate in candidates:

        if candidate in df.columns:
            return candidate


    if required:

        fail(
            f"Could not determine {description} column. "
            f"Tried: {', '.join(candidates)}"
        )

    return None


def run_command(command):

    print()
    print(
        "$ "
        +
        " ".join(
            str(x)
            for x in command
        )
    )

    subprocess.run(
        [
            str(x)
            for x in command
        ],
        check=True
    )


## ================================================================== ##
## Arguments
## ================================================================== ##

parser = argparse.ArgumentParser()


parser.add_argument(
    "--config",
    required=True,
    help="JSON module configuration file"
)


parser.add_argument(
    "--force-mmseqs",
    action="store_true",
    help="Rerun MMseqs even if cluster TSV already exists"
)


args = parser.parse_args()


## ================================================================== ##
## Configuration
## ================================================================== ##

config_path = expand_path(
    args.config
)


if not config_path.exists():
    fail(
        f"Configuration does not exist: {config_path}"
    )


with config_path.open() as handle:
    config = json.load(handle)


required_config = [
    "module_label",
    "local_prefix",
    "input_table",
    "orf_dir",
    "output_dir",
]


for key in required_config:

    if key not in config:
        fail(
            f"Configuration is missing required key: {key}"
        )


MODULE = clean(
    config[
        "module_label"
    ]
)


LOCAL_PREFIX = clean(
    config[
        "local_prefix"
    ]
)


INPUT = expand_path(
    config[
        "input_table"
    ]
)


ORF_DIR = expand_path(
    config[
        "orf_dir"
    ]
)


OUTPUT_DIR = expand_path(
    config[
        "output_dir"
    ]
)

## ================================================================== ##
## Optional stable local-family reference
##
## If supplied, newly generated MMseqs clusters must match the
## reference protein-membership sets exactly.
##
## Existing local_family_id values are then reused.
## ================================================================== ##

REFERENCE_MEMBERSHIP = None


reference_value = config.get(
    "reference_membership_table"
)


if reference_value:

    REFERENCE_MEMBERSHIP = expand_path(
        reference_value
    )

    if not REFERENCE_MEMBERSHIP.exists():

        fail(
            "Configured reference_membership_table does not exist: "
            f"{REFERENCE_MEMBERSHIP}"
        )

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


if not INPUT.exists():
    fail(
        f"Input table not found: {INPUT}"
    )


if not ORF_DIR.exists():
    fail(
        f"ORF directory not found: {ORF_DIR}"
    )


## ================================================================== ##
## Output paths
## ================================================================== ##

SEQUENCE_FASTA = (
    OUTPUT_DIR
    /
    f"{MODULE}_neighborhood_proteins.faa"
)


UNIQUE_GENE_TABLE = (
    OUTPUT_DIR
    /
    f"{MODULE}_unique_neighborhood_genes.tsv"
)


MMSEQS_DIR = (
    OUTPUT_DIR
    /
    "mmseqs"
)


MMSEQS_DIR.mkdir(
    parents=True,
    exist_ok=True
)


MMSEQS_DB = (
    MMSEQS_DIR
    /
    "proteins_db"
)


MMSEQS_CLUSTER_DB = (
    MMSEQS_DIR
    /
    "proteins_cluster"
)


MMSEQS_TMP = (
    MMSEQS_DIR
    /
    "tmp"
)


MMSEQS_CLUSTER_TSV = (
    OUTPUT_DIR
    /
    f"{MODULE}_mmseqs_cluster.tsv"
)


MEMBERSHIP_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_local_family_membership.tsv"
)


FAMILY_SUMMARY_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_local_family_summary.tsv"
)


AUGMENTED_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_gene_map_data_with_local_families.tsv"
)


REGION_PRESENCE_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_region_family_presence_matrix.tsv"
)


GENOME_PRESENCE_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_genome_family_presence_matrix.tsv"
)


ARCHITECTURE_SUMMARY_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_local_family_architecture_summary.tsv"
)


REGION_PATTERNS_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_exact_region_presence_patterns.tsv"
)


GENOME_PATTERNS_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_exact_genome_presence_patterns.tsv"
)


JACCARD_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_pairwise_region_jaccard.tsv"
)


PAIR_ARCHITECTURE_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_pairwise_local_architecture.tsv"
)


SIGNATURE_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_region_architecture_signatures.tsv"
)


QC_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_local_architecture_qc.tsv"
)

REFERENCE_VALIDATION_OUT = (
    OUTPUT_DIR
    /
    f"{MODULE}_reference_family_validation.tsv"
)

## ================================================================== ##
## Read input
## ================================================================== ##

print("=" * 100)
print(
    f"LOCAL ARCHITECTURE PIPELINE: {MODULE}"
)
print("=" * 100)
print()


genes = pd.read_csv(
    INPUT,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


print(
    f"Input neighborhood rows: {len(genes):,}"
)


## ================================================================== ##
## Resolve column names
## ================================================================== ##

genome_col = detect_column(
    genes,
    config.get(
        "genome_column"
    ),
    [
        "genome"
    ],
    "genome"
)


gene_id_col = detect_column(
    genes,
    config.get(
        "gene_id_column"
    ),
    [
        "protein_id",
        "neighbor_protein_id"
    ],
    "gene ID"
)


taxonomy_col = detect_column(
    genes,
    config.get(
        "taxonomy_column"
    ),
    [
        "taxonomy_display"
    ],
    "taxonomy",
    required=False
)


start_col = detect_column(
    genes,
    config.get(
        "start_column"
    ),
    [
        "plot_start_bp"
    ],
    "normalized start coordinate"
)


end_col = detect_column(
    genes,
    config.get(
        "end_column"
    ),
    [
        "plot_end_bp"
    ],
    "normalized end coordinate"
)


strand_col = detect_column(
    genes,
    config.get(
        "strand_column"
    ),
    [
        "plot_strand"
    ],
    "normalized strand"
)


## ================================================================== ##
## Normalize basic fields
## ================================================================== ##

genes["genome"] = genes[
    genome_col
].map(
    clean
)


genes["protein_id"] = genes[
    gene_id_col
].map(
    clean
)


genes["plot_start_bp"] = pd.to_numeric(
    genes[
        start_col
    ],
    errors="raise"
)


genes["plot_end_bp"] = pd.to_numeric(
    genes[
        end_col
    ],
    errors="raise"
)


genes["plot_strand"] = genes[
    strand_col
].map(
    clean
)


genes["midpoint_bp"] = (
    genes[
        "plot_start_bp"
    ]
    +
    genes[
        "plot_end_bp"
    ]
) / 2


if taxonomy_col:

    genes["taxonomy_display_internal"] = genes[
        taxonomy_col
    ].map(
        clean
    )

else:

    genes["taxonomy_display_internal"] = genes[
        "genome"
    ]


## ================================================================== ##
## Define focal-region identity
##
## Priority:
##
## 1. Explicit region_id_column
## 2. Explicit/detected focal_id_column
## 3. Infer one focal protein per genome from focal-marker rows
##
## If several focal regions occur in the same genome and the table
## lacks focal-region identity, the script intentionally fails.
## ================================================================== ##

region_col = config.get(
    "region_id_column"
)


focal_col = config.get(
    "focal_id_column"
)


if region_col:

    if region_col not in genes.columns:
        fail(
            f"Configured region_id_column '{region_col}' not found."
        )

    genes["region_id"] = genes[
        region_col
    ].map(
        clean
    )


else:

    if not focal_col:

        for candidate in [
            "focal_protein_id",
            "focal_gene_id",
            "focal_id"
        ]:

            if candidate in genes.columns:
                focal_col = candidate
                break


    if focal_col:

        if focal_col not in genes.columns:
            fail(
                f"Configured focal_id_column '{focal_col}' not found."
            )

        genes["focal_protein_id_internal"] = genes[
            focal_col
        ].map(
            clean
        )

        genes["region_id"] = (
            genes[
                "genome"
            ]
            +
            "|"
            +
            genes[
                "focal_protein_id_internal"
            ]
        )


    else:

        marker_col = config.get(
            "focal_marker_column"
        )


        marker_values = {
            clean(x)
            for x in config.get(
                "focal_marker_values",
                []
            )
        }


        if (
            marker_col is None
            or
            marker_col not in genes.columns
            or
            not marker_values
        ):

            fail(
                "No region/focal ID column was available and focal "
                "rows cannot be inferred from the configuration."
            )


        focal_rows = genes[
            genes[
                marker_col
            ].map(
                clean
            ).isin(
                marker_values
            )
        ].copy()


        if focal_rows.empty:
            fail(
                "No focal rows matched configured focal marker."
            )


        focal_counts = (
            focal_rows
            .groupby(
                "genome"
            )[
                "protein_id"
            ]
            .nunique()
        )


        problematic = focal_counts[
            focal_counts != 1
        ]


        if len(problematic) > 0:

            fail(
                "Multiple focal proteins occur in at least one genome, "
                "but the input table does not contain an explicit focal "
                "or region ID column. For multi-paralogue modules, "
                "provide focal_id_column or region_id_column."
            )


        genome_to_focal = (
            focal_rows
            .drop_duplicates(
                [
                    "genome",
                    "protein_id"
                ]
            )
            .set_index(
                "genome"
            )[
                "protein_id"
            ]
            .to_dict()
        )


        missing_focal = sorted(
            set(
                genes[
                    "genome"
                ]
            )
            -
            set(
                genome_to_focal
            )
        )


        if missing_focal:

            fail(
                "No focal gene could be inferred for genomes: "
                +
                ", ".join(
                    missing_focal
                )
            )


        genes["focal_protein_id_internal"] = genes[
            "genome"
        ].map(
            genome_to_focal
        )


        genes["region_id"] = (
            genes[
                "genome"
            ]
            +
            "|"
            +
            genes[
                "focal_protein_id_internal"
            ]
        )


## ================================================================== ##
## Region QC
## ================================================================== ##

n_regions = genes[
    "region_id"
].nunique()


n_genomes = genes[
    "genome"
].nunique()


print(
    f"Genomes:       {n_genomes:,}"
)


print(
    f"Focal regions: {n_regions:,}"
)


## ================================================================== ##
## Unique protein set
##
## A gene may appear in overlapping focal regions. It must be
## clustered only once, then mapped back to every region observation.
## ================================================================== ##

unique_genes = (
    genes[
        [
            "genome",
            "protein_id"
        ]
    ]
    .drop_duplicates()
)


n_unique_genes = len(
    unique_genes
)


print(
    f"Unique genes:  {n_unique_genes:,}"
)


protein_genome_counts = (
    unique_genes
    .groupby(
        "protein_id"
    )[
        "genome"
    ]
    .nunique()
)


nonunique_ids = protein_genome_counts[
    protein_genome_counts > 1
]


if len(nonunique_ids) > 0:

    fail(
        "Protein IDs are not globally unique across genomes. "
        "Synthetic MMseqs IDs would be required."
    )


## ================================================================== ##
## Extract authoritative protein sequences
## ================================================================== ##

wanted_by_genome = defaultdict(
    set
)


for row in unique_genes.itertuples(
    index=False
):

    wanted_by_genome[
        row.genome
    ].add(
        row.protein_id
    )


found = {}


print()
print(
    "Extracting authoritative protein sequences..."
)


for genome in sorted(
    wanted_by_genome
):

    fasta = find_orf_fasta(
        ORF_DIR,
        genome
    )


    wanted = wanted_by_genome[
        genome
    ]


    for header, sequence in fasta_iter(
        fasta
    ):

        protein_id = prodigal_protein_id(
            header
        )

        if protein_id in wanted:

            found[
                protein_id
            ] = sequence


missing_sequences = sorted(
    set(
        unique_genes[
            "protein_id"
        ]
    )
    -
    set(
        found
    )
)


if missing_sequences:

    fail(
        f"{len(missing_sequences)} proteins were not found in "
        "authoritative ORF FASTAs. First examples: "
        +
        ", ".join(
            missing_sequences[:10]
        )
    )


print(
    f"Sequences found: {len(found):,}/{n_unique_genes:,}"
)


## ================================================================== ##
## Write clustering FASTA
## ================================================================== ##

with SEQUENCE_FASTA.open(
    "w"
) as handle:

    for protein_id in sorted(
        found
    ):

        handle.write(
            f">{protein_id}\n"
        )

        handle.write(
            found[
                protein_id
            ]
            +
            "\n"
        )


unique_gene_metadata = (
    unique_genes
    .copy()
)


unique_gene_metadata[
    "sequence_length_aa"
] = unique_gene_metadata[
    "protein_id"
].map(
    lambda p:
        len(
            found[p]
        )
)


unique_gene_metadata.to_csv(
    UNIQUE_GENE_TABLE,
    sep="\t",
    index=False
)


## ================================================================== ##
## MMseqs local clustering
## ================================================================== ##

if shutil.which(
    "mmseqs"
) is None:

    fail(
        "mmseqs is not available in PATH. Activate the MMseqs2 "
        "environment before running this pipeline."
    )


run_mmseqs = (
    args.force_mmseqs
    or
    not MMSEQS_CLUSTER_TSV.exists()
)


if run_mmseqs:

    if args.force_mmseqs:

        if MMSEQS_DIR.exists():

            shutil.rmtree(
                MMSEQS_DIR
            )

        MMSEQS_DIR.mkdir(
            parents=True,
            exist_ok=True
        )


    min_seq_id = float(
        config.get(
            "mmseqs_min_seq_id",
            0.40
        )
    )


    coverage = float(
        config.get(
            "mmseqs_coverage",
            0.80
        )
    )


    cov_mode = int(
        config.get(
            "mmseqs_cov_mode",
            0
        )
    )


    run_command(
        [
            "mmseqs",
            "createdb",
            SEQUENCE_FASTA,
            MMSEQS_DB
        ]
    )


    run_command(
        [
            "mmseqs",
            "cluster",
            MMSEQS_DB,
            MMSEQS_CLUSTER_DB,
            MMSEQS_TMP,
            "--min-seq-id",
            str(
                min_seq_id
            ),
            "-c",
            str(
                coverage
            ),
            "--cov-mode",
            str(
                cov_mode
            )
        ]
    )


    run_command(
        [
            "mmseqs",
            "createtsv",
            MMSEQS_DB,
            MMSEQS_DB,
            MMSEQS_CLUSTER_DB,
            MMSEQS_CLUSTER_TSV
        ]
    )


else:

    print()
    print(
        f"Reusing existing MMseqs cluster TSV: "
        f"{MMSEQS_CLUSTER_TSV}"
    )


## ================================================================== ##
## Read headerless MMseqs membership
## ================================================================== ##

cluster_pairs = []


with MMSEQS_CLUSTER_TSV.open(
    "r"
) as handle:

    reader = csv.reader(
        handle,
        delimiter="\t"
    )

    for line_number, fields in enumerate(
        reader,
        start=1
    ):

        if not fields:
            continue

        if len(fields) < 2:

            fail(
                f"Malformed MMseqs TSV line {line_number}"
            )

        representative = clean(
            fields[0]
        )

        member = clean(
            fields[1]
        )

        cluster_pairs.append(
            (
                representative,
                member
            )
        )


member_to_rep = {}


for representative, member in cluster_pairs:

    if member in member_to_rep:

        fail(
            f"Protein occurs more than once in MMseqs output: {member}"
        )

    member_to_rep[
        member
    ] = representative


expected_proteins = set(
    unique_genes[
        "protein_id"
    ]
)


assigned_proteins = set(
    member_to_rep
)


if expected_proteins != assigned_proteins:

    fail(
        "MMseqs membership does not exactly match input protein set."
    )


print(
    f"MMseqs assignments: "
    f"{len(assigned_proteins):,}/{n_unique_genes:,}"
)


## ================================================================== ##
## Assign stable module-local family IDs
##
## CRITICAL DESIGN RULE:
##
## The analytical identity of a local family is its SET OF MEMBER
## PROTEINS.
##
## The protein chosen by MMseqs as cluster representative is metadata
## only and must NOT determine analytical family identity.
##
## Two modes are supported:
##
## 1. REFERENCE MODE
##
##    If reference_membership_table is provided:
##
##      - construct the member set of every new MMseqs cluster
##      - construct the member set of every reference family
##      - require exact one-to-one set identity
##      - reuse the established local_family_id
##
##    Any changed family membership causes a hard failure.
##
## 2. FIRST-RUN MODE
##
##    If no reference exists:
##
##      - assign IDs deterministically from family support and the
##        canonical sorted member IDs
##
##    The resulting membership table can then be locked as the
##    reference for future reruns.
## ================================================================== ##


## ------------------------------------------------------------------ ##
## Map MMseq representative -> exact set of member proteins
## ------------------------------------------------------------------ ##

rep_to_members = defaultdict(
    set
)


for member, representative in member_to_rep.items():

    rep_to_members[
        representative
    ].add(
        member
    )


## ------------------------------------------------------------------ ##
## Map every neighborhood observation to its current MMseq
## representative.
## ------------------------------------------------------------------ ##

genes[
    "local_family_representative"
] = genes[
    "protein_id"
].map(
    member_to_rep
)


if genes[
    "local_family_representative"
].isna().any():

    fail(
        "At least one neighborhood protein could not be mapped "
        "to an MMseqs cluster representative."
    )


## ------------------------------------------------------------------ ##
## Calculate support for each newly generated MMseq cluster.
##
## These values describe the new cluster but DO NOT determine its
## identity when a validated reference is available.
## ------------------------------------------------------------------ ##

family_stats = (
    genes[
        [
            "local_family_representative",
            "protein_id",
            "genome",
            "region_id"
        ]
    ]
    .drop_duplicates()
    .groupby(
        "local_family_representative",
        as_index=False
    )
    .agg(

        n_proteins=(
            "protein_id",
            "nunique"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

        n_regions=(
            "region_id",
            "nunique"
        )

    )
)


## ================================================================== ##
## REFERENCE MODE
## ================================================================== ##

reference_validation_rows = []


if REFERENCE_MEMBERSHIP is not None:

    print()
    print("=" * 100)
    print("VALIDATING LOCAL FAMILIES AGAINST REFERENCE MEMBERSHIP")
    print("=" * 100)

    print(
        f"Reference: {REFERENCE_MEMBERSHIP}"
    )


    reference = pd.read_csv(
        REFERENCE_MEMBERSHIP,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )


    reference_required = [
        "protein_id",
        "local_family_id"
    ]


    reference_missing = [
        col
        for col in reference_required
        if col not in reference.columns
    ]


    if reference_missing:

        fail(
            "Reference membership table is missing required columns: "
            +
            ", ".join(
                reference_missing
            )
        )


    reference[
        "protein_id"
    ] = reference[
        "protein_id"
    ].map(
        clean
    )


    reference[
        "local_family_id"
    ] = reference[
        "local_family_id"
    ].map(
        clean
    )


    reference = reference[
        (
            reference[
                "protein_id"
            ]
            !=
            ""
        )
        &
        (
            reference[
                "local_family_id"
            ]
            !=
            ""
        )
    ].copy()


    ## -------------------------------------------------------------- ##
    ## A protein may appear multiple times in a neighborhood table
    ## because overlapping focal regions can contain the same gene.
    ##
    ## Repeated observations are fine as long as that protein always
    ## belongs to the same local family.
    ## -------------------------------------------------------------- ##

    reference_protein_family = (
        reference[
            [
                "protein_id",
                "local_family_id"
            ]
        ]
        .drop_duplicates()
    )


    conflicting_reference = (
        reference_protein_family
        .groupby(
            "protein_id"
        )[
            "local_family_id"
        ]
        .nunique()
    )


    conflicting_reference = conflicting_reference[
        conflicting_reference > 1
    ]


    if len(
        conflicting_reference
    ) > 0:

        fail(
            f"{len(conflicting_reference)} proteins are assigned "
            "to multiple local families in the reference table."
        )


    ## -------------------------------------------------------------- ##
    ## The protein universe must be exactly identical.
    ## -------------------------------------------------------------- ##

    reference_proteins = set(
        reference_protein_family[
            "protein_id"
        ]
    )


    new_proteins = set(
        member_to_rep
    )


    missing_from_new = sorted(
        reference_proteins
        -
        new_proteins
    )


    new_not_in_reference = sorted(
        new_proteins
        -
        reference_proteins
    )


    if missing_from_new:

        fail(
            f"{len(missing_from_new)} reference proteins are absent "
            "from the new clustering. First examples: "
            +
            ", ".join(
                missing_from_new[:10]
            )
        )


    if new_not_in_reference:

        fail(
            f"{len(new_not_in_reference)} newly clustered proteins "
            "are absent from the reference. First examples: "
            +
            ", ".join(
                new_not_in_reference[:10]
            )
        )


    print(
        f"Protein universe: "
        f"{len(new_proteins):,}/{len(reference_proteins):,} exact"
    )


    ## -------------------------------------------------------------- ##
    ## Build reference family -> exact membership set
    ## -------------------------------------------------------------- ##

    reference_family_to_members = {

        family_id:
            frozenset(
                group[
                    "protein_id"
                ]
            )

        for family_id, group in (
            reference_protein_family
            .groupby(
                "local_family_id"
            )
        )

    }


    ## -------------------------------------------------------------- ##
    ## Build membership-signature -> reference family ID
    ##
    ## Two different family IDs must never describe exactly the same
    ## protein set.
    ## -------------------------------------------------------------- ##

    reference_signature_to_family = {}


    for family_id, members in reference_family_to_members.items():

        if members in reference_signature_to_family:

            fail(
                "Reference contains two local-family IDs with the "
                "same exact protein membership: "
                f"{reference_signature_to_family[members]} and "
                f"{family_id}"
            )

        reference_signature_to_family[
            members
        ] = family_id


    print(
        f"Reference families: "
        f"{len(reference_family_to_members):,}"
    )


    print(
        f"New MMseqs families: "
        f"{len(rep_to_members):,}"
    )


    if (
        len(
            reference_family_to_members
        )
        !=
        len(
            rep_to_members
        )
    ):

        fail(
            "Reference/new family counts differ: "
            f"{len(reference_family_to_members)} reference versus "
            f"{len(rep_to_members)} new."
        )


    ## -------------------------------------------------------------- ##
    ## Exact set matching
    ## -------------------------------------------------------------- ##

    rep_to_family = {}

    unmatched_new = []


    for representative, members_raw in rep_to_members.items():

        members = frozenset(
            members_raw
        )


        if members not in reference_signature_to_family:

            unmatched_new.append(
                (
                    representative,
                    sorted(
                        members
                    )
                )
            )

            continue


        stable_family_id = (
            reference_signature_to_family[
                members
            ]
        )


        rep_to_family[
            representative
        ] = stable_family_id


        reference_validation_rows.append(
            {

                "local_family_id":
                    stable_family_id,

                "new_mmseq_representative":
                    representative,

                "n_proteins":
                    len(
                        members
                    ),

                "membership_match":
                    "exact",

                "canonical_first_member":
                    sorted(
                        members
                    )[0],

            }
        )


    if unmatched_new:

        print()
        print(
            "The following new MMseqs families do not exactly match "
            "any reference family:"
        )

        for representative, members in unmatched_new[:10]:

            print()
            print(
                f"Representative: {representative}"
            )

            print(
                f"Members ({len(members)}): "
                +
                ", ".join(
                    members[:10]
                )
            )


        fail(
            f"{len(unmatched_new)} newly generated families failed "
            "exact reference-membership matching."
        )


    matched_reference_ids = set(
        rep_to_family.values()
    )


    missing_reference_families = sorted(
        set(
            reference_family_to_members
        )
        -
        matched_reference_ids
    )


    if missing_reference_families:

        fail(
            f"{len(missing_reference_families)} reference families "
            "were not recovered. First examples: "
            +
            ", ".join(
                missing_reference_families[:10]
            )
        )


    print()
    print(
        f"Exact family matches: "
        f"{len(rep_to_family):,}/"
        f"{len(reference_family_to_members):,}"
    )


    print(
        "Changed families:     0"
    )


    reference_validation = pd.DataFrame(
        reference_validation_rows
    )


    reference_validation = reference_validation.sort_values(
        "local_family_id"
    )


    reference_validation.to_csv(
        REFERENCE_VALIDATION_OUT,
        sep="\t",
        index=False
    )


## ================================================================== ##
## FIRST-RUN MODE
## ================================================================== ##

else:

    print()
    print(
        "No reference membership table supplied."
    )

    print(
        "Assigning deterministic first-run local-family IDs."
    )


    ## -------------------------------------------------------------- ##
    ## Add a canonical membership key that is independent of which
    ## protein MMseqs happened to choose as representative.
    ## -------------------------------------------------------------- ##

    family_stats[
        "canonical_first_member"
    ] = family_stats[
        "local_family_representative"
    ].map(
        lambda rep:
            sorted(
                rep_to_members[
                    rep
                ]
            )[0]
    )


    family_stats[
        "canonical_membership_key"
    ] = family_stats[
        "local_family_representative"
    ].map(
        lambda rep:
            "|".join(
                sorted(
                    rep_to_members[
                        rep
                    ]
                )
            )
    )


    ## -------------------------------------------------------------- ##
    ## Ranking uses biological support first, then canonical member
    ## identity as deterministic tie-breaker.
    ##
    ## It does NOT use the MMseq representative as a tie-breaker.
    ## -------------------------------------------------------------- ##

    family_stats = family_stats.sort_values(
        [
            "n_regions",
            "n_genomes",
            "n_proteins",
            "canonical_first_member",
            "canonical_membership_key"
        ],
        ascending=[
            False,
            False,
            False,
            True,
            True
        ],
        kind="stable"
    )


    rep_to_family = {}


    for index, representative in enumerate(
        family_stats[
            "local_family_representative"
        ],
        start=1
    ):

        rep_to_family[
            representative
        ] = (
            f"{LOCAL_PREFIX}_{index:03d}"
        )


## ================================================================== ##
## Apply stable family IDs
## ================================================================== ##

genes[
    "local_family_id"
] = genes[
    "local_family_representative"
].map(
    rep_to_family
)


if genes[
    "local_family_id"
].isna().any():

    fail(
        "At least one neighborhood gene was not assigned a stable "
        "local_family_id."
    )


family_stats[
    "local_family_id"
] = family_stats[
    "local_family_representative"
].map(
    rep_to_family
)


if family_stats[
    "local_family_id"
].isna().any():

    fail(
        "At least one MMseq family was not assigned a stable "
        "local_family_id."
    )


## ------------------------------------------------------------------ ##
## Restore useful output order by stable family ID.
##
## This affects presentation only, not identity.
## ------------------------------------------------------------------ ##

family_stats = family_stats.sort_values(
    "local_family_id",
    kind="stable"
)


n_local_families = len(
    family_stats
)


print()
print(
    f"Stable local families: "
    f"{n_local_families:,}"
)

## ================================================================== ##
## Add local-family support fields to every neighborhood observation
## ================================================================== ##

support_map = (
    family_stats
    .set_index(
        "local_family_id"
    )[
        [
            "n_proteins",
            "n_genomes",
            "n_regions"
        ]
    ]
    .to_dict(
        orient="index"
    )
)


genes[
    "local_family_n_proteins"
] = genes[
    "local_family_id"
].map(
    lambda x:
        support_map[x][
            "n_proteins"
        ]
)


genes[
    "local_family_n_genomes"
] = genes[
    "local_family_id"
].map(
    lambda x:
        support_map[x][
            "n_genomes"
        ]
)


genes[
    "local_family_n_regions"
] = genes[
    "local_family_id"
].map(
    lambda x:
        support_map[x][
            "n_regions"
        ]
)


## ================================================================== ##
## Save augmented neighborhood table
## ================================================================== ##

genes.to_csv(
    AUGMENTED_OUT,
    sep="\t",
    index=False
)


genes.to_csv(
    MEMBERSHIP_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Local-family annotation summary
## ================================================================== ##

annotation_fields = {

    "display_class":
        "display_class",

    "neighbor_globdb_cog":
        "dominant_globdb_cog",

    "neighbor_globdb_product":
        "dominant_globdb_product",

    "neighbor_fegenie_HMMs":
        "dominant_fegenie_HMM",

}


summary_rows = []


for family_id, g in genes.groupby(
    "local_family_id"
):

    row = {

        "local_family_id":
            family_id,

        "representative_protein":
            g[
                "local_family_representative"
            ].iloc[0],

        "n_proteins":
            g[
                "protein_id"
            ].nunique(),

        "n_regions":
            g[
                "region_id"
            ].nunique(),

        "n_genomes":
            g[
                "genome"
            ].nunique(),

    }


    for source_col, output_col in annotation_fields.items():

        if source_col in g.columns:

            value, support = dominant_nonempty(
                g[
                    source_col
                ]
            )

            row[
                output_col
            ] = value

            row[
                output_col
                +
                "_support"
            ] = support


    summary_rows.append(
        row
    )


family_summary = pd.DataFrame(
    summary_rows
)


family_summary = family_summary.merge(
    family_stats[
        [
            "local_family_id",
            "local_family_representative"
        ]
    ],
    on="local_family_id",
    how="left"
)


family_summary.to_csv(
    FAMILY_SUMMARY_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Presence matrices
##
## Region-level presence is PRIMARY.
##
## Genome-level presence is secondary because a genome may contain
## several focal regions.
## ================================================================== ##

region_presence = (
    genes[
        [
            "local_family_id",
            "region_id"
        ]
    ]
    .drop_duplicates()
    .assign(
        present=1
    )
    .pivot(
        index="local_family_id",
        columns="region_id",
        values="present"
    )
    .fillna(
        0
    )
    .astype(
        int
    )
)


region_presence.to_csv(
    REGION_PRESENCE_OUT,
    sep="\t"
)


genome_presence = (
    genes[
        [
            "local_family_id",
            "genome"
        ]
    ]
    .drop_duplicates()
    .assign(
        present=1
    )
    .pivot(
        index="local_family_id",
        columns="genome",
        values="present"
    )
    .fillna(
        0
    )
    .astype(
        int
    )
)


genome_presence.to_csv(
    GENOME_PRESENCE_OUT,
    sep="\t"
)


## ================================================================== ##
## Architecture summary
## ================================================================== ##

region_ids = sorted(
    genes[
        "region_id"
    ].unique()
)


genome_ids = sorted(
    genes[
        "genome"
    ].unique()
)


architecture_rows = []


for family_id, g in genes.groupby(
    "local_family_id"
):

    unique_region_family = g.drop_duplicates(
        [
            "region_id",
            "protein_id"
        ]
    )


    present_regions = set(
        unique_region_family[
            "region_id"
        ]
    )


    present_genomes = set(
        unique_region_family[
            "genome"
        ]
    )


    region_pattern = "".join(
        "1"
        if region in present_regions
        else "0"
        for region in region_ids
    )


    genome_pattern = "".join(
        "1"
        if genome in present_genomes
        else "0"
        for genome in genome_ids
    )


    midpoint_values = unique_region_family[
        "midpoint_bp"
    ].astype(
        float
    )


    strand_counts = Counter(
        unique_region_family[
            "plot_strand"
        ]
    )


    n_plus = strand_counts.get(
        "+",
        0
    )


    n_minus = strand_counts.get(
        "-",
        0
    )


    denominator = (
        n_plus
        +
        n_minus
    )


    strand_conservation = (
        max(
            n_plus,
            n_minus
        )
        /
        denominator
        if denominator > 0
        else np.nan
    )


    dominant_strand = ""

    if n_plus > n_minus:
        dominant_strand = "+"

    elif n_minus > n_plus:
        dominant_strand = "-"

    elif denominator > 0:
        dominant_strand = "mixed"


    architecture_rows.append(
        {

            "local_family_id":
                family_id,

            "n_proteins":
                unique_region_family[
                    "protein_id"
                ].nunique(),

            "n_regions":
                len(
                    present_regions
                ),

            "n_genomes":
                len(
                    present_genomes
                ),

            "region_fraction":
                len(
                    present_regions
                )
                /
                len(
                    region_ids
                ),

            "genome_fraction":
                len(
                    present_genomes
                )
                /
                len(
                    genome_ids
                ),

            "region_presence_pattern":
                region_pattern,

            "genome_presence_pattern":
                genome_pattern,

            "median_midpoint_bp":
                float(
                    midpoint_values.median()
                ),

            "min_midpoint_bp":
                float(
                    midpoint_values.min()
                ),

            "max_midpoint_bp":
                float(
                    midpoint_values.max()
                ),

            "position_span_bp":
                float(
                    midpoint_values.max()
                    -
                    midpoint_values.min()
                ),

            "position_sd_bp":
                safe_sd(
                    midpoint_values
                ),

            "dominant_strand":
                dominant_strand,

            "strand_conservation_fraction":
                strand_conservation,

            "n_plus":
                n_plus,

            "n_minus":
                n_minus,

        }
    )


architecture_summary = pd.DataFrame(
    architecture_rows
)


architecture_summary = architecture_summary.merge(
    family_summary,
    on=[
        "local_family_id",
        "n_proteins",
        "n_regions",
        "n_genomes"
    ],
    how="left"
)


architecture_summary = architecture_summary.sort_values(
    [
        "n_regions",
        "position_span_bp",
        "median_midpoint_bp",
        "local_family_id"
    ],
    ascending=[
        False,
        True,
        True,
        True
    ]
)


architecture_summary.to_csv(
    ARCHITECTURE_SUMMARY_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Exact presence-pattern summaries
## ================================================================== ##

def make_pattern_table(
    architecture_df,
    pattern_column
):

    rows = []


    for index, (
        pattern,
        group
    ) in enumerate(
        architecture_df.groupby(
            pattern_column
        ),
        start=1
    ):

        rows.append(
            {

                "pattern_id":
                    f"pattern_{index:03d}",

                "presence_pattern":
                    pattern,

                "n_families":
                    len(
                        group
                    ),

                "local_families":
                    "; ".join(
                        sorted(
                            group[
                                "local_family_id"
                            ]
                        )
                    ),

                "max_n_regions":
                    group[
                        "n_regions"
                    ].max(),

                "max_n_genomes":
                    group[
                        "n_genomes"
                    ].max(),

            }
        )


    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "max_n_regions",
            "n_families"
        ],
        ascending=[
            False,
            False
        ]
    )


make_pattern_table(
    architecture_summary,
    "region_presence_pattern"
).to_csv(
    REGION_PATTERNS_OUT,
    sep="\t",
    index=False
)


make_pattern_table(
    architecture_summary,
    "genome_presence_pattern"
).to_csv(
    GENOME_PATTERNS_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Exhaustive region-level pairwise Jaccard
##
## Region level is used because local architecture belongs to focal
## regions, not whole genomes.
## ================================================================== ##

family_to_regions = {

    family:
        set(
            genes.loc[
                genes[
                    "local_family_id"
                ]
                ==
                family,
                "region_id"
            ]
        )

    for family in sorted(
        genes[
            "local_family_id"
        ].unique()
    )
}


pair_rows = []


families = sorted(
    family_to_regions
)


for family_a, family_b in itertools.combinations(
    families,
    2
):

    regions_a = family_to_regions[
        family_a
    ]


    regions_b = family_to_regions[
        family_b
    ]


    shared = (
        regions_a
        &
        regions_b
    )


    union = (
        regions_a
        |
        regions_b
    )


    n_a = len(
        regions_a
    )


    n_b = len(
        regions_b
    )


    n_shared = len(
        shared
    )


    n_union = len(
        union
    )


    pair_rows.append(
        {

            "family_a":
                family_a,

            "family_b":
                family_b,

            "n_regions_a":
                n_a,

            "n_regions_b":
                n_b,

            "n_shared_regions":
                n_shared,

            "n_union_regions":
                n_union,

            "jaccard":
                (
                    n_shared
                    /
                    n_union
                    if n_union > 0
                    else np.nan
                ),

            "fraction_a_regions_with_b":
                (
                    n_shared
                    /
                    n_a
                    if n_a > 0
                    else np.nan
                ),

            "fraction_b_regions_with_a":
                (
                    n_shared
                    /
                    n_b
                    if n_b > 0
                    else np.nan
                ),

            "exact_same_region_pattern":
                int(
                    regions_a
                    ==
                    regions_b
                ),

            "shared_regions":
                "; ".join(
                    sorted(
                        shared
                    )
                ),

        }
    )


pairwise = pd.DataFrame(
    pair_rows
)


pairwise = pairwise.sort_values(
    [
        "jaccard",
        "n_shared_regions",
        "family_a",
        "family_b"
    ],
    ascending=[
        False,
        False,
        True,
        True
    ]
)


pairwise.to_csv(
    JACCARD_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Pairwise order / spacing / orientation
##
## Every actual paralogue combination is retained.
##
## We report both:
##
##   n_shared_regions
##   n_pair_observations
##
## so paralogues cannot silently inflate support.
## ================================================================== ##

pair_architecture_rows = []


for family_a, family_b in itertools.combinations(
    families,
    2
):

    a = genes[
        genes[
            "local_family_id"
        ]
        ==
        family_a
    ][
        [
            "region_id",
            "protein_id",
            "midpoint_bp",
            "plot_strand"
        ]
    ].drop_duplicates()


    b = genes[
        genes[
            "local_family_id"
        ]
        ==
        family_b
    ][
        [
            "region_id",
            "protein_id",
            "midpoint_bp",
            "plot_strand"
        ]
    ].drop_duplicates()


    ab = a.merge(
        b,
        on="region_id",
        how="inner",
        suffixes=(
            "_a",
            "_b"
        )
    )


    if ab.empty:
        continue


    shared_regions = sorted(
        ab[
            "region_id"
        ].unique()
    )


    if len(
        shared_regions
    ) < 2:
        continue


    deltas = (
        ab[
            "midpoint_bp_b"
        ]
        -
        ab[
            "midpoint_bp_a"
        ]
    )


    n_a_before_b = int(
        (
            deltas > 0
        ).sum()
    )


    n_b_before_a = int(
        (
            deltas < 0
        ).sum()
    )


    ordered_count = (
        n_a_before_b
        +
        n_b_before_a
    )


    order_conservation = (
        max(
            n_a_before_b,
            n_b_before_a
        )
        /
        ordered_count
        if ordered_count > 0
        else np.nan
    )


    if n_a_before_b > n_b_before_a:

        dominant_order = "A_before_B"

    elif n_b_before_a > n_a_before_b:

        dominant_order = "B_before_A"

    else:

        dominant_order = "mixed"


    orientation = (
        ab[
            "plot_strand_a"
        ]
        +
        "/"
        +
        ab[
            "plot_strand_b"
        ]
    )


    orientation_counts = Counter(
        orientation
    )


    dominant_orientation, orientation_support = (
        orientation_counts.most_common(
            1
        )[0]
    )


    pair_architecture_rows.append(
        {

            "family_a":
                family_a,

            "family_b":
                family_b,

            "n_shared_regions":
                len(
                    shared_regions
                ),

            "n_pair_observations":
                len(
                    ab
                ),

            "shared_regions":
                "; ".join(
                    shared_regions
                ),

            "median_delta_bp_B_minus_A":
                float(
                    np.median(
                        deltas
                    )
                ),

            "min_delta_bp_B_minus_A":
                float(
                    np.min(
                        deltas
                    )
                ),

            "max_delta_bp_B_minus_A":
                float(
                    np.max(
                        deltas
                    )
                ),

            "delta_span_bp":
                float(
                    np.max(
                        deltas
                    )
                    -
                    np.min(
                        deltas
                    )
                ),

            "delta_sd_bp":
                safe_sd(
                    deltas
                ),

            "dominant_order":
                dominant_order,

            "order_conservation_fraction":
                order_conservation,

            "n_A_before_B":
                n_a_before_b,

            "n_B_before_A":
                n_b_before_a,

            "dominant_relative_orientation":
                dominant_orientation,

            "orientation_conservation_fraction":
                (
                    orientation_support
                    /
                    len(
                        orientation
                    )
                ),

        }
    )


pair_architecture = pd.DataFrame(
    pair_architecture_rows
)


if not pair_architecture.empty:

    pair_architecture = pair_architecture.sort_values(
        [
            "n_shared_regions",
            "order_conservation_fraction",
            "delta_span_bp",
            "family_a",
            "family_b"
        ],
        ascending=[
            False,
            False,
            True,
            True,
            True
        ]
    )


pair_architecture.to_csv(
    PAIR_ARCHITECTURE_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Ordered architecture signature per focal region
##
## This will later be useful for representative-region selection.
## ================================================================== ##

signature_rows = []


for region_id, g in genes.groupby(
    "region_id"
):

    ordered = (
        g[
            [
                "protein_id",
                "local_family_id",
                "midpoint_bp",
                "plot_strand",
                "genome",
                "taxonomy_display_internal"
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "midpoint_bp",
                "protein_id"
            ]
        )
    )


    tokens = []


    for row in ordered.itertuples(
        index=False
    ):

        strand = (
            ">"
            if row.plot_strand == "+"
            else
            "<"
        )


        tokens.append(
            f"{strand}{row.local_family_id}"
        )


    signature_rows.append(
        {

            "region_id":
                region_id,

            "genome":
                ordered[
                    "genome"
                ].iloc[0],

            "taxonomy_display":
                ordered[
                    "taxonomy_display_internal"
                ].iloc[0],

            "n_genes":
                len(
                    ordered
                ),

            "n_local_families":
                ordered[
                    "local_family_id"
                ].nunique(),

            "architecture_signature":
                " | ".join(
                    tokens
                ),

        }
    )


signatures = pd.DataFrame(
    signature_rows
)


signatures.to_csv(
    SIGNATURE_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Regression / expected-value QC
## ================================================================== ##

expected = config.get(
    "expected",
    {}
)


observed = {

    "n_regions":
        n_regions,

    "n_unique_genes":
        n_unique_genes,

    "n_local_families":
        n_local_families,

}


for key, expected_value in expected.items():

    if key not in observed:
        continue

    if int(
        observed[
            key
        ]
    ) != int(
        expected_value
    ):

        fail(
            f"Regression QC failed for {key}: "
            f"expected {expected_value}, "
            f"observed {observed[key]}"
        )


## ================================================================== ##
## QC output
## ================================================================== ##

qc_rows = [
    {
        "metric":
            "genomes",
        "value":
            n_genomes
    },
    {
        "metric":
            "focal_regions",
        "value":
            n_regions
    },
    {
        "metric":
            "unique_genes",
        "value":
            n_unique_genes
    },
    {
        "metric":
            "local_families",
        "value":
            n_local_families
    },
    {
        "metric":
            "pairwise_family_comparisons",
        "value":
            len(
                pairwise
            )
    },
    {
        "metric":
            "pairs_copresent_ge2_regions",
        "value":
            len(
                pair_architecture
            )
    },
]


if REFERENCE_MEMBERSHIP is not None:

    qc_rows.extend(
        [
            {
                "metric":
                    "reference_validation_used",
                "value":
                    1
            },
            {
                "metric":
                    "reference_families",
                "value":
                    len(
                        reference_family_to_members
                    )
            },
            {
                "metric":
                    "reference_exact_family_matches",
                "value":
                    len(
                        rep_to_family
                    )
            },
            {
                "metric":
                    "reference_changed_families",
                "value":
                    0
            },
        ]
    )


qc = pd.DataFrame(
    qc_rows
)


qc.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)


## ================================================================== ##
## Console summary
## ================================================================== ##

print()
print("=" * 100)
print("SUCCESS")
print("=" * 100)


print(
    f"Module:         {MODULE}"
)


print(
    f"Genomes:        {n_genomes}"
)


print(
    f"Focal regions:  {n_regions}"
)


print(
    f"Unique genes:   {n_unique_genes}"
)


print(
    f"Local families: {n_local_families}"
)


print(
    f"Family pairs:   {len(pairwise):,}"
)


print(
    f"Pairs in >=2 regions: "
    f"{len(pair_architecture):,}"
)

if REFERENCE_MEMBERSHIP is not None:

    print()
    print(
        "Reference-family validation:"
    )

    print(
        f"  Reference families: "
        f"{len(reference_family_to_members)}"
    )

    print(
        f"  Exact matches:       "
        f"{len(rep_to_family)}/"
        f"{len(reference_family_to_members)}"
    )

    print(
        "  Changed families:    0"
    )

print()
print(
    "Most prevalent families:"
)


print(
    architecture_summary[
        [
            "local_family_id",
            "n_regions",
            "n_genomes",
            "median_midpoint_bp",
            "position_span_bp",
            "dominant_strand",
            "strand_conservation_fraction"
        ]
    ]
    .head(
        20
    )
    .to_string(
        index=False
    )
)


print()
print(
    f"Outputs written to:\n{OUTPUT_DIR}"
)
