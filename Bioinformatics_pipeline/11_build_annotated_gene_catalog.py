#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent


PRODIGAL_COORDS = (
    WORKFLOW
    / "10_genomic_coordinates"
    / "all_prodigal_gene_coordinates.tsv"
)


ANNOTATION_MAP = (
    HERE
    / "prodigal_globdb_annotation_map.tsv"
)


CLUSTER_MEMBERSHIP = (
    WORKFLOW
    / "06_clustering"
    / "cluster_membership.tsv"
)


MODULE_MEMBERSHIP = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)


## ================================================================== ##
## Output
## ================================================================== ##

OUT_CATALOG = (
    HERE
    / "annotated_gene_catalog.tsv"
)


OUT_QC = (
    HERE
    / "annotated_gene_catalog_qc.tsv"
)


## ================================================================== ##
## Expected upstream values
## ================================================================== ##

EXPECTED_GENOMES = 631
EXPECTED_GENES = 2002656

EXPECTED_CLUSTERED_PROTEINS = 14308
EXPECTED_CLUSTERS = 1390

EXPECTED_MODULE_PROTEINS = 10537
EXPECTED_MODULE_CLUSTERS = 156
EXPECTED_MODULES = 35


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr
    )

    sys.exit(1)


def require_columns(
    df,
    required,
    filename,
):

    missing = (
        set(required)
        -
        set(df.columns)
    )

    if missing:

        fail(
            f"{filename} missing required "
            f"column(s):\n"
            +
            ", ".join(
                sorted(missing)
            )
        )


## ================================================================== ##
## 1. Read authoritative Prodigal coordinates
## ================================================================== ##

print("=" * 80)
print("BUILD ANNOTATED GENE CATALOG")
print("=" * 80)


print()
print("Reading Prodigal coordinate table...")


coords = pd.read_csv(
    PRODIGAL_COORDS,
    sep="\t",
    dtype={
        "genome": str,
        "protein_id": str,
        "contig": str,
    }
)


require_columns(
    coords,
    [
        "genome",
        "protein_id",
        "contig",
        "start",
        "end",
        "strand",
    ],
    PRODIGAL_COORDS.name
)


coords["start"] = pd.to_numeric(
    coords["start"],
    errors="raise"
).astype(int)


coords["end"] = pd.to_numeric(
    coords["end"],
    errors="raise"
).astype(int)


if len(coords) != EXPECTED_GENES:

    fail(
        f"Expected {EXPECTED_GENES:,} "
        f"Prodigal genes but found "
        f"{len(coords):,}."
    )


if (
    coords[
        [
            "genome",
            "protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome + protein_id "
        "keys in Prodigal coordinate table."
    )


n_genomes = (
    coords["genome"]
    .nunique()
)


if n_genomes != EXPECTED_GENOMES:

    fail(
        f"Expected {EXPECTED_GENOMES} "
        f"genomes but found "
        f"{n_genomes}."
    )


print(
    f"  Genes:   {len(coords):,}"
)

print(
    f"  Genomes: {n_genomes:,}"
)


## ================================================================== ##
## 2. Establish genomic gene order
##
## contig_gene_rank is LEFT-TO-RIGHT genomic order, independent
## of transcriptional strand.
##
## This is intentional:
##
##     rank 1, rank 2, rank 3...
##
## describes physical order along the contig.
## ================================================================== ##

print()
print("Calculating gene order within contigs...")


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
        kind="stable"
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
        sort=False
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
        sort=False
    )[
        "protein_id"
    ]
    .transform(
        "size"
    )
)


## ================================================================== ##
## 3. Read GlobDB annotation mapping
## ================================================================== ##

print()
print("Reading GlobDB annotation mapping...")


annotation = pd.read_csv(
    ANNOTATION_MAP,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


require_columns(
    annotation,
    [
        "genome",
        "protein_id",
        "annotation_match_type",
        "annotation_accepted",
        "reciprocal_overlap",
        "gff_id",
        "globdb_cog",
        "globdb_db_xref",
        "globdb_gene",
        "globdb_product",
    ],
    ANNOTATION_MAP.name
)


if len(annotation) != EXPECTED_GENES:

    fail(
        f"Expected {EXPECTED_GENES:,} "
        f"annotation-map rows but found "
        f"{len(annotation):,}."
    )


if (
    annotation[
        [
            "genome",
            "protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome + protein_id "
        "keys in annotation map."
    )


## ================================================================== ##
## 4. IMPORTANT:
##
## Only ACCEPTED GlobDB mappings are allowed to populate the
## authoritative annotation fields.
##
## Weak-overlap candidate annotations remain available in the
## original mapping-QC table, but they are NOT propagated here.
## ================================================================== ##

accepted_mask = (
    annotation[
        "annotation_accepted"
    ]
    ==
    "1"
)


annotation[
    "globdb_gff_id"
] = ""


annotation.loc[
    accepted_mask,
    "globdb_gff_id"
] = (
    annotation.loc[
        accepted_mask,
        "gff_id"
    ]
)


for column in [
    "globdb_cog",
    "globdb_db_xref",
    "globdb_gene",
    "globdb_product",
]:

    annotation.loc[
        ~accepted_mask,
        column
    ] = ""


annotation_keep = annotation[
    [
        "genome",
        "protein_id",

        "annotation_match_type",
        "annotation_accepted",
        "reciprocal_overlap",

        "globdb_gff_id",
        "globdb_cog",
        "globdb_db_xref",
        "globdb_gene",
        "globdb_product",
    ]
].copy()


## ================================================================== ##
## 5. Join GlobDB annotation to every Prodigal gene
## ================================================================== ##

catalog = coords.merge(
    annotation_keep,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one"
)


if len(catalog) != EXPECTED_GENES:

    fail(
        "GlobDB annotation join changed "
        "the number of genes."
    )


if (
    catalog[
        "annotation_match_type"
    ]
    .isna()
    .any()
):

    fail(
        "One or more Prodigal genes are "
        "missing from the annotation map."
    )


## ================================================================== ##
## 6. Read exported-heme MMseqs cluster membership
##
## Only 14,308 / ~2 million genes belong to this clustered
## protein universe.
## ================================================================== ##

print()
print("Reading MMseqs cluster membership...")


clusters = pd.read_csv(
    CLUSTER_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


require_columns(
    clusters,
    [
        "genome",
        "protein_id",
        "cluster",
    ],
    CLUSTER_MEMBERSHIP.name
)


if len(clusters) != EXPECTED_CLUSTERED_PROTEINS:

    fail(
        f"Expected "
        f"{EXPECTED_CLUSTERED_PROTEINS:,} "
        f"clustered proteins but found "
        f"{len(clusters):,}."
    )


if (
    clusters[
        [
            "genome",
            "protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome + protein_id "
        "keys in cluster membership."
    )


if (
    clusters[
        "cluster"
    ]
    .nunique()
    !=
    EXPECTED_CLUSTERS
):

    fail(
        f"Expected "
        f"{EXPECTED_CLUSTERS:,} "
        f"MMseqs clusters."
    )


## ================================================================== ##
## 7. Join MCL module assignment onto cluster table
##
## Only 156 of the 1,390 MMseqs clusters belong to the final
## connected MCL network.
## ================================================================== ##

modules = pd.read_csv(
    MODULE_MEMBERSHIP,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


require_columns(
    modules,
    [
        "cluster",
        "module",
    ],
    MODULE_MEMBERSHIP.name
)


if (
    modules[
        "cluster"
    ]
    .nunique()
    !=
    EXPECTED_MODULE_CLUSTERS
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_CLUSTERS} "
        f"MCL-assigned clusters."
    )


if (
    modules[
        "module"
    ]
    .nunique()
    !=
    EXPECTED_MODULES
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULES} "
        f"MCL modules."
    )


clusters = clusters.merge(
    modules,
    on="cluster",
    how="left",
    validate="many_to_one"
)


## Empty string rather than NaN for non-module clusters ##

clusters[
    "module"
] = (
    clusters[
        "module"
    ]
    .fillna("")
)


n_module_proteins = int(
    (
        clusters[
            "module"
        ]
        != ""
    )
    .sum()
)


if (
    n_module_proteins
    !=
    EXPECTED_MODULE_PROTEINS
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"module proteins but found "
        f"{n_module_proteins:,}."
    )


## ================================================================== ##
## 8. Select useful protein-level analysis columns
##
## We preserve the rich protein annotation fields from the
## clustering stage when they exist.
## ================================================================== ##

cluster_columns = [
    "genome",
    "protein_id",
    "cluster",
    "module",
]


optional_cluster_columns = [
    "representative_protein",

    "length",
    "clustering_sequence_length",

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

    "export_evidence",
    "localization_class",

    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",

    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate",
]


cluster_columns += [
    column
    for column
    in optional_cluster_columns
    if column
    in clusters.columns
]


clusters_keep = clusters[
    cluster_columns
].copy()


## ================================================================== ##
## 9. Join our protein-family / module evidence onto all genes
## ================================================================== ##

catalog = catalog.merge(
    clusters_keep,
    on=[
        "genome",
        "protein_id",
    ],
    how="left",
    validate="one_to_one"
)


## Structural analysis fields should be explicit empty strings
## outside the exported-heme clustering universe. ##

for column in [
    "cluster",
    "module",
]:

    catalog[
        column
    ] = (
        catalog[
            column
        ]
        .fillna("")
    )


catalog[
    "is_exported_heme_cluster_protein"
] = (
    catalog[
        "cluster"
    ]
    != ""
).astype(int)


catalog[
    "is_mcl_module_protein"
] = (
    catalog[
        "module"
    ]
    != ""
).astype(int)


## ================================================================== ##
## 10. Reorder the authoritative structural columns
## ================================================================== ##

first_columns = [
    "genome",
    "contig",

    "contig_gene_rank",
    "contig_gene_count",

    "protein_id",

    "start",
    "end",
    "gene_length_nt",
    "strand",

    "prodigal_id",
    "prodigal_seqnum",
    "prodigal_gene_ordinal",

    "partial",
    "start_type",
    "rbs_motif",
    "rbs_spacer",
    "gc_cont",

    "annotation_match_type",
    "annotation_accepted",
    "reciprocal_overlap",

    "globdb_gff_id",
    "globdb_cog",
    "globdb_gene",
    "globdb_product",
    "globdb_db_xref",

    "is_exported_heme_cluster_protein",
    "cluster",

    "is_mcl_module_protein",
    "module",
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
    +
    remaining_columns
]


## ================================================================== ##
## 11. Final deterministic genomic ordering
## ================================================================== ##

catalog = (
    catalog
    .sort_values(
        [
            "genome",
            "contig",
            "contig_gene_rank",
        ],
        kind="stable"
    )
    .reset_index(
        drop=True
    )
)


## ================================================================== ##
## 12. Strict QC
## ================================================================== ##

if len(catalog) != EXPECTED_GENES:

    fail(
        "Final catalog does not contain "
        "the expected number of genes."
    )


if (
    catalog[
        [
            "genome",
            "protein_id",
        ]
    ]
    .duplicated()
    .any()
):

    fail(
        "Duplicate genome + protein_id "
        "keys in final gene catalog."
    )


## Check contig ranks are continuous 1..N ##

rank_check = (
    catalog
    .groupby(
        [
            "genome",
            "contig",
        ],
        sort=False
    )
    .agg(
        n_genes=(
            "protein_id",
            "size"
        ),

        min_rank=(
            "contig_gene_rank",
            "min"
        ),

        max_rank=(
            "contig_gene_rank",
            "max"
        ),
    )
)


bad_rank = rank_check[
    (rank_check[
        "min_rank"
    ] != 1)
    |
    (
        rank_check[
            "max_rank"
        ]
        !=
        rank_check[
            "n_genes"
        ]
    )
]


if len(
    bad_rank
) > 0:

    fail(
        "Non-continuous contig gene ranks "
        "detected."
    )


n_clustered_final = int(
    catalog[
        "is_exported_heme_cluster_protein"
    ]
    .sum()
)


n_module_final = int(
    catalog[
        "is_mcl_module_protein"
    ]
    .sum()
)


if (
    n_clustered_final
    !=
    EXPECTED_CLUSTERED_PROTEINS
):

    fail(
        "Final catalog does not recover "
        "all clustered proteins."
    )


if (
    n_module_final
    !=
    EXPECTED_MODULE_PROTEINS
):

    fail(
        "Final catalog does not recover "
        "all MCL module proteins."
    )


## ================================================================== ##
## 13. Write catalogue
## ================================================================== ##

print()
print("Writing annotated gene catalog...")


catalog.to_csv(
    OUT_CATALOG,
    sep="\t",
    index=False,
    na_rep=""
)


## ================================================================== ##
## 14. QC summary
## ================================================================== ##

n_contigs = (
    catalog[
        [
            "genome",
            "contig",
        ]
    ]
    .drop_duplicates()
    .shape[0]
)


n_accepted_globdb = int(
    (
        catalog[
            "annotation_accepted"
        ]
        ==
        "1"
    )
    .sum()
)


n_with_product = int(
    (
        catalog[
            "globdb_product"
        ]
        .fillna("")
        .astype(str)
        .str.strip()
        != ""
    )
    .sum()
)


n_with_cog = int(
    (
        catalog[
            "globdb_cog"
        ]
        .fillna("")
        .astype(str)
        .str.strip()
        != ""
    )
    .sum()
)


qc = pd.DataFrame(
    [
        [
            "genomes",
            n_genomes,
        ],

        [
            "contigs_with_prodigal_genes",
            n_contigs,
        ],

        [
            "total_prodigal_genes",
            len(
                catalog
            ),
        ],

        [
            "accepted_globdb_coordinate_mappings",
            n_accepted_globdb,
        ],

        [
            "genes_with_globdb_product",
            n_with_product,
        ],

        [
            "genes_with_cog_assignment",
            n_with_cog,
        ],

        [
            "exported_heme_cluster_proteins",
            n_clustered_final,
        ],

        [
            "mmseqs_clusters",
            catalog.loc[
                catalog[
                    "cluster"
                ]
                != "",
                "cluster"
            ]
            .nunique(),
        ],

        [
            "mcl_module_proteins",
            n_module_final,
        ],

        [
            "mcl_module_clusters",
            catalog.loc[
                catalog[
                    "module"
                ]
                != "",
                "cluster"
            ]
            .nunique(),
        ],

        [
            "mcl_modules",
            catalog.loc[
                catalog[
                    "module"
                ]
                != "",
                "module"
            ]
            .nunique(),
        ],
    ],
    columns=[
        "metric",
        "value",
    ]
)


qc.to_csv(
    OUT_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## Final report
## ================================================================== ##

print()
print("Annotated gene catalogue summary")

print(
    f"  Genomes:                    "
    f"{n_genomes:,}"
)

print(
    f"  Contigs:                    "
    f"{n_contigs:,}"
)

print(
    f"  Prodigal genes:             "
    f"{len(catalog):,}"
)

print(
    f"  Accepted GlobDB mappings:   "
    f"{n_accepted_globdb:,}"
)

print(
    f"  With GlobDB product:        "
    f"{n_with_product:,}"
)

print(
    f"  With COG:                   "
    f"{n_with_cog:,}"
)

print(
    f"  Exported-heme proteins:     "
    f"{n_clustered_final:,}"
)

print(
    f"  MCL module proteins:        "
    f"{n_module_final:,}"
)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(
    f"Gene catalog: "
    f"{OUT_CATALOG}"
)

print(
    f"QC:           "
    f"{OUT_QC}"
)
