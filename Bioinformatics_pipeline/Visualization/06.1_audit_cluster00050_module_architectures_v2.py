#!/usr/bin/env python3


## ================================================================== ##
## STAGE 16A - CLUSTER_00050 Cyc2 ARCHITECTURE AUDIT v2
##
## Biological-focal deduplication
## ------------------------------
##
## The first audit showed:
##
##   72 Stage-15A focal-region records
##   53 unique Cluster_00050 focal proteins
##   53 genomes
##
## Therefore Stage 15A contains multiple region representations for
## some of the SAME biological Cyc2 focal genes.
##
##
## Authoritative biological focal unit
## -----------------------------------
##
##       genome + focal_protein_id
##
## NOT:
##
##       focal_region_id
##
##
## This script:
##
##   1. resolves the MCL module containing Cluster_00050;
##
##   2. identifies every unique biological Cluster_00050 focal locus;
##
##   3. audits duplicate Stage-15A focal-region representations;
##
##   4. checks that duplicated representations do not disagree about
##      the same neighboring protein;
##
##   5. collapses duplicate representations by taking the UNION of
##      unique neighboring proteins around each biological focal locus;
##
##   6. recomputes the five-cluster community architecture using one
##      row per biological focal locus;
##
##   7. generates corrected priority tables for representative
##      selection.
##
##
## IMPORTANT
## ---------
##
## The union operation does NOT merge different biological loci.
##
## It only combines Stage-15A representations having exactly the same:
##
##       genome + focal_protein_id
##
## Shared neighboring proteins must have consistent coordinates,
## cluster assignments and strand. Any disagreement causes the script
## to stop rather than silently reconcile conflicting data.
##
##
## Primary architecture window:
##
##       +/-10 kb from Cluster_00050 Cyc2
##
## Secondary QC:
##
##       +/-20 kb
##
## ================================================================== ##


from pathlib import Path
import hashlib

import numpy as np
import pandas as pd


## ================================================================== ##
## Constants
## ================================================================== ##

WORKFLOW = Path(
    "~/methanotrophs/methanotroph_project/jeppe/"
    "genome_analysis_workflow"
).expanduser()


FOCAL_CLUSTER = "Cluster_00050"

PRIMARY_WINDOW_BP = 10_000

SECONDARY_WINDOW_BP = 20_000


## ================================================================== ##
## Inputs
## ================================================================== ##

MEMBERSHIP_FILE = (
    WORKFLOW
    / "07_cooccurrence_network"
    / "module_membership.tsv"
)


NEIGHBORHOOD_FILE = (
    WORKFLOW
    / "15_gene_level_analysis"
    / "15A_resolution_recovery"
    / "focal_gene_neighborhoods.tsv"
)


TAXONOMY_FILE = (
    WORKFLOW
    / "16_visualization"
    / "globdb_r226_taxonomy.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUTPUT_DIR = (
    WORKFLOW
    / "16_visualization"
    / "16A_cluster00050_cyc2_gene_map_framework"
)


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


OUT_DUPLICATE_QC = (
    OUTPUT_DIR
    / "Cluster_00050_duplicate_focal_qc.tsv"
)


OUT_CONFLICTS = (
    OUTPUT_DIR
    / "Cluster_00050_neighbor_conflicts.tsv"
)


OUT_UNIQUE_NEIGHBORHOODS = (
    OUTPUT_DIR
    / "Cluster_00050_unique_focal_neighborhoods.tsv"
)


OUT_REGION_CLUSTER_LONG = (
    OUTPUT_DIR
    / "Cluster_00050_region_cluster_presence_long_deduplicated.tsv"
)


OUT_REGION_ARCHITECTURE = (
    OUTPUT_DIR
    / "Cluster_00050_region_architecture_deduplicated.tsv"
)


OUT_PATTERN_SUMMARY = (
    OUTPUT_DIR
    / "Cluster_00050_architecture_patterns_deduplicated.tsv"
)


OUT_PRIORITY = (
    OUTPUT_DIR
    / "Cluster_00050_priority_candidates_deduplicated.tsv"
)


OUT_REPORT = (
    OUTPUT_DIR
    / "Cluster_00050_architecture_audit_v2.txt"
)


## ================================================================== ##
## Helpers
## ================================================================== ##

REPORT_LINES = []


def fail(message):

    raise RuntimeError(
        message
    )


def report(message=""):

    message = str(
        message
    )

    print(
        message
    )

    REPORT_LINES.append(
        message
    )


def heading(message):

    report("")

    report(
        "=" * 110
    )

    report(
        message
    )

    report(
        "=" * 110
    )


def require_file(path):

    if not path.exists():

        fail(
            f"Required file not found:\n{path}"
        )


def clean_string(series):

    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
    )


def safe_numeric(series):

    return pd.to_numeric(
        series,
        errors="coerce"
    )


def parse_taxonomy_value(
    taxonomy,
    prefix
):

    taxonomy = str(
        taxonomy
    )

    for token in taxonomy.split(";"):

        token = token.strip()

        if token.startswith(
            prefix
        ):

            value = token[
                len(prefix):
            ].strip()

            if value:

                return value

    return ""


def join_unique(values):

    values = sorted(
        {
            str(x).strip()
            for x in values
            if pd.notna(x)
            and str(x).strip()
        }
    )

    return ";".join(
        values
    )


def hash_signature(rows):

    if len(rows) == 0:

        return "EMPTY"


    values = []


    for _, row in rows.iterrows():

        values.append(
            "|".join(
                [
                    str(
                        row[
                            "neighbor_protein_id"
                        ]
                    ),
                    str(
                        row[
                            "neighbor_mmseq_cluster"
                        ]
                    ),
                    str(
                        row[
                            "oriented_midpoint_offset_bp"
                        ]
                    ),
                    str(
                        row[
                            "oriented_gene_offset"
                        ]
                    ),
                    str(
                        row[
                            "neighbor_strand"
                        ]
                    ),
                ]
            )
        )


    text = "\n".join(
        sorted(
            values
        )
    )


    return hashlib.sha256(
        text.encode(
            "utf-8"
        )
    ).hexdigest()


## ================================================================== ##
## Validate inputs
## ================================================================== ##

for path in [

    MEMBERSHIP_FILE,
    NEIGHBORHOOD_FILE,
    TAXONOMY_FILE,

]:

    require_file(
        path
    )


## ================================================================== ##
## 1. Resolve Cluster_00050 module
## ================================================================== ##

heading(
    "1. RESOLVE MCL COMMUNITY"
)


membership = pd.read_csv(
    MEMBERSHIP_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


required_membership = {
    "cluster",
    "module",
}


missing = (
    required_membership
    -
    set(
        membership.columns
    )
)


if missing:

    fail(
        "Missing module-membership columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


membership[
    "cluster"
] = clean_string(
    membership[
        "cluster"
    ]
)


membership[
    "module"
] = clean_string(
    membership[
        "module"
    ]
)


focal_modules = sorted(

    membership.loc[
        membership[
            "cluster"
        ]
        ==
        FOCAL_CLUSTER,
        "module"
    ]

    .unique()

)


if len(
    focal_modules
) != 1:

    fail(
        f"{FOCAL_CLUSTER} maps to {len(focal_modules)} modules: "
        f"{focal_modules}"
    )


MCL_MODULE = focal_modules[0]


module_clusters = sorted(

    membership.loc[
        membership[
            "module"
        ]
        ==
        MCL_MODULE,
        "cluster"
    ]

    .drop_duplicates()

    .tolist()

)


N_MODULE_CLUSTERS = len(
    module_clusters
)


report(
    f"Focal cluster: {FOCAL_CLUSTER}"
)

report(
    f"MCL module:    {MCL_MODULE}"
)

report(
    f"Members:       {N_MODULE_CLUSTERS}"
)

report("")


for cluster in module_clusters:

    suffix = (
        "  <-- Cyc2"
        if cluster == FOCAL_CLUSTER
        else ""
    )

    report(
        f"  {cluster}{suffix}"
    )


if N_MODULE_CLUSTERS != 5:

    report("")

    report(
        "WARNING: module no longer contains exactly five clusters. "
        "All calculations will use the observed membership."
    )


## ================================================================== ##
## 2. Read Stage-15A neighborhoods
## ================================================================== ##

heading(
    "2. READ STAGE-15A NEIGHBORHOODS"
)


neigh = pd.read_csv(
    NEIGHBORHOOD_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)


required_columns = {

    "genome",

    "focal_region_id",
    "focal_protein_id",
    "focal_mmseq_cluster",
    "focal_contig",

    "neighbor_protein_id",
    "neighbor_mmseq_cluster",
    "neighbor_contig",
    "neighbor_strand",

    "oriented_midpoint_offset_bp",
    "oriented_gene_offset",

}


missing = (
    required_columns
    -
    set(
        neigh.columns
    )
)


if missing:

    fail(
        "Missing Stage-15A neighborhood columns: "
        +
        ", ".join(
            sorted(
                missing
            )
        )
    )


report(
    f"Total Stage-15A rows: {len(neigh):,}"
)


## ================================================================== ##
## 3. Restrict to Cluster_00050 focal loci
## ================================================================== ##

heading(
    "3. IDENTIFY BIOLOGICAL Cluster_00050 FOCALS"
)


cyc2 = neigh.loc[
    clean_string(
        neigh[
            "focal_mmseq_cluster"
        ]
    )
    ==
    FOCAL_CLUSTER
].copy()


## ------------------------------------------------------------------ ##
## Same contig only
## ------------------------------------------------------------------ ##

cyc2 = cyc2.loc[
    clean_string(
        cyc2[
            "neighbor_contig"
        ]
    )
    ==
    clean_string(
        cyc2[
            "focal_contig"
        ]
    )
].copy()


cyc2[
    "genome"
] = clean_string(
    cyc2[
        "genome"
    ]
)


cyc2[
    "focal_protein_id"
] = clean_string(
    cyc2[
        "focal_protein_id"
    ]
)


cyc2[
    "focal_region_id"
] = clean_string(
    cyc2[
        "focal_region_id"
    ]
)


cyc2[
    "neighbor_protein_id"
] = clean_string(
    cyc2[
        "neighbor_protein_id"
    ]
)


cyc2[
    "neighbor_mmseq_cluster"
] = clean_string(
    cyc2[
        "neighbor_mmseq_cluster"
    ]
)


cyc2[
    "neighbor_strand"
] = clean_string(
    cyc2[
        "neighbor_strand"
    ]
)


cyc2[
    "oriented_midpoint_offset_bp"
] = safe_numeric(
    cyc2[
        "oriented_midpoint_offset_bp"
    ]
)


cyc2[
    "oriented_gene_offset"
] = safe_numeric(
    cyc2[
        "oriented_gene_offset"
    ]
)


cyc2[
    "biological_focal_id"
] = (

    cyc2[
        "genome"
    ]

    +
    "||"

    +
    cyc2[
        "focal_protein_id"
    ]

)


n_stage15_regions = cyc2[
    "focal_region_id"
].nunique()


n_biological_focals = cyc2[
    "biological_focal_id"
].nunique()


n_focal_genomes = cyc2[
    "genome"
].nunique()


report(
    f"Stage-15A focal-region IDs: {n_stage15_regions}"
)

report(
    f"Unique biological focals:   {n_biological_focals}"
)

report(
    f"Unique genomes:             {n_focal_genomes}"
)


## ================================================================== ##
## 4. Check focal identity consistency
## ================================================================== ##

heading(
    "4. FOCAL-ID CONSISTENCY QC"
)


focal_identity_qc = (

    cyc2

    .groupby(
        "biological_focal_id"
    )

    .agg(

        genome=(
            "genome",
            "first"
        ),

        focal_protein_id=(
            "focal_protein_id",
            "first"
        ),

        n_focal_contigs=(
            "focal_contig",
            "nunique"
        ),

        focal_contigs=(
            "focal_contig",
            join_unique
        ),

        n_stage15_regions=(
            "focal_region_id",
            "nunique"
        ),

        stage15_region_ids=(
            "focal_region_id",
            join_unique
        ),

    )

    .reset_index()

)


bad_contigs = focal_identity_qc.loc[
    focal_identity_qc[
        "n_focal_contigs"
    ]
    !=
    1
]


if not bad_contigs.empty:

    print(
        bad_contigs.to_string(
            index=False
        )
    )

    fail(
        "At least one genome + focal_protein_id maps to multiple "
        "focal contigs. Refusing to deduplicate."
    )


duplicate_focals = focal_identity_qc.loc[
    focal_identity_qc[
        "n_stage15_regions"
    ]
    >
    1
].copy()


report(
    f"Biological focals represented by >1 Stage-15A region: "
    f"{len(duplicate_focals)}"
)


## ================================================================== ##
## 5. Compare duplicate Stage-15A neighborhood signatures
## ================================================================== ##

heading(
    "5. DUPLICATE REGION SIGNATURE QC"
)


signature_rows = []


for (
    biological_focal_id,
    focal_region_id
), x in cyc2.groupby(

    [
        "biological_focal_id",
        "focal_region_id",
    ],

    sort=False

):


    x = x.copy()


    x10 = x.loc[
        x[
            "oriented_midpoint_offset_bp"
        ].abs()
        <=
        PRIMARY_WINDOW_BP
    ]


    x20 = x.loc[
        x[
            "oriented_midpoint_offset_bp"
        ].abs()
        <=
        SECONDARY_WINDOW_BP
    ]


    signature_rows.append(
        {

            "biological_focal_id":
                biological_focal_id,

            "focal_region_id":
                focal_region_id,

            "n_neighbors_total":
                x[
                    "neighbor_protein_id"
                ].nunique(),

            "n_neighbors_10kb":
                x10[
                    "neighbor_protein_id"
                ].nunique(),

            "n_neighbors_20kb":
                x20[
                    "neighbor_protein_id"
                ].nunique(),

            "signature_total":
                hash_signature(
                    x
                ),

            "signature_10kb":
                hash_signature(
                    x10
                ),

            "signature_20kb":
                hash_signature(
                    x20
                ),

        }
    )


region_signatures = pd.DataFrame(
    signature_rows
)


signature_qc = (

    region_signatures

    .groupby(
        "biological_focal_id"
    )

    .agg(

        n_stage15_regions=(
            "focal_region_id",
            "nunique"
        ),

        n_distinct_total_signatures=(
            "signature_total",
            "nunique"
        ),

        n_distinct_10kb_signatures=(
            "signature_10kb",
            "nunique"
        ),

        n_distinct_20kb_signatures=(
            "signature_20kb",
            "nunique"
        ),

    )

    .reset_index()

)


duplicate_qc = focal_identity_qc.merge(

    signature_qc,

    on=[
        "biological_focal_id",
        "n_stage15_regions",
    ],

    how="left",

    validate="one_to_one"

)


duplicate_qc[
    "duplicate_regions_identical_10kb"
] = np.where(

    duplicate_qc[
        "n_stage15_regions"
    ]
    <=
    1,

    True,

    duplicate_qc[
        "n_distinct_10kb_signatures"
    ]
    ==
    1

)


duplicate_qc.to_csv(
    OUT_DUPLICATE_QC,
    sep="\t",
    index=False
)


report(
    f"Duplicate QC written:\n{OUT_DUPLICATE_QC}"
)


if len(
    duplicate_focals
) > 0:

    report("")

    report(
        "Duplicate biological focals:"
    )

    print(

        duplicate_qc.loc[
            duplicate_qc[
                "n_stage15_regions"
            ]
            >
            1,

            [
                "genome",
                "focal_protein_id",

                "n_stage15_regions",

                "n_distinct_10kb_signatures",
                "n_distinct_20kb_signatures",
                "n_distinct_total_signatures",

                "duplicate_regions_identical_10kb",
            ]

        ].to_string(
            index=False
        )

    )


## ================================================================== ##
## 6. Check shared-neighbor consistency
##
## If the same biological focal + neighboring protein occurs through
## multiple Stage-15A regions, all biological coordinates/assignments
## must agree.
## ================================================================== ##

heading(
    "6. SHARED NEIGHBOR CONSISTENCY QC"
)


conflict_rows = []


for (
    biological_focal_id,
    neighbor_protein_id
), x in cyc2.groupby(

    [
        "biological_focal_id",
        "neighbor_protein_id",
    ],

    sort=False

):


    tests = {

        "neighbor_contig":
            x[
                "neighbor_contig"
            ].nunique(),

        "neighbor_mmseq_cluster":
            x[
                "neighbor_mmseq_cluster"
            ].nunique(),

        "neighbor_strand":
            x[
                "neighbor_strand"
            ].nunique(),

        "oriented_midpoint_offset_bp":
            x[
                "oriented_midpoint_offset_bp"
            ].dropna().nunique(),

        "oriented_gene_offset":
            x[
                "oriented_gene_offset"
            ].dropna().nunique(),

    }


    conflicting_fields = [

        field

        for field, n_unique in tests.items()

        if n_unique > 1

    ]


    if conflicting_fields:

        conflict_rows.append(
            {

                "biological_focal_id":
                    biological_focal_id,

                "genome":
                    x[
                        "genome"
                    ].iloc[0],

                "focal_protein_id":
                    x[
                        "focal_protein_id"
                    ].iloc[0],

                "neighbor_protein_id":
                    neighbor_protein_id,

                "n_source_regions":
                    x[
                        "focal_region_id"
                    ].nunique(),

                "source_regions":
                    join_unique(
                        x[
                            "focal_region_id"
                        ]
                    ),

                "conflicting_fields":
                    ";".join(
                        conflicting_fields
                    ),

            }
        )


conflicts = pd.DataFrame(
    conflict_rows
)


if conflicts.empty:

    pd.DataFrame(
        columns=[
            "biological_focal_id",
            "genome",
            "focal_protein_id",
            "neighbor_protein_id",
            "n_source_regions",
            "source_regions",
            "conflicting_fields",
        ]
    ).to_csv(
        OUT_CONFLICTS,
        sep="\t",
        index=False
    )

    report(
        "Shared-neighbor conflicts: 0"
    )


else:

    conflicts.to_csv(
        OUT_CONFLICTS,
        sep="\t",
        index=False
    )

    print(
        conflicts.to_string(
            index=False
        )
    )

    fail(
        "\nConflicting biological information was found for shared "
        "neighbor proteins.\n"
        f"Inspect:\n{OUT_CONFLICTS}"
    )


## ================================================================== ##
## 7. Collapse to one biological focal neighborhood
##
## UNION unique neighbor proteins across duplicate Stage-15A regions.
##
## Because shared-neighbor consistency has already been validated,
## dropping duplicate neighbor rows is safe.
## ================================================================== ##

heading(
    "7. BUILD UNIQUE BIOLOGICAL FOCAL NEIGHBORHOODS"
)


source_region_info = (

    cyc2

    .groupby(
        "biological_focal_id"
    )

    .agg(

        n_source_focal_regions=(
            "focal_region_id",
            "nunique"
        ),

        source_focal_region_ids=(
            "focal_region_id",
            join_unique
        ),

    )

    .reset_index()

)


unique_neigh = (

    cyc2

    .sort_values(
        [
            "biological_focal_id",
            "neighbor_protein_id",
            "focal_region_id",
        ]
    )

    .drop_duplicates(
        [
            "biological_focal_id",
            "neighbor_protein_id",
        ]
    )

    .copy()

)


unique_neigh = unique_neigh.drop(
    columns=[
        "focal_region_id",
    ]
)


unique_neigh = unique_neigh.rename(
    columns={
        "biological_focal_id":
            "plot_region_id"
    }
)


unique_neigh = unique_neigh.merge(

    source_region_info.rename(
        columns={
            "biological_focal_id":
                "plot_region_id"
        }
    ),

    on="plot_region_id",

    how="left",

    validate="many_to_one"

)


unique_neigh.to_csv(
    OUT_UNIQUE_NEIGHBORHOODS,
    sep="\t",
    index=False
)


report(
    f"Unique biological focal loci: "
    f"{unique_neigh['plot_region_id'].nunique()}"
)

report(
    f"Unique focal-neighborhood rows: {len(unique_neigh):,}"
)

report(
    f"Clean neighborhood file:\n{OUT_UNIQUE_NEIGHBORHOODS}"
)


if (
    unique_neigh[
        "plot_region_id"
    ].nunique()
    !=
    n_biological_focals
):

    fail(
        "Unexpected biological focal count after deduplication."
    )


## ================================================================== ##
## 8. Read taxonomy
## ================================================================== ##

heading(
    "8. ADD GLOBDB r226 TAXONOMY"
)


taxonomy = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t",
    header=None,
    names=[
        "genome",
        "taxonomy",
    ],
    usecols=[
        0,
        1,
    ],
    dtype=str,
    keep_default_na=False
)


taxonomy = taxonomy.drop_duplicates(
    "genome"
)


taxonomy[
    "taxonomy_family"
] = taxonomy[
    "taxonomy"
].map(
    lambda x:
    parse_taxonomy_value(
        x,
        "f__"
    )
)


taxonomy[
    "taxonomy_genus"
] = taxonomy[
    "taxonomy"
].map(
    lambda x:
    parse_taxonomy_value(
        x,
        "g__"
    )
)


taxonomy[
    "taxonomy_species"
] = taxonomy[
    "taxonomy"
].map(
    lambda x:
    parse_taxonomy_value(
        x,
        "s__"
    )
)


taxonomy[
    "taxonomy_display"
] = np.where(

    taxonomy[
        "taxonomy_species"
    ]
    !=
    "",

    taxonomy[
        "taxonomy_species"
    ],

    np.where(

        taxonomy[
            "taxonomy_genus"
        ]
        !=
        "",

        taxonomy[
            "taxonomy_genus"
        ],

        taxonomy[
            "genome"
        ]

    )

)


taxonomy[
    "priority_taxon"
] = np.select(

    [

        taxonomy[
            "taxonomy_genus"
        ].str.startswith(
            "Methylobacter",
            na=False
        ),

        taxonomy[
            "taxonomy_genus"
        ].str.startswith(
            "Crenothrix",
            na=False
        ),

    ],

    [

        "Methylobacter",
        "Crenothrix",

    ],

    default="Other"

)


## ================================================================== ##
## 9. Restrict to the five Module_10 families
## ================================================================== ##

heading(
    "9. RECOMPUTE LOCAL COMMUNITY ARCHITECTURE"
)


module_neigh = unique_neigh.loc[
    unique_neigh[
        "neighbor_mmseq_cluster"
    ].isin(
        module_clusters
    )
].copy()


module_neigh[
    "abs_midpoint_offset_bp"
] = module_neigh[
    "oriented_midpoint_offset_bp"
].abs()


module_neigh[
    "within_10kb"
] = (
    module_neigh[
        "abs_midpoint_offset_bp"
    ]
    <=
    PRIMARY_WINDOW_BP
)


module_neigh[
    "within_20kb"
] = (
    module_neigh[
        "abs_midpoint_offset_bp"
    ]
    <=
    SECONDARY_WINDOW_BP
)


## ================================================================== ##
## 10. Biological focal metadata
## ================================================================== ##

focal_meta = (

    unique_neigh

    [
        [
            "plot_region_id",
            "genome",
            "focal_protein_id",
            "focal_contig",

            "n_source_focal_regions",
            "source_focal_region_ids",
        ]
    ]

    .drop_duplicates()

)


if (
    focal_meta[
        "plot_region_id"
    ].duplicated().any()
):

    fail(
        "A biological focal locus has inconsistent focal metadata "
        "after deduplication."
    )


## ================================================================== ##
## 11. One row per biological focal x module cluster
## ================================================================== ##

region_cluster_rows = []


for _, focal in focal_meta.iterrows():

    region_id = focal[
        "plot_region_id"
    ]


    region_rows = module_neigh.loc[
        module_neigh[
            "plot_region_id"
        ]
        ==
        region_id
    ]


    for cluster in module_clusters:

        x = region_rows.loc[
            region_rows[
                "neighbor_mmseq_cluster"
            ]
            ==
            cluster
        ].copy()


        is_focal_cluster = (
            cluster
            ==
            FOCAL_CLUSTER
        )


        if not x.empty:

            x = x.sort_values(
                [
                    "abs_midpoint_offset_bp",
                    "oriented_gene_offset",
                    "neighbor_protein_id",
                ],
                na_position="last"
            )


            nearest = x.iloc[0]


            present10 = bool(
                x[
                    "within_10kb"
                ].any()
            )


            present20 = bool(
                x[
                    "within_20kb"
                ].any()
            )


            n10 = x.loc[
                x[
                    "within_10kb"
                ],
                "neighbor_protein_id"
            ].nunique()


            n20 = x.loc[
                x[
                    "within_20kb"
                ],
                "neighbor_protein_id"
            ].nunique()


            nall = x[
                "neighbor_protein_id"
            ].nunique()


            nearest_protein = nearest[
                "neighbor_protein_id"
            ]


            nearest_midpoint = nearest[
                "oriented_midpoint_offset_bp"
            ]


            nearest_gene_offset = nearest[
                "oriented_gene_offset"
            ]


            nearest_strand = nearest[
                "neighbor_strand"
            ]


        else:

            present10 = False
            present20 = False

            n10 = 0
            n20 = 0
            nall = 0

            nearest_protein = ""
            nearest_midpoint = np.nan
            nearest_gene_offset = np.nan
            nearest_strand = ""


        ## ---------------------------------------------------------- ##
        ## Cyc2 itself is present by definition.
        ## ---------------------------------------------------------- ##

        if is_focal_cluster:

            present10 = True
            present20 = True

            n10 = max(
                n10,
                1
            )

            n20 = max(
                n20,
                1
            )

            nall = max(
                nall,
                1
            )


            if nearest_protein == "":

                nearest_protein = focal[
                    "focal_protein_id"
                ]

                nearest_midpoint = 0
                nearest_gene_offset = 0


        region_cluster_rows.append(
            {

                "plot_region_id":
                    region_id,

                "genome":
                    focal[
                        "genome"
                    ],

                "focal_protein_id":
                    focal[
                        "focal_protein_id"
                    ],

                "module_cluster":
                    cluster,

                "is_focal_cyc2_cluster":
                    int(
                        is_focal_cluster
                    ),

                "present_within_10kb":
                    int(
                        present10
                    ),

                "present_within_20kb":
                    int(
                        present20
                    ),

                "n_copies_within_10kb":
                    int(
                        n10
                    ),

                "n_copies_within_20kb":
                    int(
                        n20
                    ),

                "n_copies_extracted":
                    int(
                        nall
                    ),

                "nearest_protein_id":
                    nearest_protein,

                "nearest_midpoint_offset_bp":
                    nearest_midpoint,

                "nearest_gene_offset":
                    nearest_gene_offset,

                "nearest_strand":
                    nearest_strand,

            }
        )


region_cluster_long = pd.DataFrame(
    region_cluster_rows
)


region_cluster_long.to_csv(
    OUT_REGION_CLUSTER_LONG,
    sep="\t",
    index=False
)


expected_rows = (
    n_biological_focals
    *
    N_MODULE_CLUSTERS
)


if len(
    region_cluster_long
) != expected_rows:

    fail(
        f"Expected {expected_rows} region x cluster rows, "
        f"found {len(region_cluster_long)}."
    )


report(
    f"Region x cluster rows: {len(region_cluster_long):,}"
)


## ================================================================== ##
## 12. Region-level architecture table
## ================================================================== ##

heading(
    "10. BUILD DEDUPLICATED ARCHITECTURE MATRIX"
)


architecture_rows = []


for _, focal in focal_meta.iterrows():

    region_id = focal[
        "plot_region_id"
    ]


    x = region_cluster_long.loc[
        region_cluster_long[
            "plot_region_id"
        ]
        ==
        region_id
    ]


    present10 = [

        cluster

        for cluster in module_clusters

        if int(
            x.loc[
                x[
                    "module_cluster"
                ]
                ==
                cluster,
                "present_within_10kb"
            ].iloc[0]
        )
        ==
        1

    ]


    present20 = [

        cluster

        for cluster in module_clusters

        if int(
            x.loc[
                x[
                    "module_cluster"
                ]
                ==
                cluster,
                "present_within_20kb"
            ].iloc[0]
        )
        ==
        1

    ]


    missing10 = [

        cluster

        for cluster in module_clusters

        if cluster not in present10

    ]


    row = {

        "plot_region_id":
            region_id,

        "genome":
            focal[
                "genome"
            ],

        "focal_protein_id":
            focal[
                "focal_protein_id"
            ],

        "focal_cluster":
            FOCAL_CLUSTER,

        "mcl_module":
            MCL_MODULE,

        "n_source_focal_regions":
            focal[
                "n_source_focal_regions"
            ],

        "source_focal_region_ids":
            focal[
                "source_focal_region_ids"
            ],

        "n_module_clusters_total":
            N_MODULE_CLUSTERS,

        "n_module_clusters_10kb":
            len(
                present10
            ),

        "n_module_clusters_20kb":
            len(
                present20
            ),

        "clusters_10kb":
            ";".join(
                present10
            ),

        "clusters_20kb":
            ";".join(
                present20
            ),

        "missing_clusters_10kb":
            ";".join(
                missing10
            ),

        "architecture_signature_10kb":
            "+".join(
                present10
            ),

        "is_complete_10kb":
            int(
                len(
                    present10
                )
                ==
                N_MODULE_CLUSTERS
            ),

    }


    ## -------------------------------------------------------------- ##
    ## Explicit columns for all five module members
    ## -------------------------------------------------------------- ##

    for cluster in module_clusters:

        z = x.loc[
            x[
                "module_cluster"
            ]
            ==
            cluster
        ].iloc[0]


        short = cluster.replace(
            "Cluster_",
            "C"
        )


        row[
            f"{short}_present_10kb"
        ] = int(
            z[
                "present_within_10kb"
            ]
        )


        row[
            f"{short}_copies_10kb"
        ] = int(
            z[
                "n_copies_within_10kb"
            ]
        )


        row[
            f"{short}_nearest_midpoint_bp"
        ] = z[
            "nearest_midpoint_offset_bp"
        ]


        row[
            f"{short}_nearest_gene_offset"
        ] = z[
            "nearest_gene_offset"
        ]


        row[
            f"{short}_nearest_strand"
        ] = z[
            "nearest_strand"
        ]


    architecture_rows.append(
        row
    )


architecture = pd.DataFrame(
    architecture_rows
)


## ================================================================== ##
## 13. Add taxonomy
## ================================================================== ##

architecture = architecture.merge(

    taxonomy[
        [
            "genome",

            "taxonomy_family",
            "taxonomy_genus",
            "taxonomy_species",
            "taxonomy_display",

            "priority_taxon",
        ]
    ],

    on="genome",

    how="left",

    validate="many_to_one"

)


for col in [

    "taxonomy_family",
    "taxonomy_genus",
    "taxonomy_species",
    "taxonomy_display",
    "priority_taxon",

]:

    architecture[
        col
    ] = architecture[
        col
    ].fillna("")


architecture.loc[
    architecture[
        "taxonomy_display"
    ]
    ==
    "",
    "taxonomy_display"
] = architecture[
    "genome"
]


architecture.loc[
    architecture[
        "priority_taxon"
    ]
    ==
    "",
    "priority_taxon"
] = "Other"


## ================================================================== ##
## 14. Architecture patterns
## ================================================================== ##

pattern_summary = (

    architecture

    .groupby(
        "architecture_signature_10kb",
        dropna=False
    )

    .agg(

        n_regions=(
            "plot_region_id",
            "nunique"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

        n_module_clusters=(
            "n_module_clusters_10kb",
            "first"
        ),

        n_methylobacter=(
            "priority_taxon",
            lambda x:
            int(
                (
                    x
                    ==
                    "Methylobacter"
                ).sum()
            )
        ),

        n_crenothrix=(
            "priority_taxon",
            lambda x:
            int(
                (
                    x
                    ==
                    "Crenothrix"
                ).sum()
            )
        ),

        taxa=(
            "taxonomy_display",
            join_unique
        ),

    )

    .reset_index()

    .sort_values(
        [
            "n_module_clusters",
            "n_regions",
        ],
        ascending=[
            False,
            False,
        ]
    )

    .reset_index(
        drop=True
    )

)


pattern_summary[
    "architecture_pattern_id"
] = [

    f"C00050_arch_{i:02d}"

    for i in range(
        1,
        len(
            pattern_summary
        )
        +
        1
    )

]


architecture = architecture.merge(

    pattern_summary[
        [
            "architecture_signature_10kb",
            "architecture_pattern_id",
            "n_regions",
        ]
    ].rename(
        columns={
            "n_regions":
                "architecture_pattern_n_regions"
        }
    ),

    on="architecture_signature_10kb",

    how="left",

    validate="many_to_one"

)


## ================================================================== ##
## 15. Explicit priority hierarchy
##
## P1:
##     all module clusters locally present
##
## P2:
##     Methylobacter / Crenothrix with >=2 clusters
##
## P3+:
##     remaining taxa ordered by decreasing completeness
##
## No opaque composite score.
## ================================================================== ##

def assign_priority(row):

    n = int(
        row[
            "n_module_clusters_10kb"
        ]
    )


    priority_taxon = (
        row[
            "priority_taxon"
        ]
        in
        {
            "Methylobacter",
            "Crenothrix",
        }
    )


    ## -------------------------------------------------------------- ##
    ## Complete architecture
    ## -------------------------------------------------------------- ##

    if n == N_MODULE_CLUSTERS:

        if priority_taxon:

            return pd.Series(
                [
                    1,
                    0,
                    "P1_complete_priority_taxon",
                    (
                        f"Complete {n}/{N_MODULE_CLUSTERS} architecture; "
                        f"{row['priority_taxon']}"
                    ),
                ]
            )


        return pd.Series(
            [
                1,
                1,
                "P1_complete_other_taxon",
                (
                    f"Complete {n}/{N_MODULE_CLUSTERS} architecture"
                ),
            ]
        )


    ## -------------------------------------------------------------- ##
    ## Methylobacter / Crenothrix partial architecture
    ## -------------------------------------------------------------- ##

    if (
        priority_taxon
        and
        n >= 2
    ):

        return pd.Series(
            [
                2,
                -n,
                "P2_priority_taxon_partial",
                (
                    f"{row['priority_taxon']} with "
                    f"{n}/{N_MODULE_CLUSTERS} local module clusters"
                ),
            ]
        )


    ## -------------------------------------------------------------- ##
    ## Other taxa
    ## -------------------------------------------------------------- ##

    if n >= 2:

        missing = (
            N_MODULE_CLUSTERS
            -
            n
        )


        rank = (
            2
            +
            missing
        )


        return pd.Series(
            [
                rank,
                0,
                f"P{rank}_{n}_of_{N_MODULE_CLUSTERS}_other_taxon",
                (
                    f"Other taxon with "
                    f"{n}/{N_MODULE_CLUSTERS} local module clusters"
                ),
            ]
        )


    return pd.Series(
        [
            99,
            0,
            "not_priority_lt2_clusters",
            (
                f"Only {n}/{N_MODULE_CLUSTERS} local module clusters"
            ),
        ]
    )


architecture[
    [
        "priority_rank",
        "priority_subrank",
        "priority_class",
        "priority_reason",
    ]
] = architecture.apply(
    assign_priority,
    axis=1
)


## ================================================================== ##
## 16. Candidate ordering
## ================================================================== ##

architecture = architecture.sort_values(

    [
        "priority_rank",
        "priority_subrank",

        "n_module_clusters_10kb",

        "architecture_pattern_n_regions",

        "taxonomy_family",
        "taxonomy_genus",
        "taxonomy_species",

        "genome",
        "focal_protein_id",
    ],

    ascending=[
        True,
        True,

        False,

        True,

        True,
        True,
        True,

        True,
        True,
    ]

).reset_index(
    drop=True
)


architecture[
    "candidate_order"
] = np.arange(
    1,
    len(
        architecture
    )
    +
    1
)


## ================================================================== ##
## 17. Write final corrected outputs
## ================================================================== ##

architecture.to_csv(
    OUT_REGION_ARCHITECTURE,
    sep="\t",
    index=False
)


pattern_summary.to_csv(
    OUT_PATTERN_SUMMARY,
    sep="\t",
    index=False
)


priority_candidates = architecture.loc[
    architecture[
        "n_module_clusters_10kb"
    ]
    >=
    2
].copy()


priority_candidates.to_csv(
    OUT_PRIORITY,
    sep="\t",
    index=False
)


## ================================================================== ##
## 18. Human-readable summary
## ================================================================== ##

heading(
    "11. DEDUPLICATED LOCAL COMPLETENESS"
)


completeness = (

    architecture

    .groupby(
        "n_module_clusters_10kb"
    )

    .agg(

        n_focal_loci=(
            "plot_region_id",
            "nunique"
        ),

        n_genomes=(
            "genome",
            "nunique"
        ),

    )

    .reset_index()

    .sort_values(
        "n_module_clusters_10kb",
        ascending=False
    )

)


report(
    completeness.to_string(
        index=False
    )
)


## ================================================================== ##
## Priority taxa
## ================================================================== ##

heading(
    "12. METHYLOBACTER / CRENOTHRIX AFTER DEDUPLICATION"
)


priority_taxa = architecture.loc[

    architecture[
        "priority_taxon"
    ].isin(
        [
            "Methylobacter",
            "Crenothrix",
        ]
    ),

    [
        "candidate_order",

        "priority_rank",
        "priority_class",

        "genome",
        "taxonomy_display",
        "priority_taxon",

        "focal_protein_id",

        "n_module_clusters_10kb",

        "architecture_pattern_id",

        "clusters_10kb",
        "missing_clusters_10kb",

    ]

]


if priority_taxa.empty:

    report(
        "No Methylobacter or Crenothrix loci."
    )

else:

    report(
        priority_taxa.to_string(
            index=False
        )
    )


## ================================================================== ##
## Architecture patterns
## ================================================================== ##

heading(
    "13. DEDUPLICATED ARCHITECTURE PATTERNS"
)


report(
    pattern_summary.to_string(
        index=False
    )
)


## ================================================================== ##
## Top candidates
## ================================================================== ##

heading(
    "14. TOP DEDUPLICATED REPRESENTATIVE CANDIDATES"
)


candidate_columns = [

    "candidate_order",

    "priority_rank",
    "priority_class",

    "genome",
    "taxonomy_display",
    "priority_taxon",

    "focal_protein_id",

    "n_module_clusters_10kb",

    "architecture_pattern_id",
    "architecture_pattern_n_regions",

    "clusters_10kb",
    "missing_clusters_10kb",

]


report(

    priority_candidates[
        candidate_columns
    ]

    .head(
        60
    )

    .to_string(
        index=False
    )

)


## ================================================================== ##
## 19. Final QC
## ================================================================== ##

heading(
    "15. FINAL QC"
)


report(
    f"Original Stage-15A region IDs:       {n_stage15_regions}"
)

report(
    f"Unique biological Cyc2 focal loci:  {n_biological_focals}"
)

report(
    f"Final architecture rows:            {len(architecture)}"
)

report(
    f"Final unique genomes:               {architecture['genome'].nunique()}"
)


if len(
    architecture
) != n_biological_focals:

    fail(
        "Final architecture table does not contain exactly one row "
        "per biological focal locus."
    )


if architecture[
    [
        "genome",
        "focal_protein_id",
    ]
].duplicated().any():

    fail(
        "Duplicate genome + focal_protein_id remained in final "
        "architecture table."
    )


report("")

report(
    "PASS: exactly one architecture row per genome + focal protein."
)


## ================================================================== ##
## 20. Output files
## ================================================================== ##

heading(
    "16. OUTPUT FILES"
)


for path in [

    OUT_DUPLICATE_QC,
    OUT_CONFLICTS,
    OUT_UNIQUE_NEIGHBORHOODS,

    OUT_REGION_CLUSTER_LONG,
    OUT_REGION_ARCHITECTURE,
    OUT_PATTERN_SUMMARY,
    OUT_PRIORITY,

]:

    report(
        str(
            path
        )
    )


## ================================================================== ##
## Write report
## ================================================================== ##

OUT_REPORT.write_text(

    "\n".join(
        REPORT_LINES
    )
    +
    "\n",

    encoding="utf-8"

)


print("")

print(
    "Cluster_00050 architecture audit v2 complete."
)

print(
    f"Report:\n{OUT_REPORT}"
)
