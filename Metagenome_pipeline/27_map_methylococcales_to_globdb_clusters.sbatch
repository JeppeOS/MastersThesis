#!/usr/bin/env bash

#SBATCH --job-name=MG_vs_GlobDB_clusters
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/globdb_cluster_mapping/logs/mmseqs_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/metagenome/functional_analysis/globdb_cluster_mapping/logs/mmseqs_%j.err

set -Eeo pipefail


###############################################################################
## Paths
###############################################################################

PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

WORK_DIR="${PROJECT}/functional_analysis/globdb_cluster_mapping"


QUERY_FASTA="${WORK_DIR}/Methylococcales_heme_queries.faa"

QUERY_METADATA="${WORK_DIR}/Methylococcales_heme_query_metadata.tsv"


GLOBDB_PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"

GLOBDB_CLUSTER_DIR="${GLOBDB_PROJECT}/06_clustering"


## Search against ALL proteins that formed the original clusters,
## not only the cluster representatives.
TARGET_FASTA="${GLOBDB_CLUSTER_DIR}/exported_heme_candidates.faa"

OLD_MEMBERSHIP="${GLOBDB_CLUSTER_DIR}/cluster_membership.tsv"

OLD_CLUSTER_SUMMARY="${GLOBDB_CLUSTER_DIR}/cluster_summary.tsv"


RESULTS_TSV="${WORK_DIR}/Methylococcales_vs_GlobDB_all_members.tsv"

ASSIGNMENTS="${WORK_DIR}/Methylococcales_GlobDB_cluster_assignments.tsv"

PRIORITY_OUT="${WORK_DIR}/Methylococcales_priority_cluster_assignments.tsv"

ALL_HITS_OUT="${WORK_DIR}/Methylococcales_all_passing_hits_with_clusters.tsv"

CLUSTER_COUNTS="${WORK_DIR}/Methylococcales_cluster_mapping_summary.tsv"


TMP_DIR="${SLURM_TMPDIR:-${WORK_DIR}/tmp_${SLURM_JOB_ID}}"


###############################################################################
## Search criteria
###############################################################################

MIN_SEQ_ID="0.40"

COVERAGE="0.80"

COV_MODE="0"


###############################################################################
## Helpers
###############################################################################

timestamp() {

    date '+%Y-%m-%d %H:%M:%S'

}


die() {

    echo "[$(timestamp)] ERROR: $*" >&2

    exit 1

}


###############################################################################
## Validate inputs
###############################################################################

mkdir -p \
    "$WORK_DIR"


for path in \
    "$QUERY_FASTA" \
    "$QUERY_METADATA" \
    "$TARGET_FASTA" \
    "$OLD_MEMBERSHIP" \
    "$OLD_CLUSTER_SUMMARY"
do

    [[ -s "$path" ]] \
        || die \
        "Required file missing or empty: $path"

done


###############################################################################
## Activate MMseqs2 if necessary
###############################################################################

if ! command -v mmseqs >/dev/null 2>&1; then

    if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then

        source "$HOME/miniforge3/etc/profile.d/conda.sh"

    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then

        source "$HOME/miniconda3/etc/profile.d/conda.sh"

    else

        die \
        "mmseqs is not in PATH and conda.sh could not be found."

    fi


    conda activate MMseqs2 \
        || die \
        "Could not activate conda environment 'mmseqs2'."

fi


command -v mmseqs >/dev/null 2>&1 \
    || die \
    "mmseqs is not available."


command -v python3 >/dev/null 2>&1 \
    || die \
    "python3 is not available."


###############################################################################
## Report
###############################################################################

echo \
    "[$(timestamp)] Mapping Methylococcales proteins to fixed GlobDB clusters."

echo \
    "MMseqs executable: $(command -v mmseqs)"

echo \
    "MMseqs version:"

mmseqs version


echo \
    "Query proteins:    $(grep -c '^>' "$QUERY_FASTA")"

echo \
    "Target proteins:   $(grep -c '^>' "$TARGET_FASTA")"

echo \
    "Minimum identity:  $MIN_SEQ_ID"

echo \
    "Coverage:          $COVERAGE"

echo \
    "Coverage mode:     $COV_MODE"

echo


###############################################################################
## Refuse accidental overwrite
###############################################################################

if [[ -s "$RESULTS_TSV" ]]; then

    die \
    "Existing search result found: $RESULTS_TSV"

fi


rm -rf \
    "$TMP_DIR"


mkdir -p \
    "$TMP_DIR"


###############################################################################
## Search
###############################################################################

mmseqs easy-search \
    "$QUERY_FASTA" \
    "$TARGET_FASTA" \
    "$RESULTS_TSV" \
    "$TMP_DIR" \
    --min-seq-id "$MIN_SEQ_ID" \
    -c "$COVERAGE" \
    --cov-mode "$COV_MODE" \
    --max-seqs 10000 \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"


[[ -f "$RESULTS_TSV" ]] \
    || die \
    "MMseqs2 did not create the expected search result."


###############################################################################
## Convert old GlobDB hits into fixed cluster assignments
###############################################################################

python3 - \
    "$QUERY_METADATA" \
    "$OLD_MEMBERSHIP" \
    "$OLD_CLUSTER_SUMMARY" \
    "$RESULTS_TSV" \
    "$ALL_HITS_OUT" \
    "$ASSIGNMENTS" \
    "$PRIORITY_OUT" \
    "$CLUSTER_COUNTS" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path


(
    query_metadata_path,
    old_membership_path,
    old_cluster_summary_path,
    search_path,
    all_hits_out,
    assignments_out,
    priority_out,
    cluster_counts_out,
) = map(Path, sys.argv[1:])


###############################################################################
## Helpers
###############################################################################

def read_tsv(path):

    with path.open(
        "r",
        encoding="utf-8",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        return (
            list(
                reader
            ),
            reader.fieldnames
            or []
        )


def as_float(value):

    try:

        return float(
            value
        )

    except Exception:

        return float(
            "-inf"
        )


def as_int(value):

    value = str(
        value
    ).strip()

    if value in {
        "1",
        "TRUE",
        "True",
        "true"
    }:

        return 1

    return 0


###############################################################################
## Query metadata
###############################################################################

query_rows, query_fields = read_tsv(
    query_metadata_path
)


query_meta = {}


for row in query_rows:

    qid = row[
        "query_id"
    ].strip()


    if qid in query_meta:

        raise SystemExit(
            f"ERROR: duplicate query_id: {qid}"
        )


    query_meta[
        qid
    ] = row


###############################################################################
## Existing GlobDB cluster membership
###############################################################################

membership_rows, membership_fields = read_tsv(
    old_membership_path
)


required_membership = {
    "cluster",
    "representative_protein",
    "protein_id"
}


missing = required_membership.difference(
    membership_fields
)


if missing:

    raise SystemExit(

        "ERROR: old cluster membership missing columns: "
        f"{sorted(missing)}"
    )


target_to_cluster = {}

target_to_rep = {}


for row in membership_rows:

    target = row[
        "protein_id"
    ].strip()


    if target in target_to_cluster:

        raise SystemExit(

            f"ERROR: duplicate target protein in old "
            f"cluster membership: {target}"
        )


    target_to_cluster[
        target
    ] = row[
        "cluster"
    ].strip()


    target_to_rep[
        target
    ] = row[
        "representative_protein"
    ].strip()


###############################################################################
## Existing cluster summary
###############################################################################

cluster_summary_rows, cluster_summary_fields = read_tsv(
    old_cluster_summary_path
)


cluster_summary = {

    row[
        "cluster"
    ].strip():
        row

    for row
    in cluster_summary_rows
}


###############################################################################
## Parse MMseqs hits
###############################################################################

hits_by_query = defaultdict(
    list
)


all_enriched_hits = []


if search_path.stat().st_size > 0:

    with search_path.open(
        "r",
        encoding="utf-8"
    ) as handle:

        for line_no, raw in enumerate(
            handle,
            start=1
        ):

            line = raw.rstrip(
                "\r\n"
            )


            if not line:

                continue


            fields = line.split(
                "\t"
            )


            if len(
                fields
            ) != 8:

                raise SystemExit(

                    f"ERROR: malformed MMseqs result line "
                    f"{line_no}: {line!r}"
                )


            (
                query,
                target,
                pident,
                alnlen,
                qcov,
                tcov,
                evalue,
                bits
            ) = fields


            if query not in query_meta:

                raise SystemExit(

                    f"ERROR: MMseqs query not present "
                    f"in query metadata: {query}"
                )


            if target not in target_to_cluster:

                raise SystemExit(

                    f"ERROR: MMseqs target not present "
                    f"in old GlobDB cluster membership: {target}"
                )


            cluster = target_to_cluster[
                target
            ]


            hit = {

                "query_id":
                    query,

                "target_protein":
                    target,

                "cluster":
                    cluster,

                "cluster_representative":
                    target_to_rep[
                        target
                    ],

                "pident":
                    pident,

                "alignment_length":
                    alnlen,

                "query_coverage":
                    qcov,

                "target_coverage":
                    tcov,

                "evalue":
                    evalue,

                "bits":
                    bits
            }


            hits_by_query[
                query
            ].append(
                hit
            )


            enriched = {

                **query_meta[
                    query
                ],

                **hit
            }


            all_enriched_hits.append(
                enriched
            )


###############################################################################
## Rank hits
###############################################################################

def hit_rank(
    hit
):

    return (

        as_float(
            hit[
                "bits"
            ]
        ),

        as_float(
            hit[
                "pident"
            ]
        ),

        min(
            as_float(
                hit[
                    "query_coverage"
                ]
            ),
            as_float(
                hit[
                    "target_coverage"
                ]
            )
        )
    )


###############################################################################
## Build one assignment row per query
###############################################################################

assignment_rows = []


for qid in sorted(
    query_meta
):

    metadata = query_meta[
        qid
    ]


    hits = hits_by_query.get(
        qid,
        []
    )


    if not hits:

        row = {

            **metadata,

            "assignment_status":
                "no_match",

            "n_passing_hits":
                0,

            "n_clusters_hit":
                0,

            "top_scoring_cluster":
                "",

            "top_target_protein":
                "",

            "top_cluster_representative":
                "",

            "top_pident":
                "",

            "top_query_coverage":
                "",

            "top_target_coverage":
                "",

            "top_evalue":
                "",

            "top_bits":
                "",

            "all_clusters_hit":
                ""
        }


        assignment_rows.append(
            row
        )

        continue


    clusters = defaultdict(
        list
    )


    for hit in hits:

        clusters[
            hit[
                "cluster"
            ]
        ].append(
            hit
        )


    best_per_cluster = {}


    for cluster, cluster_hits in clusters.items():

        best_per_cluster[
            cluster
        ] = max(
            cluster_hits,
            key=hit_rank
        )


    ordered_clusters = sorted(

        best_per_cluster,

        key=lambda cluster:
            hit_rank(
                best_per_cluster[
                    cluster
                ]
            ),

        reverse=True
    )


    top_cluster = ordered_clusters[
        0
    ]


    top_hit = best_per_cluster[
        top_cluster
    ]


    status = (

        "unique_cluster"

        if len(
            ordered_clusters
        ) == 1

        else "multiple_clusters"
    )


    old_summary = cluster_summary.get(
        top_cluster,
        {}
    )


    row = {

        **metadata,

        "assignment_status":
            status,

        "n_passing_hits":
            len(
                hits
            ),

        "n_clusters_hit":
            len(
                ordered_clusters
            ),

        "top_scoring_cluster":
            top_cluster,

        "top_target_protein":
            top_hit[
                "target_protein"
            ],

        "top_cluster_representative":
            top_hit[
                "cluster_representative"
            ],

        "top_pident":
            top_hit[
                "pident"
            ],

        "top_query_coverage":
            top_hit[
                "query_coverage"
            ],

        "top_target_coverage":
            top_hit[
                "target_coverage"
            ],

        "top_evalue":
            top_hit[
                "evalue"
            ],

        "top_bits":
            top_hit[
                "bits"
            ],

        "all_clusters_hit":
            ";".join(
                ordered_clusters
            ),

        "old_cluster_n_proteins":
            old_summary.get(
                "n_proteins",
                ""
            ),

        "old_cluster_n_genomes":
            old_summary.get(
                "n_genomes",
                ""
            ),

        "old_cluster_mean_heme_count":
            old_summary.get(
                "mean_heme_count",
                ""
            ),

        "old_cluster_max_heme_count":
            old_summary.get(
                "max_heme_count",
                ""
            ),

        "old_cluster_dominant_fegenie_HMM":
            old_summary.get(
                "dominant_fegenie_HMM",
                ""
            ),

        "old_cluster_dominant_signalp_prediction":
            old_summary.get(
                "dominant_signalp_prediction",
                ""
            ),

        "old_cluster_dominant_deeptmhmm_class":
            old_summary.get(
                "dominant_deeptmhmm_class",
                ""
            ),

        "old_cluster_dominant_localization_class":
            old_summary.get(
                "dominant_localization_class",
                ""
            )
    }


    assignment_rows.append(
        row
    )


###############################################################################
## Output helpers
###############################################################################

def write_rows(
    path,
    rows
):

    if not rows:

        raise SystemExit(
            f"ERROR: no rows available for {path}"
        )


    fields = []


    seen = set()


    for row in rows:

        for field in row:

            if field not in seen:

                fields.append(
                    field
                )

                seen.add(
                    field
                )


    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(

            handle,

            fieldnames=
                fields,

            delimiter="\t",

            lineterminator="\n",

            extrasaction="ignore"
        )


        writer.writeheader()

        writer.writerows(
            rows
        )


###############################################################################
## Full passing-hit table
###############################################################################

if all_enriched_hits:

    all_enriched_hits.sort(

        key=lambda row: (
            row[
                "query_id"
            ],
            -as_float(
                row[
                    "bits"
                ]
            )
        )
    )


    write_rows(
        all_hits_out,
        all_enriched_hits
    )


else:

    with all_hits_out.open(
        "w",
        encoding="utf-8"
    ) as handle:

        handle.write(
            "query_id\n"
        )


###############################################################################
## Assignment table
###############################################################################

write_rows(
    assignments_out,
    assignment_rows
)


###############################################################################
## Priority subset
###############################################################################

priority_rows = [

    row

    for row
    in assignment_rows

    if row.get(
        "priority_reason",
        ""
    ).strip()
]


write_rows(
    priority_out,
    priority_rows
)


###############################################################################
## Cluster-level metagenome summary
###############################################################################

clusters_seen = defaultdict(

    lambda: {
        "queries_with_any_hit":
            set(),

        "unique_assignments":
            set(),

        "multi_cluster_top_hits":
            set(),

        "genomes_with_any_hit":
            set(),

        "cyc2_queries":
            set(),

        "ge5_heme_queries":
            set()
    }
)


for hit in all_enriched_hits:

    cluster = hit[
        "cluster"
    ]


    qid = hit[
        "query_id"
    ]


    clusters_seen[
        cluster
    ][
        "queries_with_any_hit"
    ].add(
        qid
    )


    clusters_seen[
        cluster
    ][
        "genomes_with_any_hit"
    ].add(
        hit[
            "genome"
        ]
    )


    if as_int(
        hit.get(
            "priority_cyc2",
            "0"
        )
    ):

        clusters_seen[
            cluster
        ][
            "cyc2_queries"
        ].add(
            qid
        )


    if as_int(
        hit.get(
            "priority_ge5_hemes",
            "0"
        )
    ):

        clusters_seen[
            cluster
        ][
            "ge5_heme_queries"
        ].add(
            qid
        )


for row in assignment_rows:

    cluster = row.get(
        "top_scoring_cluster",
        ""
    )


    if not cluster:

        continue


    if row[
        "assignment_status"
    ] == "unique_cluster":

        clusters_seen[
            cluster
        ][
            "unique_assignments"
        ].add(
            row[
                "query_id"
            ]
        )


    elif row[
        "assignment_status"
    ] == "multiple_clusters":

        clusters_seen[
            cluster
        ][
            "multi_cluster_top_hits"
        ].add(
            row[
                "query_id"
            ]
        )


cluster_rows = []


for cluster in sorted(
    clusters_seen
):

    values = clusters_seen[
        cluster
    ]


    old = cluster_summary.get(
        cluster,
        {}
    )


    cluster_rows.append({

        "cluster":
            cluster,

        "queries_with_any_hit":
            len(
                values[
                    "queries_with_any_hit"
                ]
            ),

        "genomes_with_any_hit":
            len(
                values[
                    "genomes_with_any_hit"
                ]
            ),

        "unique_assignments":
            len(
                values[
                    "unique_assignments"
                ]
            ),

        "multiple_cluster_top_hits":
            len(
                values[
                    "multi_cluster_top_hits"
                ]
            ),

        "Cyc2_priority_queries":
            len(
                values[
                    "cyc2_queries"
                ]
            ),

        "ge5_heme_priority_queries":
            len(
                values[
                    "ge5_heme_queries"
                ]
            ),

        "old_n_proteins":
            old.get(
                "n_proteins",
                ""
            ),

        "old_n_genomes":
            old.get(
                "n_genomes",
                ""
            ),

        "old_mean_heme_count":
            old.get(
                "mean_heme_count",
                ""
            ),

        "old_max_heme_count":
            old.get(
                "max_heme_count",
                ""
            ),

        "old_dominant_fegenie_HMM":
            old.get(
                "dominant_fegenie_HMM",
                ""
            ),

        "old_dominant_localization_class":
            old.get(
                "dominant_localization_class",
                ""
            )
    })


write_rows(
    cluster_counts_out,
    cluster_rows
)


###############################################################################
## Report
###############################################################################

n_no_match = sum(

    row[
        "assignment_status"
    ] == "no_match"

    for row
    in assignment_rows
)


n_unique = sum(

    row[
        "assignment_status"
    ] == "unique_cluster"

    for row
    in assignment_rows
)


n_multiple = sum(

    row[
        "assignment_status"
    ] == "multiple_clusters"

    for row
    in assignment_rows
)


print()

print(
    "GlobDB fixed-cluster mapping complete."
)

print(
    f"Queries:                 {len(assignment_rows)}"
)

print(
    f"Unique-cluster matches:  {n_unique}"
)

print(
    f"Multiple-cluster hits:   {n_multiple}"
)

print(
    f"No matching cluster:     {n_no_match}"
)

print(
    f"Priority proteins:       {len(priority_rows)}"
)

print()

print(
    f"Assignments: {assignments_out}"
)

print(
    f"Priority:    {priority_out}"
)

print(
    f"All hits:    {all_hits_out}"
)

print(
    f"Clusters:    {cluster_counts_out}"
)
PY


###############################################################################
## Completion marker
###############################################################################

cat > "${WORK_DIR}/GLOBDB_CLUSTER_MAPPING_COMPLETE.txt" <<EOF
GlobDB fixed-cluster mapping completed successfully.
completion_time=$(timestamp)
job_id=${SLURM_JOB_ID}
query_fasta=${QUERY_FASTA}
target_fasta=${TARGET_FASTA}
min_seq_id=${MIN_SEQ_ID}
coverage=${COVERAGE}
cov_mode=${COV_MODE}
EOF


echo
echo \
    "[$(timestamp)] Finished mapping Methylococcales proteins to GlobDB clusters."
