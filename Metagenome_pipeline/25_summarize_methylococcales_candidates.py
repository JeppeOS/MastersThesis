#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re


###############################################################################
## Paths
###############################################################################

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

METHYLO_TABLE = (
    PROJECT
    / "semibin2_analysis"
    / "Methylococcales_MAGs_SemiBin2.tsv"
)

MASTER = (
    PROJECT
    / "functional_analysis"
    / "integrated_annotations"
    / "MASTER_PROTEIN_ANNOTATIONS.tsv"
)

OUTDIR = (
    PROJECT
    / "functional_analysis"
    / "integrated_annotations"
    / "Methylococcales"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


###############################################################################
## Read input
###############################################################################

mags = pd.read_csv(
    METHYLO_TABLE,
    sep="\t",
    dtype=str
)

master = pd.read_csv(
    MASTER,
    sep="\t",
    dtype=str
)


###############################################################################
## Determine MAG ID column
###############################################################################

possible_id_cols = [
    "Genome_ID",
    "genome",
    "MAG_ID",
    "Bin_ID"
]

id_col = None

for col in possible_id_cols:

    if col in mags.columns:

        id_col = col

        break


if id_col is None:

    raise SystemExit(
        "ERROR: Could not identify MAG ID column in "
        f"{METHYLO_TABLE}\n"
        f"Columns: {list(mags.columns)}"
    )


if "genome" not in master.columns:

    raise SystemExit(
        "ERROR: MASTER_PROTEIN_ANNOTATIONS.tsv "
        "does not contain a 'genome' column."
    )


methylococcales_ids = set(
    mags[id_col]
    .dropna()
    .astype(str)
    .str.strip()
)


if len(methylococcales_ids) != 16:

    print(
        f"WARNING: expected 16 Methylococcales MAGs, "
        f"found {len(methylococcales_ids)} in metadata table."
    )


###############################################################################
## Extract Methylococcales proteins
###############################################################################

sub = master[
    master["genome"].isin(
        methylococcales_ids
    )
].copy()


found_ids = set(
    sub["genome"]
)


missing_ids = sorted(
    methylococcales_ids
    - found_ids
)


if missing_ids:

    print()
    print(
        "WARNING: Methylococcales MAGs with no "
        "integrated candidate proteins:"
    )

    for genome in missing_ids:

        print(
            f"  {genome}"
        )


###############################################################################
## Convert numeric columns
###############################################################################

numeric_cols = [
    "length",
    "fegenie_positive",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_export_positive",
    "deeptmhmm_sp_positive",
    "export_evidence",
    "deeptmhmm_n_tm_helices"
]


for col in numeric_cols:

    if col in sub.columns:

        sub[col] = pd.to_numeric(
            sub[col],
            errors="coerce"
        ).fillna(0)


###############################################################################
## Attach MAG metadata
###############################################################################

metadata = mags.copy()

metadata = metadata.rename(
    columns={
        id_col: "genome"
    }
)


full = sub.merge(
    metadata,
    on="genome",
    how="left",
    suffixes=("", "_MAG")
)


full.to_csv(
    OUTDIR
    / "Methylococcales_MASTER_PROTEIN_ANNOTATIONS.tsv",
    sep="\t",
    index=False
)


###############################################################################
## Per-MAG summary
###############################################################################

summary_rows = []


for genome in sorted(
    methylococcales_ids
):

    d = sub[
        sub["genome"] == genome
    ].copy()


    row = {
        "genome":
            genome,

        "candidate_proteins":
            len(d),

        "fegenie_positive":
            int(
                d["fegenie_positive"].sum()
            )
            if "fegenie_positive" in d
            else 0,

        "findmehemes_positive":
            int(
                d["findmehemes_positive"].sum()
            )
            if "findmehemes_positive" in d
            else 0,

        "heme_ge_2":
            int(
                (
                    d["number_of_hemes"]
                    >= 2
                ).sum()
            )
            if "number_of_hemes" in d
            else 0,

        "heme_ge_5":
            int(
                (
                    d["number_of_hemes"]
                    >= 5
                ).sum()
            )
            if "number_of_hemes" in d
            else 0,

        "heme_ge_10":
            int(
                (
                    d["number_of_hemes"]
                    >= 10
                ).sum()
            )
            if "number_of_hemes" in d
            else 0,

        "export_evidence":
            int(
                d["export_evidence"].sum()
            )
            if "export_evidence" in d
            else 0,

        "signalp_export":
            int(
                d["signalp_export_positive"].sum()
            )
            if "signalp_export_positive" in d
            else 0,

        "deeptmhmm_SP":
            int(
                d["deeptmhmm_sp_positive"].sum()
            )
            if "deeptmhmm_sp_positive" in d
            else 0,
    }


    if "localization_class" in d:

        row[
            "soluble_exported_periplasmic_like"
        ] = int(
            (
                d["localization_class"]
                ==
                "soluble_exported_periplasmic_like_candidate"
            ).sum()
        )


        row[
            "exported_lipoprotein"
        ] = int(
            (
                d["localization_class"]
                ==
                "exported_lipoprotein_candidate"
            ).sum()
        )


        row[
            "exported_membrane"
        ] = int(
            d["localization_class"].isin(
                [
                    "exported_single_pass_membrane_candidate",
                    "exported_multipass_membrane_candidate"
                ]
            ).sum()
        )


    summary_rows.append(
        row
    )


summary = pd.DataFrame(
    summary_rows
)


###############################################################################
## Add useful MAG metadata to summary
###############################################################################

meta_keep = [
    col
    for col in metadata.columns
    if col != "genome"
]


summary = summary.merge(
    metadata[
        ["genome"]
        + meta_keep
    ],
    on="genome",
    how="left"
)


summary.to_csv(
    OUTDIR
    / "Methylococcales_MAG_candidate_summary.tsv",
    sep="\t",
    index=False
)


###############################################################################
## FeGenie HMM counts
###############################################################################

if "fegenie_HMMs" in sub.columns:

    hmm_rows = []


    for _, row in sub.iterrows():

        value = str(
            row.get(
                "fegenie_HMMs",
                ""
            )
        ).strip()


        if (
            not value
            or value.lower() == "nan"
        ):

            continue


        for hmm in value.split(
            ";"
        ):

            hmm = hmm.strip()


            if not hmm:

                continue


            hmm_rows.append(
                {
                    "genome":
                        row["genome"],

                    "HMM":
                        hmm,

                    "protein_id":
                        row["protein_id"]
                }
            )


    hmm_df = pd.DataFrame(
        hmm_rows
    )


    if not hmm_df.empty:

        hmm_counts = (
            hmm_df
            .groupby(
                [
                    "genome",
                    "HMM"
                ],
                as_index=False
            )
            .agg(
                protein_count=(
                    "protein_id",
                    "nunique"
                )
            )
            .sort_values(
                [
                    "genome",
                    "protein_count",
                    "HMM"
                ],
                ascending=[
                    True,
                    False,
                    True
                ]
            )
        )


        hmm_counts.to_csv(
            OUTDIR
            / "Methylococcales_FeGenie_HMM_counts.tsv",
            sep="\t",
            index=False
        )


        overall_hmms = (
            hmm_df
            .groupby(
                "HMM",
                as_index=False
            )
            .agg(
                protein_count=(
                    "protein_id",
                    "nunique"
                ),

                MAG_count=(
                    "genome",
                    "nunique"
                )
            )
            .sort_values(
                [
                    "MAG_count",
                    "protein_count",
                    "HMM"
                ],
                ascending=[
                    False,
                    False,
                    True
                ]
            )
        )


        overall_hmms.to_csv(
            OUTDIR
            / "Methylococcales_FeGenie_HMM_overall.tsv",
            sep="\t",
            index=False
        )


###############################################################################
## All FindMeHemes-positive proteins sorted by heme count
###############################################################################

if (
    "findmehemes_positive" in sub.columns
    and
    "number_of_hemes" in sub.columns
):

    heme = sub[
        sub[
            "findmehemes_positive"
        ] == 1
    ].copy()


    heme = heme.sort_values(
        [
            "number_of_hemes",
            "genome"
        ],
        ascending=[
            False,
            True
        ]
    )


    heme.to_csv(
        OUTDIR
        / "Methylococcales_heme_candidates.tsv",
        sep="\t",
        index=False
    )


###############################################################################
## Broad EET / cytochrome keyword subset
###############################################################################

search_cols = [
    col
    for col in [
        "fegenie_HMMs",
        "fegenie_categories"
    ]
    if col in sub.columns
]


if search_cols:

    combined_text = (
        sub[
            search_cols
        ]
        .fillna("")
        .astype(str)
        .agg(
            " ; ".join,
            axis=1
        )
    )


    eet_pattern = re.compile(
        r"cyc2|mto[a-z0-9]*|mtr[a-z0-9]*|"
        r"omc[a-z0-9]*|pio[a-z0-9]*|"
        r"porin|cytochrome|iron.oxid|iron.reduc",
        flags=re.IGNORECASE
    )


    targeted = sub[
        combined_text.str.contains(
            eet_pattern,
            regex=True
        )
    ].copy()


    targeted = targeted.sort_values(
        [
            "genome",
            "fegenie_HMMs",
            "protein_id"
        ]
    )


    targeted.to_csv(
        OUTDIR
        / "Methylococcales_targeted_EET_FeGenie_hits.tsv",
        sep="\t",
        index=False
    )


###############################################################################
## Console report
###############################################################################

print()
print(
    "=" * 80
)

print(
    "Methylococcales integrated-candidate summary"
)

print(
    "=" * 80
)

print(
    f"Methylococcales MAGs:        "
    f"{len(methylococcales_ids)}"
)

print(
    f"MAGs with candidates:        "
    f"{sub['genome'].nunique()}"
)

print(
    f"Integrated candidates:       "
    f"{len(sub)}"
)

print(
    f"FeGenie-positive proteins:   "
    f"{int(sub['fegenie_positive'].sum())}"
)

print(
    f"FindMeHemes-positive:        "
    f"{int(sub['findmehemes_positive'].sum())}"
)

print(
    f"Proteins with >=2 hemes:     "
    f"{int((sub['number_of_hemes'] >= 2).sum())}"
)

print(
    f"Proteins with >=5 hemes:     "
    f"{int((sub['number_of_hemes'] >= 5).sum())}"
)

print(
    f"Proteins with >=10 hemes:    "
    f"{int((sub['number_of_hemes'] >= 10).sum())}"
)

print(
    f"Export evidence:             "
    f"{int(sub['export_evidence'].sum())}"
)

print()

print(
    "Per-MAG summary:"
)

display_cols = [
    "genome",
    "candidate_proteins",
    "fegenie_positive",
    "findmehemes_positive",
    "heme_ge_5",
    "heme_ge_10",
    "soluble_exported_periplasmic_like",
    "exported_lipoprotein",
    "exported_membrane"
]


display_cols = [
    col
    for col in display_cols
    if col in summary.columns
]


print(
    summary[
        display_cols
    ].to_string(
        index=False
    )
)


if (
    "fegenie_HMMs" in sub.columns
    and
    "hmm_df" in globals()
    and
    not hmm_df.empty
):

    print()
    print(
        "FeGenie HMMs present across the "
        "Methylococcales MAGs:"
    )

    print(
        overall_hmms.to_string(
            index=False
        )
    )


print()
print(
    f"Outputs written to:\n{OUTDIR}"
)

print(
    "=" * 80
)
