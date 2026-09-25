#!/usr/bin/env python3

from __future__ import annotations

import csv
import re

from collections import Counter
from pathlib import Path


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


DEEPTMHMM_DIR = (
    PROJECT
    / "functional_analysis"
    / "deeptmhmm"
)


MANIFEST = (
    DEEPTMHMM_DIR
    / "deeptmhmm_candidate_manifest.tsv"
)


RESULTS_DIR = (
    DEEPTMHMM_DIR
    / "results_by_genome"
)


OUT_DIR = (
    PROJECT
    / "functional_analysis"
    / "integrated_annotations"
)


MASTER_OUT = (
    OUT_DIR
    / "MASTER_PROTEIN_ANNOTATIONS.tsv"
)


UNKNOWN_EXPORTED_OUT = (
    OUT_DIR
    / "unknown_exported_cytochrome_candidates.tsv"
)


UNKNOWN_SOLUBLE_OUT = (
    OUT_DIR
    / "unknown_soluble_periplasmic_like_cytochromes.tsv"
)


UNKNOWN_LIPO_OUT = (
    OUT_DIR
    / "unknown_lipoprotein_cytochromes.tsv"
)


UNKNOWN_MEMBRANE_OUT = (
    OUT_DIR
    / "unknown_membrane_associated_cytochromes.tsv"
)


SUMMARY_OUT = (
    OUT_DIR
    / "integration_summary.tsv"
)


CLASS_COUNTS_OUT = (
    OUT_DIR
    / "localization_class_counts.tsv"
)


###############################################################################
## SignalP classes
###############################################################################

SIGNALP_EXPORT_CLASSES = {
    "SP",
    "LIPO",
    "TAT",
    "TATLIPO",
    "PILIN"
}


SIGNALP_SOLUBLE_EXPORT_CLASSES = {
    "SP",
    "TAT"
}


SIGNALP_LIPO_CLASSES = {
    "LIPO",
    "TATLIPO"
}


###############################################################################
## Helper functions
###############################################################################

def count_runs(
    text: str,
    char: str
) -> int:

    return sum(
        1
        for _ in re.finditer(
            f"{re.escape(char)}+",
            text
        )
    )


def first_non_s(
    topology: str
) -> str:

    stripped = topology.lstrip(
        "S"
    )

    return (
        stripped[0]
        if stripped
        else ""
    )


def last_side(
    topology: str
) -> str:

    for char in reversed(
        topology
    ):

        if char in {
            "I",
            "O"
        }:

            return char

    return ""


###############################################################################
## Parse DeepTMHMM three-line output
###############################################################################

def parse_deeptmhmm_3line(
    path: Path
):

    lines = []


    with path.open(
        "r",
        encoding="utf-8"
    ) as handle:

        for raw in handle:

            line = raw.rstrip(
                "\r\n"
            )

            if line.strip():

                lines.append(
                    line
                )


    if len(lines) % 3 != 0:

        raise ValueError(

            f"{path} contains "
            f"{len(lines)} non-empty lines; "
            f"expected a multiple of 3."
        )


    for i in range(
        0,
        len(lines),
        3
    ):

        header = lines[
            i
        ]

        sequence = lines[
            i + 1
        ]

        topology = lines[
            i + 2
        ]


        if not header.startswith(
            ">"
        ):

            raise ValueError(

                f"Unexpected DeepTMHMM header "
                f"in {path}: {header!r}"
            )


        header_text = (
            header[1:]
            .strip()
        )


        if "|" not in header_text:

            raise ValueError(

                f"DeepTMHMM header lacks "
                f"class separator '|': "
                f"{header!r}"
            )


        protein_part, class_part = (
            header_text.rsplit(
                "|",
                1
            )
        )


        protein_id = (
            protein_part
            .strip()
            .split()[0]
        )


        deeptmhmm_class = (
            class_part
            .strip()
        )


        sequence_clean = (
            sequence.strip()
        )


        topology_clean = (
            topology.strip()
        )


        ## DeepTMHMM output from this installation retains the terminal
        ## "*" from the Prodigal translation in both the emitted sequence
        ## and topology strings.

        if (
            len(
                sequence_clean
            )
            !=
            len(
                topology_clean
            )
        ):

            raise ValueError(

                f"Length mismatch for "
                f"{protein_id} in {path}: "
                f"sequence={len(sequence_clean)}, "
                f"topology={len(topology_clean)}"
            )


        yield {

            "protein_id":
                protein_id,

            "deeptmhmm_class":
                deeptmhmm_class,

            "deeptmhmm_sequence":
                sequence_clean,

            "deeptmhmm_topology":
                topology_clean,

            "deeptmhmm_n_tm_helices":
                count_runs(
                    topology_clean,
                    "M"
                ),

            "deeptmhmm_sp_length":
                (
                    len(
                        topology_clean
                    )
                    -
                    len(
                        topology_clean.lstrip(
                            "S"
                        )
                    )
                ),

            "deeptmhmm_side_after_sp":
                first_non_s(
                    topology_clean
                ),

            "deeptmhmm_c_terminal_side":
                last_side(
                    topology_clean
                )
        }


###############################################################################
## Localization classification
###############################################################################

def localization_class(
    signalp_prediction: str,
    deeptmhmm_class: str,
    n_tm: int
) -> str:

    sp = (
        signalp_prediction
        .strip()
        .upper()
    )


    dt = (
        deeptmhmm_class
        .strip()
        .upper()
    )


    signalp_export = (
        sp
        in SIGNALP_EXPORT_CLASSES
    )


    deep_sp = (
        dt == "SP"
    )


    ## Lipoproteins are handled first because their membrane association
    ## is encoded by lipidation rather than by a conventional TM helix.
    if sp in SIGNALP_LIPO_CLASSES:

        return (
            "exported_lipoprotein_candidate"
        )


    if sp == "PILIN":

        return (
            "pilin_like_exported_candidate"
        )


    ## IMPORTANT:
    ## DeepTMHMM BETA denotes a beta-barrel membrane protein.
    ## These proteins generally have no alpha-helical "M" segments, so
    ## deeptmhmm_n_tm_helices == 0 must NOT cause them to be classified
    ## as soluble/periplasmic.
    if dt == "BETA":

        return (
            "beta_barrel_membrane_candidate"
        )


    ## Soluble exported/periplasmic-like proteins.
    if (
        (
            sp
            in SIGNALP_SOLUBLE_EXPORT_CLASSES

            or deep_sp
        )
        and n_tm == 0
    ):

        return (
            "soluble_exported_periplasmic_like_candidate"
        )


    ## Exported proteins that also contain alpha-helical TM segments.
    if (
        signalp_export
        or deep_sp
    ):

        if n_tm == 1:

            return (
                "exported_single_pass_membrane_candidate"
            )


        if n_tm >= 2:

            return (
                "exported_multipass_membrane_candidate"
            )


        return (
            "exported_candidate_uncertain_topology"
        )


    ## Membrane proteins without export evidence.
    if n_tm == 1:

        return (
            "single_pass_membrane_no_export_signal"
        )


    if n_tm >= 2:

        return (
            "multipass_membrane_no_export_signal"
        )


    if dt == "GLOB":

        return (
            "globular_no_export_signal"
        )


    return "uncertain"


def bool01(
    value
) -> int:

    return (

        1

        if str(
            value
        ).strip()
        in {
            "1",
            "TRUE",
            "True",
            "true"
        }

        else 0
    )


###############################################################################
## Validate input
###############################################################################

if not MANIFEST.is_file():

    raise SystemExit(

        f"ERROR: candidate manifest "
        f"not found: {MANIFEST}"
    )


if not RESULTS_DIR.is_dir():

    raise SystemExit(

        f"ERROR: DeepTMHMM results "
        f"directory not found: "
        f"{RESULTS_DIR}"
    )


OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


###############################################################################
## Read candidate manifest
###############################################################################

with MANIFEST.open(
    "r",
    encoding="utf-8",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )


    manifest_rows = list(
        reader
    )


    manifest_fields = (
        reader.fieldnames
        or []
    )


required_manifest = {

    "genome",
    "protein_id",
    "length",
    "candidate_source",
    "fegenie_positive",
    "fegenie_HMMs",
    "findmehemes_positive",
    "number_of_hemes",
    "signalp_prediction"
}


missing = (
    required_manifest
    .difference(
        manifest_fields
    )
)


if missing:

    raise SystemExit(

        f"ERROR: candidate manifest "
        f"missing required columns: "
        f"{sorted(missing)}"
    )


###############################################################################
## Index manifest by genome + protein
###############################################################################

manifest_by_key = {}


for row in manifest_rows:

    key = (

        row[
            "genome"
        ].strip(),

        row[
            "protein_id"
        ].strip()
    )


    if key in manifest_by_key:

        raise SystemExit(

            f"ERROR: duplicate "
            f"manifest key: {key}"
        )


    manifest_by_key[
        key
    ] = row


###############################################################################
## Parse DeepTMHMM results
###############################################################################

deep_by_key = {}


result_files = sorted(

    RESULTS_DIR.glob(
        "*/predicted_topologies.3line"
    )
)


## Some MAGs may have zero FeGenie/FindMeHemes candidates.
## Those receive DEEPTMHMM_NO_CANDIDATES.txt instead of a
## predicted_topologies.3line file.

if len(result_files) != 187:

    print(

        f"WARNING: found "
        f"{len(result_files)} "
        f"predicted_topologies.3line "
        f"files; expected 187."
    )


for result_file in result_files:

    genome = (
        result_file
        .parent
        .name
    )


    for deep_row in parse_deeptmhmm_3line(
        result_file
    ):

        key = (

            genome,

            deep_row[
                "protein_id"
            ]
        )


        if key in deep_by_key:

            raise SystemExit(

                f"ERROR: duplicate "
                f"DeepTMHMM result key: "
                f"{key}"
            )


        deep_by_key[
            key
        ] = deep_row


###############################################################################
## Verify exact candidate correspondence
###############################################################################

manifest_keys = set(
    manifest_by_key
)


deep_keys = set(
    deep_by_key
)


missing_deep = sorted(

    manifest_keys
    - deep_keys
)


unexpected_deep = sorted(

    deep_keys
    - manifest_keys
)


if missing_deep:

    preview = "\n".join(

        f"  {genome}\t{protein}"

        for genome, protein
        in missing_deep[:20]
    )


    raise SystemExit(

        "ERROR: candidate proteins are "
        "missing DeepTMHMM results.\n"

        f"Missing count: "
        f"{len(missing_deep)}\n"

        f"{preview}"
    )


if unexpected_deep:

    preview = "\n".join(

        f"  {genome}\t{protein}"

        for genome, protein
        in unexpected_deep[:20]
    )


    raise SystemExit(

        "ERROR: DeepTMHMM results contain "
        "proteins not present in the "
        "candidate manifest.\n"

        f"Unexpected count: "
        f"{len(unexpected_deep)}\n"

        f"{preview}"
    )


###############################################################################
## Integrate annotations
###############################################################################

final_rows = []


for row in manifest_rows:

    genome = (
        row[
            "genome"
        ].strip()
    )


    protein_id = (
        row[
            "protein_id"
        ].strip()
    )


    deep = deep_by_key[
        (
            genome,
            protein_id
        )
    ]


    signalp_prediction = (
        row
        .get(
            "signalp_prediction",
            ""
        )
        .strip()
        .upper()
    )


    deeptmhmm_class = (
        deep[
            "deeptmhmm_class"
        ]
        .strip()
        .upper()
    )


    n_tm = int(

        deep[
            "deeptmhmm_n_tm_helices"
        ]
    )


    signalp_export = int(

        signalp_prediction
        in SIGNALP_EXPORT_CLASSES
    )


    deeptmhmm_sp = int(

        deeptmhmm_class
        == "SP"
    )


    export_evidence = int(

        signalp_export
        or deeptmhmm_sp
    )


    loc_class = localization_class(

        signalp_prediction=
            signalp_prediction,

        deeptmhmm_class=
            deeptmhmm_class,

        n_tm=
            n_tm
    )


    fegenie_positive = bool01(

        row.get(
            "fegenie_positive",
            "0"
        )
    )


    heme_positive = bool01(

        row.get(
            "findmehemes_positive",
            "0"
        )
    )


    ## FindMeHemes-positive but FeGenie-negative.
    unknown_fegenie_cytochrome = int(

        fegenie_positive == 0

        and heme_positive == 1
    )


    unknown_exported_cytochrome = int(

        unknown_fegenie_cytochrome == 1

        and export_evidence == 1
    )


    ## Strong soluble/periplasmic-like candidate:
    ##
    ## - FindMeHemes positive
    ## - FeGenie negative
    ## - SignalP predicts SP/TAT
    ## - DeepTMHMM independently predicts SP
    ## - no alpha-helical TM segments
    ##
    ## DeepTMHMM BETA proteins do not enter this category because
    ## deeptmhmm_class must explicitly equal SP.
    strong_unknown_soluble = int(

        unknown_fegenie_cytochrome == 1

        and (
            signalp_prediction
            in SIGNALP_SOLUBLE_EXPORT_CLASSES
        )

        and deeptmhmm_class == "SP"

        and n_tm == 0
    )


    merged = dict(
        row
    )


    merged.update(
        deep
    )


    merged.update({

        "signalp_export_positive":
            signalp_export,

        "deeptmhmm_sp_positive":
            deeptmhmm_sp,

        "export_evidence":
            export_evidence,

        "localization_class":
            loc_class,

        "fegenie_unannotated_heme_candidate":
            unknown_fegenie_cytochrome,

        "fegenie_unannotated_exported_heme_candidate":
            unknown_exported_cytochrome,

        "strong_fegenie_unannotated_soluble_periplasmic_like_candidate":
            strong_unknown_soluble
    })


    final_rows.append(
        merged
    )


###############################################################################
## Output columns
###############################################################################

deep_fields = [

    "deeptmhmm_class",
    "deeptmhmm_sequence",
    "deeptmhmm_topology",
    "deeptmhmm_n_tm_helices",
    "deeptmhmm_sp_length",
    "deeptmhmm_side_after_sp",
    "deeptmhmm_c_terminal_side"
]


derived_fields = [

    "signalp_export_positive",
    "deeptmhmm_sp_positive",
    "export_evidence",
    "localization_class",
    "fegenie_unannotated_heme_candidate",
    "fegenie_unannotated_exported_heme_candidate",
    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate"
]


output_fields = (

    manifest_fields

    + [

        field

        for field
        in (
            deep_fields
            + derived_fields
        )

        if field
        not in manifest_fields
    ]
)


###############################################################################
## TSV writer
###############################################################################

def write_tsv(
    path: Path,
    rows
):

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(

            handle,

            fieldnames=
                output_fields,

            delimiter="\t",

            lineterminator="\n",

            extrasaction="ignore"
        )


        writer.writeheader()

        writer.writerows(
            rows
        )


###############################################################################
## Write master table
###############################################################################

write_tsv(
    MASTER_OUT,
    final_rows
)


###############################################################################
## Candidate subsets
###############################################################################

unknown_exported = [

    row

    for row
    in final_rows

    if int(
        row[
            "fegenie_unannotated_exported_heme_candidate"
        ]
    ) == 1
]


unknown_soluble = [

    row

    for row
    in final_rows

    if (
        row[
            "localization_class"
        ]
        ==
        "soluble_exported_periplasmic_like_candidate"

        and int(
            row[
                "fegenie_unannotated_heme_candidate"
            ]
        ) == 1
    )
]


unknown_lipo = [

    row

    for row
    in final_rows

    if (
        row[
            "localization_class"
        ]
        ==
        "exported_lipoprotein_candidate"

        and int(
            row[
                "fegenie_unannotated_heme_candidate"
            ]
        ) == 1
    )
]


## Membrane-associated unknown cytochrome candidates.
##
## This now includes beta-barrel proteins in addition to
## alpha-helical exported membrane proteins.
unknown_membrane = [

    row

    for row
    in final_rows

    if (
        row[
            "localization_class"
        ]
        in {
            "beta_barrel_membrane_candidate",
            "exported_single_pass_membrane_candidate",
            "exported_multipass_membrane_candidate"
        }

        and int(
            row[
                "fegenie_unannotated_heme_candidate"
            ]
        ) == 1
    )
]


write_tsv(
    UNKNOWN_EXPORTED_OUT,
    unknown_exported
)


write_tsv(
    UNKNOWN_SOLUBLE_OUT,
    unknown_soluble
)


write_tsv(
    UNKNOWN_LIPO_OUT,
    unknown_lipo
)


write_tsv(
    UNKNOWN_MEMBRANE_OUT,
    unknown_membrane
)


###############################################################################
## Localization counts
###############################################################################

localization_counts = Counter(

    row[
        "localization_class"
    ]

    for row
    in final_rows
)


with CLASS_COUNTS_OUT.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writerow(
        [
            "localization_class",
            "protein_count"
        ]
    )


    for class_name, count in sorted(
        localization_counts.items()
    ):

        writer.writerow(
            [
                class_name,
                count
            ]
        )


###############################################################################
## Summary
###############################################################################

summary_rows = [

    (
        "total_integrated_candidates",
        len(
            final_rows
        )
    ),

    (
        "fegenie_positive",

        sum(
            bool01(
                row[
                    "fegenie_positive"
                ]
            )

            for row
            in final_rows
        )
    ),

    (
        "findmehemes_positive",

        sum(
            bool01(
                row[
                    "findmehemes_positive"
                ]
            )

            for row
            in final_rows
        )
    ),

    (
        "signalp_export_positive",

        sum(
            int(
                row[
                    "signalp_export_positive"
                ]
            )

            for row
            in final_rows
        )
    ),

    (
        "deeptmhmm_sp_positive",

        sum(
            int(
                row[
                    "deeptmhmm_sp_positive"
                ]
            )

            for row
            in final_rows
        )
    ),

    (
        "fegenie_unannotated_heme_candidates",

        sum(
            int(
                row[
                    "fegenie_unannotated_heme_candidate"
                ]
            )

            for row
            in final_rows
        )
    ),

    (
        "fegenie_unannotated_exported_heme_candidates",
        len(
            unknown_exported
        )
    ),

    (
        "fegenie_unannotated_soluble_periplasmic_like_candidates",
        len(
            unknown_soluble
        )
    ),

    (
        "strong_fegenie_unannotated_soluble_periplasmic_like_candidates",

        sum(
            int(
                row[
                    "strong_fegenie_unannotated_soluble_periplasmic_like_candidate"
                ]
            )

            for row
            in final_rows
        )
    ),

    (
        "fegenie_unannotated_lipoprotein_candidates",
        len(
            unknown_lipo
        )
    ),

    (
        "fegenie_unannotated_membrane_associated_candidates",
        len(
            unknown_membrane
        )
    )
]


with SUMMARY_OUT.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writerow(
        [
            "metric",
            "value"
        ]
    )


    writer.writerows(
        summary_rows
    )


###############################################################################
## Report
###############################################################################

print(
    "Integrated annotation build complete."
)


print(
    f"DeepTMHMM result files:                              "
    f"{len(result_files)}"
)


print(
    f"Integrated candidate proteins:                       "
    f"{len(final_rows)}"
)


print(
    f"FeGenie-unannotated heme candidates:                 "
    f"{summary_rows[5][1]}"
)


print(
    f"FeGenie-unannotated exported heme candidates:        "
    f"{len(unknown_exported)}"
)


print(
    f"FeGenie-unannotated soluble/periplasmic-like:        "
    f"{len(unknown_soluble)}"
)


print(
    "Strong FeGenie-unannotated "
    "soluble/periplasmic-like: "
    f"{summary_rows[8][1]}"
)


print(
    f"FeGenie-unannotated lipoprotein candidates:          "
    f"{len(unknown_lipo)}"
)


print(
    f"FeGenie-unannotated membrane-associated candidates:  "
    f"{len(unknown_membrane)}"
)


print()


print(
    f"Master table: "
    f"{MASTER_OUT}"
)


print(
    f"Summary:      "
    f"{SUMMARY_OUT}"
)


print(
    f"Class counts: "
    f"{CLASS_COUNTS_OUT}"
)
