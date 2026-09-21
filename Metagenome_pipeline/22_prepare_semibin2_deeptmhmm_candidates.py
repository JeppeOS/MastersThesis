#!/usr/bin/env python3

from __future__ import annotations

import csv
import re

from collections import defaultdict
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


FEGENIE_RESULTS = (
    PROJECT
    / "functional_analysis"
    / "FeGenie_SemiBin2_full_run_20260918"
    / "fegenie_results"
)


ORF_DIR = (
    FEGENIE_RESULTS
    / "ORF_calls"
)


FEGENIE_CSV = (
    FEGENIE_RESULTS
    / "FeGenie-geneSummary.csv"
)


FINDMEHEMES_TSV = (
    PROJECT
    / "functional_analysis"
    / "findmehemes"
    / "findmehemes_all_hits.tsv"
)


SIGNALP_ROOT = (
    PROJECT
    / "functional_analysis"
    / "signalp"
    / "results_by_genome"
)


OUT_ROOT = (
    PROJECT
    / "functional_analysis"
    / "deeptmhmm"
)


CANDIDATE_DIR = (
    OUT_ROOT
    / "candidates_by_genome"
)


MANIFEST_TSV = (
    OUT_ROOT
    / "deeptmhmm_candidate_manifest.tsv"
)


GENOME_MANIFEST_TSV = (
    OUT_ROOT
    / "deeptmhmm_genome_manifest.tsv"
)


QC_TSV = (
    OUT_ROOT
    / "deeptmhmm_candidate_qc.tsv"
)


###############################################################################
## Helpers
###############################################################################

def normalise_genome(value: str) -> str:

    """
    Return the bare genome ID from FeGenie/FASTA-derived names.

    Examples:

        semibin2__flye__barcode10_seqs__SemiBin_0.fa-proteins.faa
        ->
        semibin2__flye__barcode10_seqs__SemiBin_0

        semibin2__flye__barcode10_seqs__SemiBin_0.fa
        ->
        semibin2__flye__barcode10_seqs__SemiBin_0
    """

    name = Path(
        value.strip()
    ).name


    suffixes = (
        "-proteins.faa",
        ".fasta",
        ".fna",
        ".faa",
        ".fa"
    )


    changed = True


    while changed:

        changed = False


        for suffix in suffixes:

            if name.endswith(
                suffix
            ):

                name = name[
                    :-len(
                        suffix
                    )
                ]

                changed = True

                break


    return name


def fasta_records(
    path: Path
):

    header = None

    seq_parts = []


    with path.open(
        "r",
        encoding="utf-8"
    ) as handle:

        for raw in handle:

            line = raw.rstrip(
                "\r\n"
            )


            if line.startswith(
                ">"
            ):

                if header is not None:

                    yield (
                        header,
                        "".join(
                            seq_parts
                        )
                    )


                header = line[1:]

                seq_parts = []


            else:

                seq_parts.append(
                    line.strip()
                )


    if header is not None:

        yield (
            header,
            "".join(
                seq_parts
            )
        )


def parse_signalp(
    path: Path
):

    """
    Parse SignalP 6 prediction_results.txt from Prodigal-derived FASTAs.

    SignalP preserves the complete Prodigal FASTA header before its
    prediction columns. Because that header contains spaces, Prediction
    is not simply the second whitespace-delimited field.
    """

    rows = {}


    prediction_pattern = re.compile(

        r"\s"
        r"(OTHER|SP|LIPO|TAT|TATLIPO|PILIN)"
        r"\s+([0-9.eE+-]+)"
        r"\s+([0-9.eE+-]+)"
        r"\s+([0-9.eE+-]+)"
        r"\s+([0-9.eE+-]+)"
        r"\s+([0-9.eE+-]+)"
        r"\s+([0-9.eE+-]+)"
        r"(?:\s+(.*))?$"
    )


    with path.open(
        "r",
        encoding="utf-8"
    ) as handle:

        for raw in handle:

            line = raw.rstrip(
                "\r\n"
            )


            if (
                not line
                or line.startswith(
                    "#"
                )
            ):

                continue


            protein_id = line.split()[0]


            match = prediction_pattern.search(
                line
            )


            if match is None:

                raise ValueError(

                    f"Could not parse SignalP line in "
                    f"{path}: {line[:200]!r}"
                )


            prediction = match.group(
                1
            )


            trailing = (
                match.group(
                    8
                )
                or ""
            ).strip()


            cs_match = re.search(

                r"CS pos:\s*(.*?)(?=(?:\s+Pr:)|$)",

                trailing
            )


            rows[
                protein_id
            ] = {

                "signalp_prediction":
                    prediction,

                "signalp_cs_position":
                    (
                        cs_match.group(
                            1
                        ).strip()
                        if cs_match
                        else ""
                    )
            }


    if not rows:

        raise ValueError(
            f"No SignalP predictions parsed from {path}"
        )


    return rows


###############################################################################
## Validate required inputs
###############################################################################

for path in (
    FEGENIE_CSV,
    FINDMEHEMES_TSV
):

    if not path.is_file():

        raise SystemExit(
            f"ERROR: missing required file: {path}"
        )


if not ORF_DIR.is_dir():

    raise SystemExit(
        f"ERROR: missing ORF directory: {ORF_DIR}"
    )


OUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)


CANDIDATE_DIR.mkdir(
    parents=True,
    exist_ok=True
)


###############################################################################
## Genome-name normalization self-test
###############################################################################

_normalization_tests = {

    "semibin2__flye__barcode10_seqs__SemiBin_0.fa-proteins.faa":
        "semibin2__flye__barcode10_seqs__SemiBin_0",

    "semibin2__flye__barcode10_seqs__SemiBin_0.fa":
        "semibin2__flye__barcode10_seqs__SemiBin_0",

    "semibin2__flye__barcode10_seqs__SemiBin_0":
        "semibin2__flye__barcode10_seqs__SemiBin_0"
}


for raw_name, expected_name in _normalization_tests.items():

    observed_name = normalise_genome(
        raw_name
    )


    if observed_name != expected_name:

        raise SystemExit(

            "ERROR: genome-name normalization self-test failed: "
            f"{raw_name!r} -> {observed_name!r}, "
            f"expected {expected_name!r}"
        )


###############################################################################
## Read FeGenie results
###############################################################################

fegenie = defaultdict(
    lambda: defaultdict(
        list
    )
)


with FEGENIE_CSV.open(
    "r",
    encoding="utf-8",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle
    )


    req = {
        "category",
        "genome/assembly",
        "orf",
        "HMM",
        "bitscore",
        "bitscore_cutoff"
    }


    missing = req.difference(
        reader.fieldnames
        or []
    )


    if missing:

        raise SystemExit(

            f"ERROR: FeGenie CSV missing columns: "
            f"{sorted(missing)}"
        )


    for row in reader:

        genome_raw = row[
            "genome/assembly"
        ].strip()


        pid = row[
            "orf"
        ].strip()


        ## FeGenie output can contain a repeated CSV header row.
        if (
            genome_raw
            in {
                "genome/assembly",
                "assembly"
            }
            and pid == "orf"
        ):

            continue


        if (
            not genome_raw
            or not pid
        ):

            continue


        genome = normalise_genome(
            genome_raw
        )


        fegenie[
            genome
        ][
            pid
        ].append(
            row
        )


###############################################################################
## Read FindMeHemes results
###############################################################################

hemes = defaultdict(
    dict
)


with FINDMEHEMES_TSV.open(
    "r",
    encoding="utf-8",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )


    req = {
        "genome",
        "protein_id",
        "number_of_hemes"
    }


    missing = req.difference(
        reader.fieldnames
        or []
    )


    if missing:

        raise SystemExit(

            f"ERROR: FindMeHemes TSV missing columns: "
            f"{sorted(missing)}"
        )


    for row in reader:

        genome = normalise_genome(
            row[
                "genome"
            ]
        )


        hemes[
            genome
        ][
            row[
                "protein_id"
            ].strip()
        ] = int(
            row[
                "number_of_hemes"
            ]
        )


###############################################################################
## Discover FeGenie proteomes
###############################################################################

proteomes = sorted(

    path

    for path in ORF_DIR.glob(
        "*-proteins.faa"
    )

    if (
        path.is_file()
        and path.stat().st_size > 0
    )
)


if len(proteomes) != 187:

    print(

        f"WARNING: expected 187 proteomes, "
        f"detected {len(proteomes)}"
    )


###############################################################################
## Prepare candidate FASTAs
###############################################################################

manifest_rows = []

genome_rows = []

qc_rows = []


seen_fg = set()

seen_heme = set()


for proteome in proteomes:

    genome = normalise_genome(
        proteome.name
    )


    sp_file = (
        SIGNALP_ROOT
        / genome
        / "prediction_results.txt"
    )


    sp = (
        parse_signalp(
            sp_file
        )
        if sp_file.is_file()
        else {}
    )


    out_faa = (
        CANDIDATE_DIR
        / f"{genome}.faa"
    )


    n_input = 0

    n_candidates = 0

    n_fg = 0

    n_heme = 0

    n_both = 0

    n_sp = 0


    with out_faa.open(
        "w",
        encoding="utf-8"
    ) as out:

        for full_header, seq in fasta_records(
            proteome
        ):

            n_input += 1


            pid = full_header.split()[0]


            fg_hits = (
                fegenie
                .get(
                    genome,
                    {}
                )
                .get(
                    pid,
                    []
                )
            )


            heme_count = (
                hemes
                .get(
                    genome,
                    {}
                )
                .get(
                    pid,
                    0
                )
            )


            is_fg = bool(
                fg_hits
            )


            is_heme = (
                heme_count > 0
            )


            ## Candidate if positive in FeGenie OR FindMeHemes.
            if not (
                is_fg
                or is_heme
            ):

                continue


            n_candidates += 1

            n_fg += int(
                is_fg
            )

            n_heme += int(
                is_heme
            )

            n_both += int(
                is_fg
                and is_heme
            )


            if is_fg:

                seen_fg.add(
                    (
                        genome,
                        pid
                    )
                )


            if is_heme:

                seen_heme.add(
                    (
                        genome,
                        pid
                    )
                )


            sp_row = sp.get(

                pid,

                {
                    "signalp_prediction":
                        "",

                    "signalp_cs_position":
                        ""
                }
            )


            n_sp += int(
                bool(
                    sp_row[
                        "signalp_prediction"
                    ]
                )
            )


            hmms = sorted(
                {
                    hit[
                        "HMM"
                    ].strip()

                    for hit
                    in fg_hits

                    if hit[
                        "HMM"
                    ].strip()
                }
            )


            cats = sorted(
                {
                    hit[
                        "category"
                    ].strip()

                    for hit
                    in fg_hits

                    if hit[
                        "category"
                    ].strip()
                }
            )


            bits = []

            cuts = []


            for hit in fg_hits:

                try:

                    bits.append(
                        float(
                            hit[
                                "bitscore"
                            ]
                        )
                    )

                except Exception:

                    pass


                try:

                    cuts.append(
                        float(
                            hit[
                                "bitscore_cutoff"
                            ]
                        )
                    )

                except Exception:

                    pass


            source = (

                "FeGenie+FindMeHemes"

                if (
                    is_fg
                    and is_heme
                )

                else (
                    "FeGenie"
                    if is_fg
                    else "FindMeHemes"
                )
            )


            out.write(
                f">{full_header}\n"
            )


            for i in range(
                0,
                len(
                    seq
                ),
                80
            ):

                out.write(
                    seq[
                        i:i + 80
                    ]
                    + "\n"
                )


            manifest_rows.append({

                "genome":
                    genome,

                "protein_id":
                    pid,

                "length":
                    len(
                        seq
                    ),

                "candidate_source":
                    source,

                "fegenie_positive":
                    int(
                        is_fg
                    ),

                "fegenie_HMMs":
                    ";".join(
                        hmms
                    ),

                "fegenie_categories":
                    ";".join(
                        cats
                    ),

                "fegenie_max_bitscore":
                    (
                        max(
                            bits
                        )
                        if bits
                        else ""
                    ),

                "fegenie_max_cutoff":
                    (
                        max(
                            cuts
                        )
                        if cuts
                        else ""
                    ),

                "findmehemes_positive":
                    int(
                        is_heme
                    ),

                "number_of_hemes":
                    heme_count,

                "signalp_prediction":
                    sp_row[
                        "signalp_prediction"
                    ],

                "signalp_cs_position":
                    sp_row[
                        "signalp_cs_position"
                    ],

                "candidate_fasta":
                    str(
                        out_faa
                    ),

                "source_proteome":
                    str(
                        proteome
                    )
            })


    genome_rows.append({

        "genome":
            genome,

        "candidate_fasta":
            str(
                out_faa
            ),

        "input_proteins":
            n_input,

        "candidate_proteins":
            n_candidates,

        "fegenie_candidates":
            n_fg,

        "findmehemes_candidates":
            n_heme,

        "both":
            n_both
    })


    qc_rows.append({

        "genome":
            genome,

        "signalp_file_found":
            int(
                sp_file.is_file()
            ),

        "candidate_proteins":
            n_candidates,

        "candidates_with_signalp_result":
            n_sp
    })


    if not sp_file.is_file():

        raise SystemExit(

            f"ERROR: SignalP result file missing "
            f"for {genome}: {sp_file}"
        )


    if n_sp != n_candidates:

        raise SystemExit(

            f"ERROR: SignalP join incomplete for "
            f"{genome}: "
            f"{n_sp}/{n_candidates} candidates matched."
        )


###############################################################################
## QC: make sure all FeGenie/FindMeHemes hits were recovered
###############################################################################

expected_fg = {

    (
        genome,
        protein
    )

    for genome, data
    in fegenie.items()

    for protein
    in data
}


expected_heme = {

    (
        genome,
        protein
    )

    for genome, data
    in hemes.items()

    for protein
    in data
}


missing_fg = sorted(
    expected_fg
    - seen_fg
)


missing_heme = sorted(
    expected_heme
    - seen_heme
)


if missing_fg:

    raise SystemExit(

        f"ERROR: {len(missing_fg)} FeGenie IDs "
        f"were not found in authoritative FASTAs. "
        f"First: {missing_fg[:10]}"
    )


if missing_heme:

    raise SystemExit(

        f"ERROR: {len(missing_heme)} FindMeHemes IDs "
        f"were not found in authoritative FASTAs. "
        f"First: {missing_heme[:10]}"
    )


###############################################################################
## Write manifests
###############################################################################

fields = [

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
    "candidate_fasta",
    "source_proteome"
]


with MANIFEST_TSV.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writeheader()

    writer.writerows(
        manifest_rows
    )


fields2 = [

    "genome",
    "candidate_fasta",
    "input_proteins",
    "candidate_proteins",
    "fegenie_candidates",
    "findmehemes_candidates",
    "both"
]


with GENOME_MANIFEST_TSV.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fields2,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writeheader()

    writer.writerows(
        genome_rows
    )


fields3 = [

    "genome",
    "signalp_file_found",
    "candidate_proteins",
    "candidates_with_signalp_result"
]


with QC_TSV.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fields3,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writeheader()

    writer.writerows(
        qc_rows
    )


###############################################################################
## Report
###############################################################################

print(
    "DeepTMHMM candidate preparation complete."
)


print(
    f"Proteomes:                    "
    f"{len(proteomes)}"
)


print(
    f"Total candidate proteins:     "
    f"{len(manifest_rows)}"
)


print(
    f"FeGenie only:                 "
    f"{sum(r['candidate_source'] == 'FeGenie' for r in manifest_rows)}"
)


print(
    f"FindMeHemes only:             "
    f"{sum(r['candidate_source'] == 'FindMeHemes' for r in manifest_rows)}"
)


print(
    f"Both:                         "
    f"{sum(r['candidate_source'] == 'FeGenie+FindMeHemes' for r in manifest_rows)}"
)


print(
    f"Candidates with SignalP data: "
    f"{sum(bool(r['signalp_prediction']) for r in manifest_rows)}"
)


print(
    f"Candidate FASTAs:             "
    f"{CANDIDATE_DIR}"
)


print(
    f"Protein manifest:             "
    f"{MANIFEST_TSV}"
)


print(
    f"Genome manifest:              "
    f"{GENOME_MANIFEST_TSV}"
)


print(
    f"QC table:                     "
    f"{QC_TSV}"
)
