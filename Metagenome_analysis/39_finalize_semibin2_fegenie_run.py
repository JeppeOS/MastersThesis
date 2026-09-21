#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import shutil
from datetime import datetime


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

RUN_ROOT = (
    PROJECT
    / "functional_analysis"
    / "FeGenie_SemiBin2_full_run_20260918"
)

FEGENIE_OUT = (
    RUN_ROOT
    / "fegenie_results"
)

ORF_DIR = (
    FEGENIE_OUT
    / "ORF_calls"
)

FAA_OUT = (
    RUN_ROOT
    / "predicted_proteins_faa"
)

METADATA_OUT = (
    RUN_ROOT
    / "run_metadata"
)

MAG_MANIFEST = (
    PROJECT
    / "functional_analysis"
    / "SemiBin2_187_MAG_manifest.tsv"
)

PROTEIN_MANIFEST = (
    METADATA_OUT
    / "predicted_protein_manifest.tsv"
)

COMPLETE = (
    RUN_ROOT
    / "RUN_COMPLETED.txt"
)

EXPECTED = 187


## ------------------------------------------------------------
## Validate FeGenie outputs
## ------------------------------------------------------------

required_outputs = [
    FEGENIE_OUT / "FeGenie-geneSummary.csv",
    FEGENIE_OUT / "FeGenie-geneSummary-clusters.csv",
    FEGENIE_OUT / "FeGenie-heatmap-data.csv",
]


if not ORF_DIR.is_dir():

    raise RuntimeError(
        f"Missing FeGenie ORF_calls directory:\n"
        f"{ORF_DIR}"
    )


for path in required_outputs:

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Missing or empty FeGenie output:\n"
            f"{path}"
        )


## ------------------------------------------------------------
## Read authoritative SemiBin2 manifest
## ------------------------------------------------------------

manifest = pd.read_csv(
    MAG_MANIFEST,
    sep="\t",
    dtype=str
)


if len(manifest) != EXPECTED:

    raise RuntimeError(
        f"Expected {EXPECTED} rows in "
        f"{MAG_MANIFEST}, found {len(manifest)}."
    )


required_columns = {
    "Genome_ID",
    "Functional_input_FASTA"
}


missing_columns = (
    required_columns
    - set(
        manifest.columns
    )
)


if missing_columns:

    raise RuntimeError(
        "SemiBin2 manifest is missing columns: "
        + ", ".join(
            sorted(
                missing_columns
            )
        )
    )


if manifest["Genome_ID"].duplicated().any():

    duplicates = (
        manifest.loc[
            manifest["Genome_ID"].duplicated(
                keep=False
            ),
            "Genome_ID"
        ]
        .sort_values()
        .tolist()
    )

    raise RuntimeError(
        "Duplicate Genome_ID values in manifest:\n"
        + "\n".join(
            duplicates
        )
    )


manifest_by_genome = (
    manifest
    .set_index(
        "Genome_ID"
    )
)


## ------------------------------------------------------------
## Discover FeGenie proteomes
## ------------------------------------------------------------

faa_files = sorted(
    ORF_DIR.glob(
        "*-proteins.faa"
    )
)


if len(faa_files) != EXPECTED:

    raise RuntimeError(
        f"Expected {EXPECTED} FeGenie proteomes, "
        f"found {len(faa_files)}."
    )


## ------------------------------------------------------------
## Parse FeGenie protein filename
## ------------------------------------------------------------

def genome_from_fegenie_faa(path):

    name = path.name


    suffix = "-proteins.faa"


    if not name.endswith(
        suffix
    ):

        raise RuntimeError(
            f"Unexpected FeGenie protein filename:\n"
            f"{name}"
        )


    staged_name = name[
        :-len(
            suffix
        )
    ]


    ## Example:
    ##
    ## semibin2__flye__barcode10_seqs__SemiBin_0.fa
    ## ->
    ## semibin2__flye__barcode10_seqs__SemiBin_0

    for fasta_suffix in [
        ".fa",
        ".fasta",
        ".fna"
    ]:

        if staged_name.endswith(
            fasta_suffix
        ):

            staged_name = staged_name[
                :-len(
                    fasta_suffix
                )
            ]

            break


    return staged_name


## ------------------------------------------------------------
## Validate exact correspondence
## ------------------------------------------------------------

fegenie_genomes = {
    genome_from_fegenie_faa(
        path
    )
    for path in faa_files
}


expected_genomes = set(
    manifest[
        "Genome_ID"
    ]
)


missing_from_fegenie = (
    expected_genomes
    - fegenie_genomes
)


unexpected_in_fegenie = (
    fegenie_genomes
    - expected_genomes
)


if missing_from_fegenie:

    raise RuntimeError(
        "MAGs missing from FeGenie ORF_calls:\n"
        + "\n".join(
            sorted(
                missing_from_fegenie
            )
        )
    )


if unexpected_in_fegenie:

    raise RuntimeError(
        "Unexpected proteomes in FeGenie ORF_calls:\n"
        + "\n".join(
            sorted(
                unexpected_in_fegenie
            )
        )
    )


## ------------------------------------------------------------
## Rebuild collected-protein directory
## ------------------------------------------------------------

if FAA_OUT.exists():

    shutil.rmtree(
        FAA_OUT
    )


FAA_OUT.mkdir(
    parents=True,
    exist_ok=True
)


METADATA_OUT.mkdir(
    parents=True,
    exist_ok=True
)


## ------------------------------------------------------------
## Copy proteins and build provenance table
## ------------------------------------------------------------

rows = []


for faa in faa_files:

    genome = genome_from_fegenie_faa(
        faa
    )


    row = manifest_by_genome.loc[
        genome
    ]


    functional_input = Path(
        row[
            "Functional_input_FASTA"
        ]
    )


    if not functional_input.is_file():

        raise RuntimeError(
            f"Functional input FASTA for {genome} "
            f"does not exist:\n"
            f"{functional_input}"
        )


    collected = (
        FAA_OUT
        / f"{genome}.faa"
    )


    shutil.copy2(
        faa,
        collected
    )


    rows.append({

        "Genome_ID":
            genome,

        "Functional_input_FASTA":
            str(
                functional_input
            ),

        "FeGenie_ORF_file":
            str(
                faa
            ),

        "Collected_FAA":
            str(
                collected
            )
    })


protein_manifest = pd.DataFrame(
    rows
)


protein_manifest = (
    protein_manifest
    .sort_values(
        "Genome_ID"
    )
    .reset_index(
        drop=True
    )
)


protein_manifest.to_csv(
    PROTEIN_MANIFEST,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Final validation
## ------------------------------------------------------------

collected_files = sorted(
    FAA_OUT.glob(
        "*.faa"
    )
)


if len(collected_files) != EXPECTED:

    raise RuntimeError(
        f"Expected {EXPECTED} collected proteomes, "
        f"found {len(collected_files)}."
    )


if len(protein_manifest) != EXPECTED:

    raise RuntimeError(
        f"Expected {EXPECTED} protein-manifest rows, "
        f"found {len(protein_manifest)}."
    )


## ------------------------------------------------------------
## Completion marker
## ------------------------------------------------------------

with COMPLETE.open(
    "w",
    encoding="utf-8"
) as handle:

    handle.write(
        f"run_finalized="
        f"{datetime.now().astimezone().isoformat()}\n"
    )

    handle.write(
        "status=FeGenie_completed_successfully\n"
    )

    handle.write(
        "note=Original wrapper stopped after FeGenie because "
        "it incorrectly expected FinalSummary.csv\n"
    )

    handle.write(
        f"input_MAGs={EXPECTED}\n"
    )

    handle.write(
        f"predicted_proteomes="
        f"{len(faa_files)}\n"
    )

    handle.write(
        f"collected_proteomes="
        f"{len(collected_files)}\n"
    )

    handle.write(
        f"fegenie_results="
        f"{FEGENIE_OUT}\n"
    )

    handle.write(
        f"gene_summary="
        f"{FEGENIE_OUT / 'FeGenie-geneSummary.csv'}\n"
    )

    handle.write(
        f"cluster_summary="
        f"{FEGENIE_OUT / 'FeGenie-geneSummary-clusters.csv'}\n"
    )

    handle.write(
        f"heatmap_data="
        f"{FEGENIE_OUT / 'FeGenie-heatmap-data.csv'}\n"
    )

    handle.write(
        f"collected_proteins="
        f"{FAA_OUT}\n"
    )

    handle.write(
        f"manifest="
        f"{PROTEIN_MANIFEST}\n"
    )


## ------------------------------------------------------------
## Report
## ------------------------------------------------------------

print()

print(
    "=" * 72
)

print(
    "SemiBin2 FeGenie run finalized"
)

print(
    "=" * 72
)

print(
    f"FeGenie proteomes:      "
    f"{len(faa_files)}"
)

print(
    f"Collected proteomes:    "
    f"{len(collected_files)}"
)

print(
    f"Protein manifest rows:  "
    f"{len(protein_manifest)}"
)

print()

print(
    "FeGenie gene summary:"
)

print(
    FEGENIE_OUT
    / "FeGenie-geneSummary.csv"
)

print()

print(
    "Collected proteins:"
)

print(
    FAA_OUT
)

print()

print(
    "Completion marker:"
)

print(
    COMPLETE
)

print(
    "=" * 72
)
