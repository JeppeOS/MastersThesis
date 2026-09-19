#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re


PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)

WF = PROJECT / "genome_analysis_workflow"

C63 = (
    WF
    / "08_fegenie_composition"
    / "Cluster_00063_MCA0421"
)

GENOMES_FILE = (
    C63
    / "cluster_00063_genomes.tsv"
)

TAXONOMY_FILE = (
    WF
    / "16_visualization"
    / "globdb_r226_taxonomy.tsv"
)

GFFROOT = (
    PROJECT
    / "reference_data"
    / "globdb_r226_gff_cog"
)

OUTDIR = (
    C63
    / "sMMO_screen"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


## ================================================================ ##
## Helpers
## ================================================================ ##

def find_column(df, candidates):

    lookup = {
        x.lower(): x
        for x in df.columns
    }

    for candidate in candidates:

        if candidate.lower() in lookup:
            return lookup[candidate.lower()]

    raise RuntimeError(
        f"Could not identify genome column. "
        f"Columns = {list(df.columns)}"
    )


## ================================================================ ##
## Cluster_00063 genomes
## ================================================================ ##

c63 = pd.read_csv(
    GENOMES_FILE,
    sep="\t",
    dtype=str
)

genome_col = find_column(
    c63,
    [
        "genome_id",
        "genome",
        "assembly",
        "assembly_id",
    ]
)

c63_genomes = set(
    c63[genome_col]
    .dropna()
    .astype(str)
)


## ================================================================ ##
## Taxonomy
## ================================================================ ##

tax = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t",
    header=None,
    names=[
        "genome_id",
        "taxonomy"
    ],
    usecols=[
        0,
        1
    ],
    dtype=str
)


tax = tax[
    tax["genome_id"].isin(
        c63_genomes
    )
].copy()


methylobacter = tax[
    tax["taxonomy"]
    .str.contains(
        r"g__Methylobacter",
        case=False,
        regex=True,
        na=False
    )
].copy()


methylobacter_genomes = sorted(
    methylobacter[
        "genome_id"
    ].unique()
)


print()
print("Cluster_00063 Methylobacter genomes:")
print()

for genome in methylobacter_genomes:
    print(genome)

print()
print(
    f"n = {len(methylobacter_genomes)}"
)


## ================================================================ ##
## Locate GFF files
## ================================================================ ##

all_gffs = list(
    GFFROOT.rglob("*.gff*")
)


## ================================================================ ##
## sMMO annotation patterns
## ================================================================ ##

patterns = {

    "mmoX":
        re.compile(
            r"\bmmoX\b|"
            r"methane monooxygenase.*hydroxylase.*alpha|"
            r"soluble methane monooxygenase.*alpha",
            re.I
        ),

    "mmoY":
        re.compile(
            r"\bmmoY\b|"
            r"methane monooxygenase.*hydroxylase.*beta|"
            r"soluble methane monooxygenase.*beta",
            re.I
        ),

    "mmoZ":
        re.compile(
            r"\bmmoZ\b|"
            r"methane monooxygenase.*hydroxylase.*gamma|"
            r"soluble methane monooxygenase.*gamma",
            re.I
        ),

    "mmoB":
        re.compile(
            r"\bmmoB\b|"
            r"methane monooxygenase.*regulatory protein",
            re.I
        ),

    "mmoC":
        re.compile(
            r"\bmmoC\b|"
            r"methane monooxygenase.*reductase",
            re.I
        ),

    "mmoD":
        re.compile(
            r"\bmmoD\b",
            re.I
        ),
}


## ================================================================ ##
## Search
## ================================================================ ##

rows = []


for genome in methylobacter_genomes:

    matching_files = [
        path
        for path in all_gffs
        if genome in path.name
    ]


    if not matching_files:

        ## Some GFF names may contain genome ID elsewhere in path.
        matching_files = [
            path
            for path in all_gffs
            if genome in str(path)
        ]


    if not matching_files:

        print(
            f"WARNING: no GFF found for {genome}"
        )

        rows.append(
            {
                "genome_id": genome,
                "gene": "NO_GFF_FOUND",
                "gff_file": "",
                "matching_line": "",
            }
        )

        continue


    for gff in matching_files:

        with open(
            gff,
            errors="replace"
        ) as fh:

            for line in fh:

                if line.startswith("#"):
                    continue

                for gene, pattern in patterns.items():

                    if pattern.search(line):

                        rows.append(
                            {
                                "genome_id":
                                    genome,

                                "gene":
                                    gene,

                                "gff_file":
                                    str(gff),

                                "matching_line":
                                    line.rstrip(),
                            }
                        )


hits = pd.DataFrame(
    rows,
    columns=[
        "genome_id",
        "gene",
        "gff_file",
        "matching_line",
    ]
)


HITS_FILE = (
    OUTDIR
    / "methylobacter_cluster00063_smmo_annotation_hits.tsv"
)


hits.to_csv(
    HITS_FILE,
    sep="\t",
    index=False
)


## ================================================================ ##
## Genome-level summary
## ================================================================ ##

summary_rows = []


for genome in methylobacter_genomes:

    genes = set(
        hits.loc[
            (
                hits["genome_id"] == genome
            )
            &
            (
                hits["gene"] != "NO_GFF_FOUND"
            ),
            "gene"
        ]
    )


    summary_rows.append(
        {
            "genome_id":
                genome,

            "mmoX":
                "mmoX" in genes,

            "mmoY":
                "mmoY" in genes,

            "mmoZ":
                "mmoZ" in genes,

            "mmoB":
                "mmoB" in genes,

            "mmoC":
                "mmoC" in genes,

            "mmoD":
                "mmoD" in genes,

            "n_smmo_components":
                len(genes),

            "smmo_annotation_status":
                (
                    "strong_candidate"
                    if {
                        "mmoX",
                        "mmoY",
                        "mmoZ"
                    }.issubset(genes)

                    else
                    (
                        "partial_candidate"
                        if len(genes) > 0
                        else "no_annotation_hit"
                    )
                )
        }
    )


summary = pd.DataFrame(
    summary_rows
)


SUMMARY_FILE = (
    OUTDIR
    / "methylobacter_cluster00063_smmo_summary.tsv"
)


summary.to_csv(
    SUMMARY_FILE,
    sep="\t",
    index=False
)


print()
print(
    "======================================================"
)

print(
    "sMMO annotation screen"
)

print(
    "======================================================"
)

print()

print(
    summary.to_string(
        index=False
    )
)

print()

print("Written:")
print(SUMMARY_FILE)
print(HITS_FILE)
