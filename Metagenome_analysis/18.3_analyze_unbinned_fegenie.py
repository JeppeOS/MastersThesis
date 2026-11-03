#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict


## ============================================================
## PATHS
## ============================================================

ROOT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

FEGENIE_DIR = (
    ROOT
    / "functional_analysis"
    / "unbinned_contigs"
    / "FeGenie_unbinned_full_run_20260918"
    / "fegenie_results"
)

GENE_SUMMARY = (
    FEGENIE_DIR
    / "FeGenie-geneSummary.csv"
)

OUTDIR = (
    ROOT
    / "functional_analysis"
    / "unbinned_contigs"
    / "quick_EET_analysis"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


## ============================================================
## SAMPLE LABELS
## ============================================================

sample_labels = {
    "barcode10": "BRt5_2",
    "barcode11": "ULt4_2",
    "barcode12": "FHt_2",
    "barcode13": "UL_3 (Anoxic)",
    "barcode14": "UL_SAt_2",
    "barcode15": "UL_SAt_3",
    "barcode16": "UL_SA_3",
}


## ============================================================
## HELPERS
## ============================================================

def normalise_assembly(value):

    name = Path(
        value.strip()
    ).name

    for suffix in (
        "-proteins.faa",
        ".fasta",
        ".fna",
        ".faa",
        ".fa"
    ):

        if name.endswith(
            suffix
        ):
            name = name[
                :-len(suffix)
            ]
            break

    return name


def extract_sample(assembly):

    match = re.search(
        r"(barcode\d+)",
        assembly,
        flags=re.IGNORECASE
    )

    if match:
        return match.group(1).lower()

    return ""


## ============================================================
## EET-RELEVANT FEGENIE HMM CLASSIFICATION
##
## Classification is based ONLY on the actual HMM field.
## ============================================================

def classify_hmm(hmm):

    hmm = hmm.strip()


    if re.fullmatch(
        r"Cyc2_repCluster[123]",
        hmm,
        flags=re.IGNORECASE
    ):
        return "Cyc2"


    if re.fullmatch(
        r"MtoA",
        hmm,
        flags=re.IGNORECASE
    ):
        return "MtoA"


    if re.fullmatch(
        r"MtrA",
        hmm,
        flags=re.IGNORECASE
    ):
        return "MtrA"


    if re.fullmatch(
        r"MtrB(?:_TIGR03509)?",
        hmm,
        flags=re.IGNORECASE
    ):
        return "MtrB"


    if re.fullmatch(
        r"MtrC(?:_TIGR03507)?",
        hmm,
        flags=re.IGNORECASE
    ):
        return "MtrC"


    if re.fullmatch(
        r"CymA",
        hmm,
        flags=re.IGNORECASE
    ):
        return "CymA"


    if re.fullmatch(
        r"OmcF",
        hmm,
        flags=re.IGNORECASE
    ):
        return "OmcF"


    if re.fullmatch(
        r"OmcS",
        hmm,
        flags=re.IGNORECASE
    ):
        return "OmcS"


    if re.fullmatch(
        r"OmcZ",
        hmm,
        flags=re.IGNORECASE
    ):
        return "OmcZ"


    if re.fullmatch(
        r"EetA",
        hmm,
        flags=re.IGNORECASE
    ):
        return "EetA"


    if re.fullmatch(
        r"EetB",
        hmm,
        flags=re.IGNORECASE
    ):
        return "EetB"


    if re.fullmatch(
        r"Cyc1",
        hmm,
        flags=re.IGNORECASE
    ):
        return "Cyc1"


    if re.fullmatch(
        r"DFE_0(?:448|449|450|451|461|462|463|464|465)",
        hmm,
        flags=re.IGNORECASE
    ):
        return "DFE system"


    if re.fullmatch(
        r"GACE_184[3-7]",
        hmm,
        flags=re.IGNORECASE
    ):
        return "GACE system"


    if re.fullmatch(
        r"FoxE",
        hmm,
        flags=re.IGNORECASE
    ):
        return "FoxE"


    return ""


## ============================================================
## READ FEGENIE GENE SUMMARY
## ============================================================

if not GENE_SUMMARY.is_file():

    raise SystemExit(
        f"ERROR: FeGenie gene summary not found:\n{GENE_SUMMARY}"
    )


with GENE_SUMMARY.open(
    "r",
    encoding="utf-8",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle
    )

    required = {
        "category",
        "genome/assembly",
        "orf",
        "HMM",
        "bitscore",
        "bitscore_cutoff",
        "clusterID",
        "heme_c_binding_motifs",
        "protein_sequence",
    }

    missing = required.difference(
        reader.fieldnames or []
    )

    if missing:

        raise SystemExit(
            "ERROR: missing FeGenie columns: "
            +
            ", ".join(
                sorted(
                    missing
                )
            )
        )


    rows = []

    extra_field_rows = 0


    for row in reader:

        ## ----------------------------------------------------
        ## FeGenie output can occasionally contain repeated
        ## header rows.
        ## ----------------------------------------------------

        if (
            row["genome/assembly"].strip()
            in {
                "genome/assembly",
                "assembly"
            }
            and
            row["orf"].strip()
            ==
            "orf"
        ):

            continue


        ## ----------------------------------------------------
        ## csv.DictReader puts genuinely extra CSV fields under
        ## key None. Track this as QC.
        ## ----------------------------------------------------

        if (
            None in row
            and
            row[None]
        ):

            extra_field_rows += 1


        assembly = normalise_assembly(
            row["genome/assembly"]
        )

        sample = extract_sample(
            assembly
        )


        cleaned = {
            "sample":
                sample,

            "sample_label":
                sample_labels.get(
                    sample,
                    ""
                ),

            "category":
                row["category"].strip(),

            "assembly":
                assembly,

            "orf":
                row["orf"].strip(),

            "HMM":
                row["HMM"].strip(),

            "bitscore":
                row["bitscore"].strip(),

            "bitscore_cutoff":
                row["bitscore_cutoff"].strip(),

            "clusterID":
                row["clusterID"].strip(),

            "heme_c_binding_motifs":
                row["heme_c_binding_motifs"].strip(),

            "protein_sequence":
                row["protein_sequence"].strip(),
        }


        cleaned["EET_target"] = classify_hmm(
            cleaned["HMM"]
        )


        rows.append(
            cleaned
        )


## ============================================================
## BASIC QC
## ============================================================

print("=" * 90)
print("UNBINNED FEGENIE AUDIT - CORRECTED PARSER")
print("=" * 90)
print()

print(
    f"Accepted FeGenie rows: {len(rows):,}"
)

print(
    f"Rows with unexpected extra CSV fields: {extra_field_rows:,}"
)

print()


## ============================================================
## SAMPLE QC
## ============================================================

sample_counts = Counter(
    row["sample"]
    for row in rows
)


print("FEGENIE HITS BY SAMPLE")
print("-" * 90)


for sample in sorted(
    sample_counts
):

    print(
        f"{sample:12s} "
        f"{sample_labels.get(sample, 'UNKNOWN'):18s} "
        f"{sample_counts[sample]:5d}"
    )


unknown_samples = [
    row
    for row in rows
    if not row["sample"]
]


if unknown_samples:

    print()
    print(
        "WARNING: FeGenie rows without recoverable barcode:",
        len(
            unknown_samples
        )
    )


## ============================================================
## ALL CATEGORY COUNTS
## ============================================================

category_counts = Counter(
    row["category"]
    for row in rows
)


print()
print("=" * 90)
print("FEGENIE CATEGORY COUNTS")
print("=" * 90)


for category, count in sorted(
    category_counts.items(),
    key=lambda x: (
        -x[1],
        x[0]
    )
):

    print(
        f"{category:40s} {count:5d}"
    )


## ============================================================
## TARGETED EET HITS
## ============================================================

target_rows = [
    row
    for row in rows
    if row["EET_target"]
]


print()
print("=" * 90)
print("EET-RELATED FEGENIE HITS")
print("=" * 90)
print()

print(
    f"Targeted hit rows: {len(target_rows):,}"
)


## ============================================================
## COUNTS BY SAMPLE + TARGET
## ============================================================

target_counts = Counter(
    (
        row["sample"],
        row["EET_target"]
    )
    for row in target_rows
)


print()
print(
    f"{'Sample':12s} "
    f"{'Experimental sample':20s} "
    f"{'Target':16s} "
    f"{'Hits':>5s}"
)

print(
    "-" * 60
)


for (
    sample,
    target
), count in sorted(
    target_counts.items()
):

    print(
        f"{sample:12s} "
        f"{sample_labels.get(sample, ''):20s} "
        f"{target:16s} "
        f"{count:5d}"
    )


## ============================================================
## DETAILED TARGET TABLE
## ============================================================

print()
print("=" * 90)
print("DETAILED TARGETED HITS")
print("=" * 90)
print()


for row in sorted(
    target_rows,
    key=lambda r: (
        r["sample"],
        r["EET_target"],
        r["assembly"],
        r["orf"]
    )
):

    print(
        "\t".join(
            [
                row["sample"],
                row["sample_label"],
                row["category"],
                row["assembly"],
                row["orf"],
                row["HMM"],
                row["bitscore"],
                row["bitscore_cutoff"],
                row["heme_c_binding_motifs"],
            ]
        )
    )


## ============================================================
## OUTPUT TABLES
## ============================================================

fields = [
    "sample",
    "sample_label",
    "category",
    "assembly",
    "orf",
    "HMM",
    "bitscore",
    "bitscore_cutoff",
    "clusterID",
    "heme_c_binding_motifs",
    "EET_target",
    "protein_sequence",
]


all_out = (
    OUTDIR
    / "unbinned_all_FeGenie_hits.tsv"
)


with all_out.open(
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
        rows
    )


target_out = (
    OUTDIR
    / "unbinned_EET_FeGenie_hits.tsv"
)


with target_out.open(
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
        target_rows
    )


counts_out = (
    OUTDIR
    / "unbinned_EET_FeGenie_counts_by_sample.tsv"
)


with counts_out.open(
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
            "sample",
            "sample_label",
            "target",
            "n_hits"
        ]
    )


    for (
        sample,
        target
    ), count in sorted(
        target_counts.items()
    ):

        writer.writerow(
            [
                sample,
                sample_labels.get(
                    sample,
                    ""
                ),
                target,
                count
            ]
        )


## ============================================================
## FASTA OF TARGETED PROTEINS
##
## Deduplicated by assembly + ORF.
## ============================================================

fasta_out = (
    OUTDIR
    / "unbinned_EET_FeGenie_proteins.faa"
)


seen = set()


with fasta_out.open(
    "w",
    encoding="utf-8"
) as handle:

    for row in target_rows:

        key = (
            row["assembly"],
            row["orf"]
        )

        if key in seen:
            continue

        seen.add(
            key
        )


        sequence = (
            row["protein_sequence"]
            .rstrip(
                "*"
            )
        )


        handle.write(
            f">{row['sample']}|"
            f"{row['assembly']}|"
            f"{row['orf']}|"
            f"{row['EET_target']}|"
            f"{row['HMM']}\n"
        )


        for start in range(
            0,
            len(sequence),
            80
        ):

            handle.write(
                sequence[
                    start:start + 80
                ]
                +
                "\n"
            )


print()
print("=" * 90)
print("OUTPUTS")
print("=" * 90)

for path in (
    all_out,
    target_out,
    counts_out,
    fasta_out
):

    print(
        path
    )

print()
