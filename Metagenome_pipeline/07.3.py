#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import re

PROJECT = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
    / "metagenome"
)

INPUT = (
    PROJECT
    / "abundance_comparison"
    / "sylph_sequence_abundance_merged.tsv"
)

OUTPUT = (
    PROJECT
    / "abundance_comparison"
    / "Sylph_sequence_abundance_by_category.tsv"
)

EXCLUDE_SAMPLES = {
    "barcode09_seqs"
}


with INPUT.open() as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    sample_columns = {}

    for column in reader.fieldnames:

        if column == "clade_name":
            continue

        match = re.search(
            r"(barcode\d+_seqs)",
            column
        )

        if match:
            sample = match.group(1)

            if sample not in EXCLUDE_SAMPLES:
                sample_columns[column] = sample


    abundance = defaultdict(float)
    classified_total = defaultdict(float)


    for row in reader:

        ranks = row["clade_name"].split("|")

        ## Use genus-level rows only to avoid counting the
        ## same taxon at domain/phylum/class/order/etc.
        if not ranks[-1].startswith("g__"):
            continue

        genus = ranks[-1][3:]

        if "o__Methylococcales" in ranks:

            category = (
                genus
                if genus
                else "Methylococcales_unclassified"
            )

        else:

            category = "Other"


        for column, sample in sample_columns.items():

            value = row[column].strip()

            x = (
                float(value)
                if value
                else 0.0
            )

            abundance[
                (sample, category)
            ] += x

            classified_total[
                sample
            ] += x


## Add the unclassified read fraction to Other.
##
## Sylph Sequence_abundance sums to the percentage of
## classified reads, so 100 - classified represents sequence
## not taxonomically assigned by Sylph.
for sample in sorted(
    set(sample_columns.values())
):

    unknown = max(
        0.0,
        100.0 - classified_total[sample]
    )

    abundance[
        (sample, "Other")
    ] += unknown


rows = []

for (
    sample,
    category
), value in sorted(
    abundance.items()
):

    if value <= 0:
        continue

    rows.append({
        "Sample": sample,
        "Category": category,
        "Relative_sequence_abundance_pct": value
    })


with OUTPUT.open(
    "w",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "Sample",
            "Category",
            "Relative_sequence_abundance_pct"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


print(f"Written: {OUTPUT}")
