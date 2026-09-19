#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import sys

import pandas as pd


## ================================================================== ##
## Paths
## ================================================================== ##

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE.parent

ORF_DIR = (
    WORKFLOW
    / "FeGenie_conda_full_run_20260803"
    / "fegenie_results"
    / "ORF_calls"
)

MODULE_PROTEINS = (
    WORKFLOW
    / "09_module_occurrence"
    / "protein_module_membership.tsv"
)


## ================================================================== ##
## Outputs
## ================================================================== ##

OUT_ALL_COORDS = (
    HERE
    / "all_prodigal_gene_coordinates.tsv"
)

OUT_MODULE_COORDS = (
    HERE
    / "module_protein_coordinates.tsv"
)

OUT_GENOME_QC = (
    HERE
    / "genome_coordinate_qc.tsv"
)

OUT_MAPPING_QC = (
    HERE
    / "coordinate_mapping_qc.tsv"
)

OUT_UNMATCHED = (
    HERE
    / "unmatched_module_proteins.tsv"
)


## ================================================================== ##
## Expected values from upstream locked stages
## ================================================================== ##

EXPECTED_GENOMES = 631
EXPECTED_MODULE_PROTEINS = 10537


## ================================================================== ##
## Prodigal FASTA header
##
## Example:
##
## >GCA_002134785_000000000094_46
##     # 44022
##     # 44987
##     # -1
##     # ID=94_46;partial=00;start_type=ATG;...
## ================================================================== ##

HEADER_RE = re.compile(
    r"^>(\S+)"
    r"\s*#\s*(\d+)"
    r"\s*#\s*(\d+)"
    r"\s*#\s*(-?1)"
    r"\s*#\s*(.*)$"
)


## ================================================================== ##
## Utilities
## ================================================================== ##

def fail(message):

    print(
        f"\nERROR: {message}",
        file=sys.stderr
    )

    sys.exit(1)


def parse_attributes(text):

    attributes = {}

    for field in text.split(";"):

        field = field.strip()

        if not field:
            continue

        if "=" not in field:
            continue

        key, value = field.split(
            "=",
            1
        )

        attributes[
            key.strip()
        ] = value.strip()

    return attributes


def parse_header(header, genome):

    match = HEADER_RE.match(
        header.rstrip()
    )

    if match is None:

        fail(
            "Could not parse Prodigal FASTA header:\n"
            f"{header}"
        )

    protein_id = match.group(1)

    start = int(
        match.group(2)
    )

    end = int(
        match.group(3)
    )

    strand_numeric = int(
        match.group(4)
    )

    attribute_string = (
        match.group(5)
    )

    attributes = parse_attributes(
        attribute_string
    )


    ## -------------------------------------------------------------- ##
    ## Validate coordinates
    ## -------------------------------------------------------------- ##

    if start < 1:

        fail(
            f"Invalid start coordinate for "
            f"{protein_id}: {start}"
        )

    if end < start:

        fail(
            f"End < start for "
            f"{protein_id}: "
            f"{start}-{end}"
        )


    ## -------------------------------------------------------------- ##
    ## Parse Prodigal's internal ID
    ##
    ## ID=94_46 means:
    ##
    ##     sequence ordinal = 94
    ##     gene ordinal     = 46
    ## -------------------------------------------------------------- ##

    prodigal_id = attributes.get(
        "ID",
        ""
    )

    if not prodigal_id:

        fail(
            f"No Prodigal ID= field for "
            f"{protein_id}"
        )

    if "_" not in prodigal_id:

        fail(
            f"Unexpected Prodigal ID for "
            f"{protein_id}: "
            f"{prodigal_id}"
        )

    seqnum, gene_ordinal = (
        prodigal_id.rsplit(
            "_",
            1
        )
    )

    if (
        not seqnum.isdigit()
        or
        not gene_ordinal.isdigit()
    ):

        fail(
            f"Unexpected Prodigal ID for "
            f"{protein_id}: "
            f"{prodigal_id}"
        )


    ## -------------------------------------------------------------- ##
    ## Recover original contig ID
    ##
    ## Prodigal protein identifiers are:
    ##
    ##     original_contig_id + "_" + gene_ordinal
    ##
    ## Example:
    ##
    ##     protein:
    ##       GCA_002134785_000000000094_46
    ##
    ##     gene ordinal:
    ##       46
    ##
    ##     contig:
    ##       GCA_002134785_000000000094
    ##
    ## We validate this instead of blindly stripping the suffix.
    ## -------------------------------------------------------------- ##

    expected_suffix = (
        "_"
        + gene_ordinal
    )

    if not protein_id.endswith(
        expected_suffix
    ):

        fail(
            "Protein ID does not end in its "
            "Prodigal gene ordinal:\n"
            f"protein_id = {protein_id}\n"
            f"Prodigal ID = {prodigal_id}"
        )

    contig = protein_id[
        :-len(expected_suffix)
    ]


    ## -------------------------------------------------------------- ##
    ## Convert strand
    ## -------------------------------------------------------------- ##

    if strand_numeric == 1:
        strand = "+"
    elif strand_numeric == -1:
        strand = "-"
    else:
        fail(
            f"Unexpected strand for "
            f"{protein_id}: "
            f"{strand_numeric}"
        )


    return {
        "genome": genome,
        "protein_id": protein_id,
        "contig": contig,

        "start": start,
        "end": end,

        "gene_length_nt": (
            end
            - start
            + 1
        ),

        "strand": strand,

        "prodigal_strand": (
            strand_numeric
        ),

        "prodigal_id": (
            prodigal_id
        ),

        "prodigal_seqnum": (
            int(seqnum)
        ),

        "prodigal_gene_ordinal": (
            int(gene_ordinal)
        ),

        "partial": attributes.get(
            "partial",
            ""
        ),

        "start_type": attributes.get(
            "start_type",
            ""
        ),

        "rbs_motif": attributes.get(
            "rbs_motif",
            ""
        ),

        "rbs_spacer": attributes.get(
            "rbs_spacer",
            ""
        ),

        "gc_cont": attributes.get(
            "gc_cont",
            ""
        ),
    }


## ================================================================== ##
## 1. Read module proteins
## ================================================================== ##

print("=" * 80)
print("BUILD GENOMIC COORDINATE TABLES")
print("=" * 80)


if not MODULE_PROTEINS.exists():

    fail(
        f"Cannot find:\n"
        f"{MODULE_PROTEINS}"
    )


module_df = pd.read_csv(
    MODULE_PROTEINS,
    sep="\t",
    dtype=str
)


required = {
    "genome",
    "protein_id",
    "cluster",
    "module",
}


missing = (
    required
    - set(
        module_df.columns
    )
)


if missing:

    fail(
        "Missing required columns in "
        "protein_module_membership.tsv:\n"
        + ", ".join(
            sorted(missing)
        )
    )


## genome + protein_id must remain authoritative ##

duplicates = module_df.duplicated(
    [
        "genome",
        "protein_id",
    ],
    keep=False
)


if duplicates.any():

    fail(
        "Duplicate genome + protein_id keys "
        "in protein_module_membership.tsv:\n"
        +
        module_df.loc[
            duplicates,
            [
                "genome",
                "protein_id",
                "cluster",
                "module",
            ]
        ]
        .head(20)
        .to_string(
            index=False
        )
    )


n_module_proteins = len(
    module_df
)


print()
print(
    f"Module proteins expected: "
    f"{n_module_proteins:,}"
)


if (
    n_module_proteins
    != EXPECTED_MODULE_PROTEINS
):

    fail(
        f"Expected "
        f"{EXPECTED_MODULE_PROTEINS:,} "
        f"module proteins but found "
        f"{n_module_proteins:,}."
    )


target_keys = set(
    zip(
        module_df["genome"],
        module_df["protein_id"],
    )
)


## ================================================================== ##
## 2. Locate authoritative FeGenie protein FASTAs
##
## Important:
##
## *.fa-proteins.faa
##
## does NOT match the auxiliary:
##
## *.pdb
## *.pot
## *.pto
## *.ptf
## *.pjs
## ================================================================== ##

faa_files = sorted(
    ORF_DIR.glob(
        "*.fa-proteins.faa"
    )
)


print()
print(
    f"FeGenie protein FASTAs: "
    f"{len(faa_files):,}"
)


if len(faa_files) != EXPECTED_GENOMES:

    fail(
        f"Expected {EXPECTED_GENOMES} "
        f"FeGenie protein FASTAs but found "
        f"{len(faa_files)}."
    )


## Validate genome names from filenames ##

genome_to_faa = {}


suffix = (
    ".fa-proteins.faa"
)


for path in faa_files:

    name = path.name

    if not name.endswith(
        suffix
    ):

        fail(
            f"Unexpected FASTA filename:\n"
            f"{name}"
        )

    genome = name[
        :-len(suffix)
    ]

    if genome in genome_to_faa:

        fail(
            f"More than one authoritative "
            f"protein FASTA for genome "
            f"{genome}"
        )

    genome_to_faa[
        genome
    ] = path


module_genomes = set(
    module_df["genome"]
)


missing_genome_fastas = sorted(
    module_genomes
    - set(
        genome_to_faa
    )
)


if missing_genome_fastas:

    fail(
        "Module-associated genomes lacking "
        "a FeGenie protein FASTA:\n"
        + "\n".join(
            missing_genome_fastas[:20]
        )
    )


## ================================================================== ##
## 3. Parse every predicted ORF
##
## We stream the complete coordinate table to disk rather than
## keeping ~all proteins from 631 genomes in memory.
## ================================================================== ##

all_fields = [
    "genome",
    "protein_id",
    "contig",
    "start",
    "end",
    "gene_length_nt",
    "strand",
    "prodigal_strand",
    "prodigal_id",
    "prodigal_seqnum",
    "prodigal_gene_ordinal",
    "partial",
    "start_type",
    "rbs_motif",
    "rbs_spacer",
    "gc_cont",
    "source_protein_fasta",
]


temp_all_coords = Path(
    str(
        OUT_ALL_COORDS
    )
    + ".tmp"
)


module_coordinate_records = []

matched_module_keys = set()

genome_qc_records = []

total_predicted_proteins = 0


with open(
    temp_all_coords,
    "w",
    newline=""
) as output_handle:

    writer = csv.DictWriter(
        output_handle,
        fieldnames=all_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()


    for file_number, faa in enumerate(
        faa_files,
        start=1
    ):

        genome = faa.name[
            :-len(suffix)
        ]

        seen_protein_ids = set()

        contigs = set()

        genome_protein_count = 0

        genome_module_expected = int(
            (
                module_df[
                    "genome"
                ]
                == genome
            )
            .sum()
        )

        genome_module_mapped = 0


        with open(
            faa
        ) as handle:

            for line in handle:

                if not line.startswith(
                    ">"
                ):
                    continue


                record = parse_header(
                    line,
                    genome
                )


                protein_id = record[
                    "protein_id"
                ]


                ## Within-genome protein IDs must be unique ##

                if (
                    protein_id
                    in seen_protein_ids
                ):

                    fail(
                        "Duplicate protein ID "
                        f"inside {faa.name}:\n"
                        f"{protein_id}"
                    )


                seen_protein_ids.add(
                    protein_id
                )

                contigs.add(
                    record[
                        "contig"
                    ]
                )

                genome_protein_count += 1

                total_predicted_proteins += 1


                record[
                    "source_protein_fasta"
                ] = str(
                    faa
                )


                writer.writerow(
                    record
                )


                key = (
                    genome,
                    protein_id,
                )


                if key in target_keys:

                    if (
                        key
                        in matched_module_keys
                    ):

                        fail(
                            "Module protein received "
                            "more than one coordinate "
                            "match:\n"
                            f"{genome}\t"
                            f"{protein_id}"
                        )


                    matched_module_keys.add(
                        key
                    )

                    genome_module_mapped += 1

                    module_coordinate_records.append(
                        record.copy()
                    )


        genome_qc_records.append(
            {
                "genome": genome,

                "n_predicted_proteins":
                    genome_protein_count,

                "n_contigs_with_predicted_proteins":
                    len(
                        contigs
                    ),

                "n_module_proteins_expected":
                    genome_module_expected,

                "n_module_proteins_mapped":
                    genome_module_mapped,

                "module_coordinate_mapping_complete":
                    int(
                        genome_module_expected
                        ==
                        genome_module_mapped
                    ),

                "source_protein_fasta":
                    str(
                        faa
                    ),
            }
        )


        if (
            file_number % 50 == 0
            or
            file_number
            == len(faa_files)
        ):

            print(
                f"  Parsed "
                f"{file_number:>3}/"
                f"{len(faa_files)} genomes"
            )


## ================================================================== ##
## 4. Check module-protein coordinate recovery
## ================================================================== ##

unmatched_keys = (
    target_keys
    - matched_module_keys
)


print()
print("Coordinate mapping")
print(
    f"  Total predicted proteins parsed: "
    f"{total_predicted_proteins:,}"
)
print(
    f"  Module proteins expected:        "
    f"{len(target_keys):,}"
)
print(
    f"  Module proteins matched:         "
    f"{len(matched_module_keys):,}"
)
print(
    f"  Module proteins unmatched:       "
    f"{len(unmatched_keys):,}"
)


## Write unmatched file even if empty ##

if unmatched_keys:

    unmatched_df = pd.DataFrame(
        sorted(
            unmatched_keys
        ),
        columns=[
            "genome",
            "protein_id",
        ]
    )

else:

    unmatched_df = pd.DataFrame(
        columns=[
            "genome",
            "protein_id",
        ]
    )


unmatched_df.to_csv(
    OUT_UNMATCHED,
    sep="\t",
    index=False
)


if unmatched_keys:

    fail(
        "Not every module protein could "
        "be mapped to its authoritative "
        "FeGenie Prodigal FASTA.\n\n"
        "See:\n"
        f"{OUT_UNMATCHED}"
    )


## ================================================================== ##
## 5. Coordinate table passed strict QC
##
## Promote temporary complete coordinate table to authoritative file.
## ================================================================== ##

temp_all_coords.replace(
    OUT_ALL_COORDS
)


## ================================================================== ##
## 6. Build module_protein_coordinates.tsv
## ================================================================== ##

coords_df = pd.DataFrame(
    module_coordinate_records
)


## Coordinate subset must itself be unique ##

if coords_df.duplicated(
    [
        "genome",
        "protein_id",
    ]
).any():

    fail(
        "Duplicate coordinate rows found "
        "for module proteins."
    )


## Remove source path from right table temporarily only if needed.
## Keep it in final output because it is useful for traceability. ##

module_coordinates = (
    module_df
    .merge(
        coords_df,
        on=[
            "genome",
            "protein_id",
        ],
        how="left",
        validate="one_to_one"
    )
)


if len(
    module_coordinates
) != EXPECTED_MODULE_PROTEINS:

    fail(
        "Joined module coordinate table "
        "has unexpected number of rows."
    )


## Put structural identifiers first ##

coordinate_columns = [
    "genome",
    "protein_id",
    "cluster",
    "module",
    "contig",
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
]


remaining_columns = [
    col
    for col in module_coordinates.columns
    if col not in coordinate_columns
]


module_coordinates = (
    module_coordinates[
        coordinate_columns
        + remaining_columns
    ]
    .sort_values(
        [
            "module",
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


module_coordinates.to_csv(
    OUT_MODULE_COORDS,
    sep="\t",
    index=False
)


## ================================================================== ##
## 7. Genome-level QC
## ================================================================== ##

genome_qc = pd.DataFrame(
    genome_qc_records
)


genome_qc = (
    genome_qc
    .sort_values(
        "genome"
    )
    .reset_index(
        drop=True
    )
)


genome_qc.to_csv(
    OUT_GENOME_QC,
    sep="\t",
    index=False
)


if (
    genome_qc[
        "module_coordinate_mapping_complete"
    ]
    != 1
).any():

    fail(
        "One or more genomes failed "
        "per-genome coordinate QC."
    )


## ================================================================== ##
## 8. Global QC
## ================================================================== ##

qc = pd.DataFrame(
    [
        [
            "fegenie_protein_fastas",
            len(
                faa_files
            ),
        ],
        [
            "genomes_expected",
            EXPECTED_GENOMES,
        ],
        [
            "total_predicted_proteins_parsed",
            total_predicted_proteins,
        ],
        [
            "module_proteins_expected",
            len(
                target_keys
            ),
        ],
        [
            "module_proteins_uniquely_mapped",
            len(
                matched_module_keys
            ),
        ],
        [
            "module_proteins_unmatched",
            len(
                unmatched_keys
            ),
        ],
        [
            "genomes_with_complete_module_coordinate_mapping",
            int(
                genome_qc[
                    "module_coordinate_mapping_complete"
                ].sum()
            ),
        ],
    ],
    columns=[
        "metric",
        "value",
    ]
)


qc.to_csv(
    OUT_MAPPING_QC,
    sep="\t",
    index=False
)


## ================================================================== ##
## Final report
## ================================================================== ##

print()
print("Strict coordinate QC")
print(
    f"  FeGenie FASTAs:           "
    f"{len(faa_files):,}"
)
print(
    f"  Genomes passing mapping:  "
    f"{int(genome_qc['module_coordinate_mapping_complete'].sum()):,}"
)
print(
    f"  Unique module matches:    "
    f"{len(matched_module_keys):,}/"
    f"{EXPECTED_MODULE_PROTEINS:,}"
)


print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)

print(
    f"All Prodigal coordinates:   "
    f"{OUT_ALL_COORDS}"
)

print(
    f"Module protein coordinates: "
    f"{OUT_MODULE_COORDS}"
)

print(
    f"Genome QC:                  "
    f"{OUT_GENOME_QC}"
)

print(
    f"Mapping QC:                 "
    f"{OUT_MAPPING_QC}"
)

print(
    f"Unmatched proteins:         "
    f"{OUT_UNMATCHED}"
)
