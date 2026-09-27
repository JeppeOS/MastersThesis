#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import gzip
import shutil


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

MASTER = (
    PROJECT
    / "semibin2_analysis"
    / "MAG_master_table_SemiBin2.tsv"
)

DEPTH_DIR = (
    PROJECT
    / "semibin2_analysis"
    / "depth"
)

OUTROOT = (
    PROJECT
    / "functional_analysis"
    / "unbinned_contigs"
)

INPUT_ROOT = (
    OUTROOT
    / "input"
)

MANIFEST_OUT = (
    OUTROOT
    / "unbinned_contig_manifest.tsv"
)

SUMMARY_OUT = (
    OUTROOT
    / "unbinned_contig_summary.tsv"
)


## ------------------------------------------------------------
## Input MAG table
## ------------------------------------------------------------

master = pd.read_csv(
    MASTER,
    sep="\t"
)


if len(master) != 187:

    raise RuntimeError(
        f"Expected 187 SemiBin2 MAGs, "
        f"found {len(master)}."
    )


samples = sorted(
    master["Sample"].unique()
)


if len(samples) != 7:

    raise RuntimeError(
        f"Expected 7 samples, "
        f"found {len(samples)}."
    )


## ------------------------------------------------------------
## Recreate output input directory
## ------------------------------------------------------------

if INPUT_ROOT.exists():

    shutil.rmtree(
        INPUT_ROOT
    )


INPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)


## ------------------------------------------------------------
## Open plain or gzip FASTA
## ------------------------------------------------------------

def open_fasta(path):

    path = Path(path)

    with open(
        path,
        "rb"
    ) as handle:

        magic = handle.read(2)


    if magic == b"\x1f\x8b":

        return gzip.open(
            path,
            "rt"
        )


    return open(
        path,
        "rt"
    )


## ------------------------------------------------------------
## FASTA iterator
## ------------------------------------------------------------

def fasta_records(path):

    name = None
    sequence = []


    with open_fasta(
        path
    ) as handle:

        for raw in handle:

            line = raw.rstrip()


            if line.startswith(">"):

                if name is not None:

                    yield (
                        name,
                        "".join(sequence)
                    )


                name = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

                sequence = []


            else:

                sequence.append(
                    line.strip()
                )


        if name is not None:

            yield (
                name,
                "".join(sequence)
            )


## ------------------------------------------------------------
## FASTA writer
## ------------------------------------------------------------

def write_record(
    handle,
    name,
    sequence
):

    handle.write(
        f">{name}\n"
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
            + "\n"
        )


## ------------------------------------------------------------
## Outputs
## ------------------------------------------------------------

manifest_rows = []

summary_rows = []


## ------------------------------------------------------------
## Process each barcode
## ------------------------------------------------------------

for sample in samples:

    print()
    print(
        "=" * 70
    )

    print(
        f"Processing {sample}"
    )


    sample_master = master[
        master["Sample"]
        == sample
    ].copy()


    ## --------------------------------------------------------
    ## Find every contig assigned to a SemiBin2 MAG
    ## --------------------------------------------------------

    assigned = {}


    for _, mag in sample_master.iterrows():

        genome_id = str(
            mag["Genome_ID"]
        )

        fasta = Path(
            str(
                mag["Original_path"]
            )
        )


        if not fasta.exists():

            raise RuntimeError(
                f"Missing MAG FASTA:\n"
                f"{fasta}"
            )


        for contig, sequence in fasta_records(
            fasta
        ):

            if contig in assigned:

                raise RuntimeError(
                    f"{sample}: contig "
                    f"{contig} occurs in more "
                    f"than one SemiBin2 MAG:\n"
                    f"{assigned[contig]}\n"
                    f"{genome_id}"
                )


            assigned[
                contig
            ] = genome_id


    ## --------------------------------------------------------
    ## Original Flye assembly
    ## --------------------------------------------------------

    assembly = (
        PROJECT
        / "flye_results"
        / sample
        / "assembly.fasta"
    )


    if not assembly.exists():

        raise RuntimeError(
            f"Missing Flye assembly:\n"
            f"{assembly}"
        )


    ## --------------------------------------------------------
    ## Existing JGI depth table
    ## --------------------------------------------------------

    depth_file = (
        DEPTH_DIR
        / f"{sample}_depth.txt"
    )


    if not depth_file.exists():

        raise RuntimeError(
            f"Missing depth table:\n"
            f"{depth_file}"
        )


    depth = pd.read_csv(
        depth_file,
        sep="\t"
    )


    required = {
        "contigName",
        "contigLen",
        "totalAvgDepth"
    }


    missing = (
        required
        - set(
            depth.columns
        )
    )


    if missing:

        raise RuntimeError(
            f"{depth_file} missing columns: "
            + ", ".join(
                sorted(missing)
            )
        )


    depth = depth.set_index(
        "contigName"
    )


    ## --------------------------------------------------------
    ## FeGenie input directory
    ## --------------------------------------------------------

    sample_input = (
        INPUT_ROOT
        / sample
    )


    sample_input.mkdir(
        parents=True,
        exist_ok=True
    )


    output_fasta = (
        sample_input
        / f"{sample}__unbinned_ge1000.fa"
    )


    ## --------------------------------------------------------
    ## Counters
    ## --------------------------------------------------------

    total_n = 0
    total_bp = 0

    binned_n = 0
    binned_bp = 0

    unbinned_all_n = 0
    unbinned_all_bp = 0

    unbinned_1k_n = 0
    unbinned_1k_bp = 0

    unbinned_2500_n = 0
    unbinned_2500_bp = 0

    seen_assembly_contigs = set()


    ## --------------------------------------------------------
    ## Extract
    ## --------------------------------------------------------

    with open(
        output_fasta,
        "w"
    ) as out:

        for contig, sequence in fasta_records(
            assembly
        ):

            length = len(
                sequence
            )


            seen_assembly_contigs.add(
                contig
            )


            total_n += 1
            total_bp += length


            ## ------------------------------------------------
            ## Assigned to SemiBin2 MAG
            ## ------------------------------------------------

            if contig in assigned:

                binned_n += 1
                binned_bp += length

                continue


            ## ------------------------------------------------
            ## Unbinned
            ## ------------------------------------------------

            unbinned_all_n += 1
            unbinned_all_bp += length


            if contig not in depth.index:

                raise RuntimeError(
                    f"{sample}: {contig} missing "
                    f"from depth table."
                )


            mean_depth = float(
                depth.loc[
                    contig,
                    "totalAvgDepth"
                ]
            )


            manifest_rows.append({

                "Sample":
                    sample,

                "Contig":
                    contig,

                "Length_bp":
                    length,

                "Mean_coverage_x":
                    mean_depth,

                "Binning_eligible_ge2500":
                    int(
                        length >= 2500
                    ),

                "FeGenie_screen_ge1000":
                    int(
                        length >= 1000
                    )
            })


            ## ------------------------------------------------
            ## ≥1 kb FeGenie pool
            ## ------------------------------------------------

            if length >= 1000:

                write_record(
                    out,
                    contig,
                    sequence
                )


                unbinned_1k_n += 1
                unbinned_1k_bp += length


            ## ------------------------------------------------
            ## SemiBin2-eligible but remained unbinned
            ## ------------------------------------------------

            if length >= 2500:

                unbinned_2500_n += 1
                unbinned_2500_bp += length


    ## --------------------------------------------------------
    ## Verify every binned contig belongs to the assembly
    ## --------------------------------------------------------

    missing_from_assembly = (
        set(
            assigned.keys()
        )
        - seen_assembly_contigs
    )


    if missing_from_assembly:

        raise RuntimeError(
            f"{sample}: "
            f"{len(missing_from_assembly)} "
            f"SemiBin2 contigs were not found "
            f"in the Flye assembly."
        )


    if unbinned_1k_n == 0:

        raise RuntimeError(
            f"{sample}: no unbinned ≥1 kb "
            f"contigs were recovered."
        )


    ## --------------------------------------------------------
    ## Summary
    ## --------------------------------------------------------

    summary_rows.append({

        "Sample":
            sample,

        "SemiBin2_MAGs":
            len(
                sample_master
            ),

        "Assembly_contigs":
            total_n,

        "Assembly_bp":
            total_bp,

        "Binned_contigs":
            binned_n,

        "Binned_bp":
            binned_bp,

        "Unbinned_contigs_all":
            unbinned_all_n,

        "Unbinned_bp_all":
            unbinned_all_bp,

        "Unbinned_contigs_ge1000":
            unbinned_1k_n,

        "Unbinned_bp_ge1000":
            unbinned_1k_bp,

        "Unbinned_contigs_ge2500":
            unbinned_2500_n,

        "Unbinned_bp_ge2500":
            unbinned_2500_bp,

        "Pct_assembly_bp_in_MAGs":
            (
                100
                * binned_bp
                / total_bp
            ),

        "FeGenie_input_FASTA":
            str(
                output_fasta
            )
    })


    print(
        f"Assembly:             "
        f"{total_n:,} contigs / "
        f"{total_bp:,} bp"
    )

    print(
        f"In SemiBin2 MAGs:     "
        f"{binned_n:,} contigs / "
        f"{binned_bp:,} bp"
    )

    print(
        f"Unbinned total:       "
        f"{unbinned_all_n:,} contigs / "
        f"{unbinned_all_bp:,} bp"
    )

    print(
        f"Unbinned >=1 kb:      "
        f"{unbinned_1k_n:,} contigs / "
        f"{unbinned_1k_bp:,} bp"
    )

    print(
        f"Unbinned >=2.5 kb:    "
        f"{unbinned_2500_n:,} contigs / "
        f"{unbinned_2500_bp:,} bp"
    )


## ------------------------------------------------------------
## Save tables
## ------------------------------------------------------------

manifest = pd.DataFrame(
    manifest_rows
)


manifest.to_csv(
    MANIFEST_OUT,
    sep="\t",
    index=False
)


summary = pd.DataFrame(
    summary_rows
)


summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)


print()
print(
    "=" * 70
)

print(
    "Unbinned-contig extraction complete"
)

print(
    "=" * 70
)

print(
    summary[
        [
            "Sample",
            "Unbinned_contigs_ge1000",
            "Unbinned_bp_ge1000",
            "Unbinned_contigs_ge2500",
            "Pct_assembly_bp_in_MAGs"
        ]
    ].to_string(
        index=False
    )
)

print()

print(
    f"FeGenie inputs:\n"
    f"{INPUT_ROOT}"
)

print()

print(
    f"Contig manifest:\n"
    f"{MANIFEST_OUT}"
)

print()

print(
    f"Summary:\n"
    f"{SUMMARY_OUT}"
)
