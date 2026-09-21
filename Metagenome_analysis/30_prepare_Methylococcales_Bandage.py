#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import pandas as pd
import gzip
import re
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

METH_TABLE = (
    PROJECT
    / "semibin2_analysis"
    / "Methylococcales_MAGs_SemiBin2.tsv"
)

OUTDIR = (
    PROJECT
    / "semibin2_analysis"
    / "Bandage_Methylococcales"
)


## ------------------------------------------------------------
## Start clean
## ------------------------------------------------------------

if OUTDIR.exists():

    shutil.rmtree(
        OUTDIR
    )


OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


## ------------------------------------------------------------
## Read Methylococcales table
## ------------------------------------------------------------

mag_table = pd.read_csv(
    METH_TABLE,
    sep="\t"
)


if len(mag_table) != 16:

    raise RuntimeError(
        f"Expected 16 Methylococcales MAGs, "
        f"found {len(mag_table)}."
    )


## ------------------------------------------------------------
## Open plain/gzip FASTA
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

        for line in handle:

            line = line.rstrip()


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
## Write decompressed MAG FASTA
## ------------------------------------------------------------

def write_plain_fasta(
    source,
    destination
):

    with open(
        destination,
        "w"
    ) as out:

        for name, sequence in fasta_records(
            source
        ):

            out.write(
                f">{name}\n"
            )


            for start in range(
                0,
                len(sequence),
                80
            ):

                out.write(
                    sequence[
                        start:start + 80
                    ]
                    + "\n"
                )


## ------------------------------------------------------------
## Read MAG contig IDs
## ------------------------------------------------------------

def get_mag_contigs(path):

    return [
        name
        for name, sequence
        in fasta_records(path)
    ]


## ------------------------------------------------------------
## Read Flye assembly_info.txt
##
## Expected Flye columns include:
##
## #seq_name
## length
## cov.
## circ.
## repeat
## mult.
## alt_group
## graph_path
##
## We parse this manually because the first column begins '#'.
## ------------------------------------------------------------

def read_assembly_info(path):

    path = Path(path)


    if not path.exists():

        raise RuntimeError(
            f"Missing Flye assembly_info.txt:\n"
            f"{path}"
        )


    info = pd.read_csv(
        path,
        sep="\t",
        dtype=str
    )


    ## --------------------------------------------------------
    ## Normalize first column name
    ## --------------------------------------------------------

    first_col = info.columns[0]


    if first_col.startswith("#"):

        info = info.rename(
            columns={
                first_col:
                    first_col.lstrip("#")
            }
        )


    ## --------------------------------------------------------
    ## Find sequence-name column
    ## --------------------------------------------------------

    possible_seq_cols = [
        "seq_name",
        "sequence",
        "contig",
        "name"
    ]


    seq_col = None


    for col in possible_seq_cols:

        if col in info.columns:

            seq_col = col
            break


    if seq_col is None:

        raise RuntimeError(
            f"Could not identify sequence-name column "
            f"in:\n{path}\n"
            f"Columns:\n{list(info.columns)}"
        )


    ## --------------------------------------------------------
    ## Find graph_path column
    ## --------------------------------------------------------

    possible_path_cols = [
        "graph_path",
        "graph-path",
        "path"
    ]


    path_col = None


    for col in possible_path_cols:

        if col in info.columns:

            path_col = col
            break


    if path_col is None:

        raise RuntimeError(
            f"Could not identify graph_path column "
            f"in:\n{path}\n"
            f"Columns:\n{list(info.columns)}"
        )


    mapping = {}


    for _, row in info.iterrows():

        contig = str(
            row[
                seq_col
            ]
        )

        graph_path = str(
            row[
                path_col
            ]
        )


        mapping[
            contig
        ] = graph_path


    return mapping


## ------------------------------------------------------------
## Read Flye GFA
## ------------------------------------------------------------

def read_gfa(path):

    path = Path(path)


    if not path.exists():

        raise RuntimeError(
            f"Missing Flye GFA:\n"
            f"{path}"
        )


    headers = []

    segments = {}

    links = []


    with open(
        path,
        "r"
    ) as handle:

        for line in handle:

            line = line.rstrip(
                "\n"
            )


            if not line:

                continue


            fields = line.split(
                "\t"
            )


            if fields[0] == "H":

                headers.append(
                    line
                )


            elif fields[0] == "S":

                if len(fields) < 2:

                    continue


                segment_id = fields[1]


                segments[
                    segment_id
                ] = line


            elif fields[0] == "L":

                if len(fields) < 5:

                    continue


                from_segment = fields[1]

                to_segment = fields[3]


                links.append(
                    (
                        from_segment,
                        to_segment,
                        line
                    )
                )


    return (
        headers,
        segments,
        links
    )


## ------------------------------------------------------------
## Normalize GFA segment IDs
##
## assembly_info graph paths and GFA segment names may differ
## slightly in formatting, e.g.
##
##   graph_path:  12,-14,15
##   GFA IDs:     12,14,15
##
## or occasionally edge_12 style naming.
## ------------------------------------------------------------

def segment_aliases(segment_id):

    aliases = {
        segment_id
    }


    x = segment_id.lstrip(
        "+-"
    )


    aliases.add(
        x
    )


    if x.startswith(
        "edge_"
    ):

        aliases.add(
            x[len("edge_"):]
        )


    if x.startswith(
        "contig_"
    ):

        aliases.add(
            x[len("contig_"):]
        )


    return aliases


## ------------------------------------------------------------
## Build alias -> actual GFA segment lookup
## ------------------------------------------------------------

def make_segment_lookup(segments):

    lookup = {}


    for segment_id in segments:

        for alias in segment_aliases(
            segment_id
        ):

            if alias in lookup:

                if lookup[alias] != segment_id:

                    ## Ambiguous alias:
                    ## remove it rather than guess.
                    lookup[
                        alias
                    ] = None

            else:

                lookup[
                    alias
                ] = segment_id


    return lookup


## ------------------------------------------------------------
## Parse graph_path from assembly_info.txt
##
## Flye normally stores comma-separated oriented edge IDs:
##
## 1,-2,3
##
## We only need segment membership here, so orientation signs
## are removed.
## ------------------------------------------------------------

## ------------------------------------------------------------
## Parse graph_path from assembly_info.txt
##
## Examples from Flye:
##
##   127
##   369,1495
##   *,105,*
##
## "*" indicates an unresolved/open end and is NOT a GFA
## segment. It is therefore ignored.
## ------------------------------------------------------------

def parse_graph_path(
    graph_path,
    segment_lookup
):

    graph_path = str(
        graph_path
    ).strip()


    if (
        graph_path == ""
        or graph_path.lower() == "nan"
        or graph_path == "*"
    ):

        return []


    tokens = re.split(
        r"[,;\s]+",
        graph_path
    )


    resolved = []


    for token in tokens:

        token = token.strip()


        ## ----------------------------------------------------
        ## Ignore empty tokens and Flye "*" placeholders
        ## ----------------------------------------------------

        if (
            not token
            or token == "*"
        ):

            continue


        ## ----------------------------------------------------
        ## Remove orientation signs
        ## ----------------------------------------------------

        clean = token.lstrip(
            "+-"
        )


        if (
            clean == ""
            or clean == "*"
        ):

            continue


        ## ----------------------------------------------------
        ## Try direct aliases
        ## ----------------------------------------------------

        candidates = [
            token,
            clean,
            f"edge_{clean}"
        ]


        actual = None


        for candidate in candidates:

            if candidate in segment_lookup:

                value = segment_lookup[
                    candidate
                ]


                if value is not None:

                    actual = value
                    break


        ## ----------------------------------------------------
        ## Fallback:
        ## extract numeric graph edge ID
        ## ----------------------------------------------------

        if actual is None:

            match = re.search(
                r"(\d+)",
                clean
            )


            if match:

                numeric = match.group(
                    1
                )


                for candidate in [
                    numeric,
                    f"edge_{numeric}"
                ]:

                    if candidate in segment_lookup:

                        value = segment_lookup[
                            candidate
                        ]


                        if value is not None:

                            actual = value
                            break


        ## ----------------------------------------------------
        ## Fail on genuine unresolved graph nodes
        ## ----------------------------------------------------

        if actual is None:

            raise RuntimeError(
                f"Could not resolve graph-path token "
                f"'{token}' to a GFA segment."
            )


        resolved.append(
            actual
        )


    return resolved

    ## --------------------------------------------------------
    ## Standard Flye format is comma-separated.
    ##
    ## Also tolerate whitespace and semicolons.
    ## --------------------------------------------------------

    tokens = re.split(
        r"[,;\s]+",
        graph_path
    )


    resolved = []


    for token in tokens:

        token = token.strip()


        if not token:

            continue


        ## Remove orientation symbols

        clean = token.lstrip(
            "+-"
        )


        candidates = [
            token,
            clean,
            f"edge_{clean}"
        ]


        actual = None


        for candidate in candidates:

            if candidate in segment_lookup:

                value = segment_lookup[
                    candidate
                ]


                if value is not None:

                    actual = value
                    break


        if actual is None:

            ## ------------------------------------------------
            ## Final fallback:
            ## extract an integer edge ID from strings such as
            ## edge_12+ or -12.
            ## ------------------------------------------------

            match = re.search(
                r"(\d+)",
                clean
            )


            if match:

                numeric = match.group(
                    1
                )


                for candidate in [
                    numeric,
                    f"edge_{numeric}"
                ]:

                    if candidate in segment_lookup:

                        value = segment_lookup[
                            candidate
                        ]


                        if value is not None:

                            actual = value
                            break


        if actual is None:

            raise RuntimeError(
                f"Could not resolve graph-path token "
                f"'{token}' to a GFA segment."
            )


        resolved.append(
            actual
        )


    return resolved


## ------------------------------------------------------------
## Write selected GFA subgraph
## ------------------------------------------------------------

def write_gfa_subset(
    destination,
    headers,
    segments,
    links,
    selected
):

    selected = set(
        selected
    )


    with open(
        destination,
        "w"
    ) as out:

        if headers:

            for line in headers:

                out.write(
                    line + "\n"
                )

        else:

            out.write(
                "H\tVN:Z:1.0\n"
            )


        ## ----------------------------------------------------
        ## Segments
        ## ----------------------------------------------------

        for segment_id in sorted(
            selected
        ):

            if segment_id in segments:

                out.write(
                    segments[
                        segment_id
                    ]
                    + "\n"
                )


        ## ----------------------------------------------------
        ## Links
        ## ----------------------------------------------------

        for (
            node1,
            node2,
            line
        ) in links:

            if (
                node1 in selected
                and
                node2 in selected
            ):

                out.write(
                    line + "\n"
                )


## ------------------------------------------------------------
## Adjacency
## ------------------------------------------------------------

def make_adjacency(
    links
):

    adjacency = defaultdict(
        set
    )


    for (
        node1,
        node2,
        line
    ) in links:

        adjacency[
            node1
        ].add(
            node2
        )

        adjacency[
            node2
        ].add(
            node1
        )


    return adjacency


## ------------------------------------------------------------
## Cache Flye information per barcode
## ------------------------------------------------------------

sample_cache = {}


def load_sample(sample):

    if sample in sample_cache:

        return sample_cache[
            sample
        ]


    flye_dir = (
        PROJECT
        / "flye_results"
        / sample
    )


    info_file = (
        flye_dir
        / "assembly_info.txt"
    )


    gfa_file = (
        flye_dir
        / "assembly_graph.gfa"
    )


    assembly_mapping = read_assembly_info(
        info_file
    )


    (
        headers,
        segments,
        links
    ) = read_gfa(
        gfa_file
    )


    segment_lookup = make_segment_lookup(
        segments
    )


    adjacency = make_adjacency(
        links
    )


    sample_cache[
        sample
    ] = {

        "assembly_mapping":
            assembly_mapping,

        "headers":
            headers,

        "segments":
            segments,

        "links":
            links,

        "segment_lookup":
            segment_lookup,

        "adjacency":
            adjacency,

        "info_file":
            info_file,

        "gfa_file":
            gfa_file
    }


    return sample_cache[
        sample
    ]


## ------------------------------------------------------------
## Main analysis
## ------------------------------------------------------------

summary_rows = []


for _, row in mag_table.iterrows():

    genome_id = str(
        row[
            "Genome_ID"
        ]
    )


    sample = str(
        row[
            "Sample"
        ]
    )


    genus = (
        "Unclassified"
        if pd.isna(
            row.get(
                "Genus"
            )
        )
        else str(
            row.get(
                "Genus"
            )
        )
    )


    species = (
        ""
        if pd.isna(
            row.get(
                "Species"
            )
        )
        else str(
            row.get(
                "Species"
            )
        )
    )


    source_fasta = Path(
        str(
            row[
                "Original_path"
            ]
        )
    )


    if not source_fasta.exists():

        raise RuntimeError(
            f"Missing MAG FASTA:\n"
            f"{source_fasta}"
        )


    ## --------------------------------------------------------
    ## Output directory
    ## --------------------------------------------------------

    mag_dir = (
        OUTDIR
        / sample
        / genome_id
    )


    mag_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    ## --------------------------------------------------------
    ## Plain FASTA copy
    ## --------------------------------------------------------

    output_fasta = (
        mag_dir
        / f"{genome_id}.fa"
    )


    write_plain_fasta(
        source_fasta,
        output_fasta
    )


    ## --------------------------------------------------------
    ## MAG contigs
    ## --------------------------------------------------------

    mag_contigs = get_mag_contigs(
        source_fasta
    )


    ## --------------------------------------------------------
    ## Load Flye graph + mapping
    ## --------------------------------------------------------

    flye = load_sample(
        sample
    )


    assembly_mapping = flye[
        "assembly_mapping"
    ]

    headers = flye[
        "headers"
    ]

    segments = flye[
        "segments"
    ]

    links = flye[
        "links"
    ]

    segment_lookup = flye[
        "segment_lookup"
    ]

    adjacency = flye[
        "adjacency"
    ]


    ## --------------------------------------------------------
    ## Convert MAG contigs to graph segments
    ## --------------------------------------------------------

    mag_segments = set()

    contig_rows = []


    missing_from_info = []


    for contig in mag_contigs:

        if contig not in assembly_mapping:

            missing_from_info.append(
                contig
            )

            continue


        graph_path = assembly_mapping[
            contig
        ]


        path_segments = parse_graph_path(
            graph_path,
            segment_lookup
        )


        mag_segments.update(
            path_segments
        )


        contig_rows.append({

            "Contig":
                contig,

            "Graph_path":
                graph_path,

            "N_graph_segments":
                len(
                    path_segments
                ),

            "Graph_segments":
                ",".join(
                    path_segments
                )
        })


    ## --------------------------------------------------------
    ## Validate mapping
    ## --------------------------------------------------------

    if missing_from_info:

        raise RuntimeError(
            f"{genome_id}: "
            f"{len(missing_from_info)} MAG contigs "
            f"were missing from assembly_info.txt.\n"
            f"Examples:\n"
            + "\n".join(
                missing_from_info[:10]
            )
        )


    if len(mag_segments) == 0:

        raise RuntimeError(
            f"{genome_id}: no graph segments "
            f"could be recovered."
        )


    ## --------------------------------------------------------
    ## One-hop neighbour graph
    ## --------------------------------------------------------

    neighbour_segments = set()


    for segment in mag_segments:

        neighbour_segments.update(
            adjacency.get(
                segment,
                set()
            )
        )


    neighbour_only = (
        neighbour_segments
        - mag_segments
    )


    mag_plus_neighbours = (
        mag_segments
        | neighbour_segments
    )


    ## --------------------------------------------------------
    ## MAG-only GFA
    ## --------------------------------------------------------

    mag_only_gfa = (
        mag_dir
        / f"{genome_id}_MAG_only.gfa"
    )


    write_gfa_subset(
        mag_only_gfa,
        headers,
        segments,
        links,
        mag_segments
    )


    ## --------------------------------------------------------
    ## MAG + one-hop GFA
    ## --------------------------------------------------------

    neighbour_gfa = (
        mag_dir
        / f"{genome_id}_plus_1hop.gfa"
    )


    write_gfa_subset(
        neighbour_gfa,
        headers,
        segments,
        links,
        mag_plus_neighbours
    )


    ## --------------------------------------------------------
    ## Contig -> graph-path table
    ## --------------------------------------------------------

    contig_table = (
        mag_dir
        / "MAG_contig_graph_paths.tsv"
    )


    pd.DataFrame(
        contig_rows
    ).to_csv(
        contig_table,
        sep="\t",
        index=False
    )


    ## --------------------------------------------------------
    ## Bandage node classification
    ## --------------------------------------------------------

    node_rows = []


    for segment in sorted(
        mag_plus_neighbours
    ):

        node_rows.append({

            "Segment":
                segment,

            "Status":
                (
                    "MAG"
                    if segment in mag_segments
                    else "Neighbour"
                )
        })


    node_table = (
        mag_dir
        / "Bandage_nodes.tsv"
    )


    pd.DataFrame(
        node_rows
    ).to_csv(
        node_table,
        sep="\t",
        index=False
    )


    ## --------------------------------------------------------
    ## Metadata
    ## --------------------------------------------------------

    metadata = pd.DataFrame(
        [
            {
                "Genome_ID":
                    genome_id,

                "Sample":
                    sample,

                "Genus":
                    genus,

                "Species":
                    species,

                "Completeness":
                    row.get(
                        "Completeness"
                    ),

                "Contamination":
                    row.get(
                        "Contamination"
                    ),

                "N_MAG_contigs":
                    len(
                        mag_contigs
                    ),

                "N_MAG_graph_segments":
                    len(
                        mag_segments
                    ),

                "N_1hop_neighbours":
                    len(
                        neighbour_only
                    )
            }
        ]
    )


    metadata_file = (
        mag_dir
        / "metadata.tsv"
    )


    metadata.to_csv(
        metadata_file,
        sep="\t",
        index=False
    )


    ## --------------------------------------------------------
    ## Summary
    ## --------------------------------------------------------

    summary_rows.append({

        "Sample":
            sample,

        "Genome_ID":
            genome_id,

        "Genus":
            genus,

        "Species":
            species,

        "Completeness":
            row.get(
                "Completeness"
            ),

        "Contamination":
            row.get(
                "Contamination"
            ),

        "N_MAG_contigs":
            len(
                mag_contigs
            ),

        "N_MAG_graph_segments":
            len(
                mag_segments
            ),

        "N_1hop_neighbours":
            len(
                neighbour_only
            ),

        "FASTA":
            str(
                output_fasta
            ),

        "MAG_only_GFA":
            str(
                mag_only_gfa
            ),

        "Neighbour_GFA":
            str(
                neighbour_gfa
            )
    })


## ------------------------------------------------------------
## Master manifest
## ------------------------------------------------------------

summary = pd.DataFrame(
    summary_rows
)


summary = summary.sort_values(
    [
        "Sample",
        "Genus",
        "Genome_ID"
    ]
)


manifest = (
    OUTDIR
    / "Methylococcales_Bandage_manifest.tsv"
)


summary.to_csv(
    manifest,
    sep="\t",
    index=False
)


## ------------------------------------------------------------
## Print
## ------------------------------------------------------------

print()

print(
    "Bandage-ready Methylococcales dataset"
)

print(
    "=" * 120
)


print(

    summary[
        [
            "Sample",
            "Genome_ID",
            "Genus",
            "N_MAG_contigs",
            "N_MAG_graph_segments",
            "N_1hop_neighbours"
        ]
    ].to_string(
        index=False
    )
)


print()

print(
    f"MAGs prepared: "
    f"{len(summary)}"
)

print()

print(
    f"Output directory:\n"
    f"{OUTDIR}"
)

print()

print(
    f"Manifest:\n"
    f"{manifest}"
)
