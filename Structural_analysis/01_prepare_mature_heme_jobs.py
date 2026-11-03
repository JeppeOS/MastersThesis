#!/usr/bin/env python3

import csv
import json
import random
import re
from pathlib import Path


## ================================================================
## Paths
## ================================================================

HOME = Path.home()

PROJECT = (
    HOME
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)

WF = PROJECT / "genome_analysis_workflow"

S18 = WF / "18_structure_prediction"

OLD_FIGURE = (
    S18
    / "08_apo_structure_figure"
)

OLD_MANIFEST = (
    OLD_FIGURE
    / "apo_pair_manifest.tsv"
)

SELECTED_RESULTS = (
    S18
    / "06_quantitative_comparison"
    / "selected_results"
)

ORF_DIR = (
    WF
    / "FeGenie_conda_full_run_20260803"
    / "fegenie_results"
    / "ORF_calls"
)

SIGNALP_DIR = (
    WF
    / "02_signalp"
    / "results_by_genome"
)

OUTROOT = (
    S18
    / "09_mature_heme_structure_figure"
)

FASTA_OUT = OUTROOT / "01_inputs" / "fasta"
JSON_OUT = OUTROOT / "01_inputs" / "json"

MANIFEST_OUT = (
    OUTROOT
    / "01_inputs"
    / "mature_heme_manifest.tsv"
)

ALL_FASTA_OUT = (
    OUTROOT
    / "01_inputs"
    / "all_mature_sequences.faa"
)

BATCH_JSON_OUT = (
    OUTROOT
    / "01_inputs"
    / "alphafold_server_mature_heme_BATCH.json"
)


FASTA_OUT.mkdir(parents=True, exist_ok=True)
JSON_OUT.mkdir(parents=True, exist_ok=True)


## ================================================================
## FASTA parser
## ================================================================

def read_fasta(path):

    header = None
    seq = []

    with open(path) as fh:

        for line in fh:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    yield header, "".join(seq).rstrip("*")

                header = line[1:]
                seq = []

            else:
                seq.append(line.strip())

        if header is not None:
            yield header, "".join(seq).rstrip("*")


## ================================================================
## Read an AlphaFold Server job
## ================================================================

def read_old_job(path):

    raw = json.loads(path.read_text())

    if isinstance(raw, list):

        if len(raw) != 1:
            raise RuntimeError(
                f"Expected exactly 1 job in {path}, found {len(raw)}"
            )

        job = raw[0]

    else:
        job = raw

    proteins = [
        x["proteinChain"]
        for x in job.get("sequences", [])
        if "proteinChain" in x
    ]

    if len(proteins) != 1:
        raise RuntimeError(
            f"Expected one protein chain in {path}; "
            f"found {len(proteins)}"
        )

    seq = proteins[0]["sequence"].rstrip("*")

    return job, seq


## ================================================================
## Read exact 10 old models from figure manifest
## ================================================================

targets = []

with open(OLD_MANIFEST) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        family = row["family"]

        for pair_no in (1, 2):

            model_path = Path(
                row[f"model_{pair_no}"]
            )

            model_stem = model_path.stem

            old_job_path = (
                SELECTED_RESULTS
                / model_stem
                / "job_request.json"
            )

            if not old_job_path.exists():

                raise FileNotFoundError(
                    f"Missing old job request:\n"
                    f"{old_job_path}"
                )

            old_job, sequence = read_old_job(
                old_job_path
            )

            targets.append(
                {
                    "family": family,
                    "pair_no": pair_no,
                    "old_model": model_stem,
                    "old_label": row[f"label_{pair_no}"],
                    "old_job_path": str(old_job_path),
                    "full_sequence": sequence,
                }
            )


if len(targets) != 10:

    raise RuntimeError(
        f"Expected exactly 10 figure proteins; "
        f"found {len(targets)}"
    )


print(
    f"Recovered {len(targets)} exact proteins "
    f"from old figure manifest."
)


## ================================================================
## Map exact old AlphaFold sequences back to ORF protein IDs
## ================================================================

wanted_sequences = {
    x["full_sequence"]
    for x in targets
}

sequence_hits = {
    seq: []
    for seq in wanted_sequences
}


proteome_files = sorted(
    ORF_DIR.glob("*-proteins.faa")
)

if not proteome_files:

    raise FileNotFoundError(
        f"No protein FASTAs found in:\n{ORF_DIR}"
    )


print(
    f"Searching {len(proteome_files)} FeGenie proteomes "
    f"for exact sequence matches..."
)


for fasta in proteome_files:

    for header, seq in read_fasta(fasta):

        if seq in sequence_hits:

            protein_id = header.split()[0]

            sequence_hits[seq].append(
                {
                    "protein_id": protein_id,
                    "source_fasta": str(fasta),
                }
            )


## ================================================================
## Resolve exact protein IDs
## ================================================================

for t in targets:

    hits = sequence_hits[
        t["full_sequence"]
    ]

    if len(hits) == 0:

        raise RuntimeError(
            "\nNo exact ORF match for:\n"
            f"  {t['old_model']}\n"
            f"  length={len(t['full_sequence'])}"
        )

    if len(hits) > 1:

        print()
        print(
            f"WARNING: {t['old_model']} has "
            f"{len(hits)} exact sequence matches:"
        )

        for h in hits:
            print(
                " ",
                h["protein_id"]
            )

        ## Try to use genome token encoded in model name.
        normalized_model = (
            t["old_model"]
            .replace("_", "")
            .lower()
        )

        filtered = []

        for h in hits:

            pid_norm = (
                h["protein_id"]
                .replace("_", "")
                .lower()
            )

            genome_part = (
                h["protein_id"]
                .split("_000000", 1)[0]
            )

            genome_norm = (
                genome_part
                .replace("_", "")
                .lower()
            )

            if genome_norm in normalized_model:
                filtered.append(h)

        if len(filtered) == 1:
            hits = filtered

        else:
            raise RuntimeError(
                f"Could not uniquely resolve sequence for "
                f"{t['old_model']}."
            )

    t["protein_id"] = hits[0]["protein_id"]
    t["source_fasta"] = hits[0]["source_fasta"]


## ================================================================
## Genome ID helper
## ================================================================

def genome_from_protein_id(pid):

    if "_000000" in pid:
        return pid.split("_000000", 1)[0]

    raise RuntimeError(
        f"Cannot determine genome from protein ID: {pid}"
    )


## ================================================================
## SignalP lookup
## ================================================================

def get_signalp_record(genome, protein_id):

    pred_file = (
        SIGNALP_DIR
        / genome
        / "prediction_results.txt"
    )

    if not pred_file.exists():

        raise FileNotFoundError(
            f"SignalP file missing:\n{pred_file}"
        )

    matches = []

    with open(pred_file) as fh:

        for line in fh:

            ## Avoid matching IDs such as _4 against _40.
            if (
                line.startswith(protein_id + " ")
                or
                line.startswith(protein_id + "\t")
            ):
                matches.append(
                    line.rstrip()
                )

    if len(matches) != 1:

        raise RuntimeError(
            f"Expected exactly one SignalP record for "
            f"{protein_id}, found {len(matches)}"
        )

    line = matches[0]

    if re.search(
        r"\sSP\s",
        line
    ):
        prediction = "SP"

    elif re.search(
        r"\sTAT\s",
        line
    ):
        prediction = "TAT"

    elif re.search(
        r"\sLIPO\s",
        line
    ):
        prediction = "LIPO"

    elif re.search(
        r"\sOTHER\s",
        line
    ):
        prediction = "OTHER"

    else:
        prediction = "UNKNOWN"

    cs = re.search(
        r"CS pos:\s*(\d+)-(\d+)",
        line
    )

    if cs:

        cleavage_after = int(
            cs.group(1)
        )

    else:

        cleavage_after = 0

    return (
        prediction,
        cleavage_after,
        line
    )


## ================================================================
## CXXCH motif finder
## ================================================================

def find_cxxch(seq):

    motifs = []

    for m in re.finditer(
        r"C..CH",
        seq
    ):

        ## 1-based mature-sequence coordinates.
        start = m.start() + 1
        end = m.end()

        motifs.append(
            {
                "start": start,
                "end": end,
                "sequence": m.group(0),
            }
        )

    return motifs


## ================================================================
## Process all ten targets
## ================================================================

all_jobs = []

random.seed(20260902)


for t in targets:

    pid = t["protein_id"]

    genome = genome_from_protein_id(
        pid
    )

    full_seq = t["full_sequence"]

    (
        signalp_prediction,
        cleavage_after,
        signalp_line
    ) = get_signalp_record(
        genome,
        pid
    )


    ## ------------------------------------------------------------
    ## Trim ONLY if SignalP gives an explicit cleavage site.
    ## ------------------------------------------------------------

    if cleavage_after > 0:

        mature_seq = full_seq[
            cleavage_after:
        ]

        trimming = (
            f"remove residues 1-{cleavage_after}"
        )

    else:

        mature_seq = full_seq

        trimming = "none"


    motifs = find_cxxch(
        mature_seq
    )

    n_hec = len(
        motifs
    )


    motif_positions = ",".join(
        f"{x['start']}-{x['end']}:{x['sequence']}"
        for x in motifs
    )


    ## ------------------------------------------------------------
    ## New job name
    ## ------------------------------------------------------------

    base = t["old_model"]

    if base.endswith("_apo"):
        base = base[:-4]

    if n_hec > 0:

        job_name = (
            f"{base}_mature_hemeC{n_hec}"
        )

    else:

        job_name = (
            f"{base}_mature_noHEC"
        )


    ## ------------------------------------------------------------
    ## AlphaFold Server input
    ##
    ## Reproduces the SAME HEC setup used previously:
    ##
    ##   CCD_HEC
    ##
    ## one copy per CXXCH motif.
    ## ------------------------------------------------------------

    sequences = [
        {
            "proteinChain": {
                "sequence": mature_seq,
                "count": 1,
                "useStructureTemplate": False,
            }
        }
    ]


    if n_hec > 0:

        sequences.append(
            {
                "ligand": {
                    "ligand": "CCD_HEC",
                    "count": n_hec,
                }
            }
        )


    job = {
        "name": job_name,
        "modelSeeds": [
            str(
                random.randint(
                    1,
                    2_000_000_000
                )
            )
        ],
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 3,
    }


    all_jobs.append(
        job
    )


    ## ------------------------------------------------------------
    ## Individual FASTA
    ## ------------------------------------------------------------

    fasta_path = (
        FASTA_OUT
        / f"{job_name}.faa"
    )

    with open(
        fasta_path,
        "w"
    ) as out:

        out.write(
            f">{job_name} "
            f"protein_id={pid} "
            f"SignalP={signalp_prediction} "
            f"cleavage_after={cleavage_after} "
            f"HEC={n_hec}\n"
        )

        for i in range(
            0,
            len(mature_seq),
            80
        ):
            out.write(
                mature_seq[i:i+80]
                + "\n"
            )


    ## ------------------------------------------------------------
    ## Individual AlphaFold Server JSON
    ## ------------------------------------------------------------

    json_path = (
        JSON_OUT
        / f"{job_name}.json"
    )

    json_path.write_text(
        json.dumps(
            [job],
            indent=2
        )
        + "\n"
    )


    ## ------------------------------------------------------------
    ## Save processed metadata
    ## ------------------------------------------------------------

    t["genome"] = genome
    t["signalp_prediction"] = signalp_prediction
    t["cleavage_after"] = cleavage_after
    t["trimming"] = trimming
    t["mature_sequence"] = mature_seq
    t["mature_length"] = len(mature_seq)
    t["cxxch_count"] = n_hec
    t["cxxch_positions_mature"] = motif_positions
    t["job_name"] = job_name
    t["json_file"] = str(json_path)
    t["fasta_file"] = str(fasta_path)
    t["signalp_line"] = signalp_line


## ================================================================
## Combined FASTA
## ================================================================

with open(
    ALL_FASTA_OUT,
    "w"
) as out:

    for t in targets:

        out.write(
            f">{t['job_name']} "
            f"protein_id={t['protein_id']} "
            f"SignalP={t['signalp_prediction']} "
            f"cleavage_after={t['cleavage_after']} "
            f"HEC={t['cxxch_count']}\n"
        )

        seq = t[
            "mature_sequence"
        ]

        for i in range(
            0,
            len(seq),
            80
        ):
            out.write(
                seq[i:i+80]
                + "\n"
            )


## ================================================================
## Batch AlphaFold Server JSON
## ================================================================

BATCH_JSON_OUT.write_text(
    json.dumps(
        all_jobs,
        indent=2
    )
    + "\n"
)


## ================================================================
## Manifest
## ================================================================

fields = [
    "family",
    "pair_no",
    "old_label",
    "old_model",
    "protein_id",
    "genome",
    "signalp_prediction",
    "cleavage_after",
    "trimming",
    "full_length",
    "mature_length",
    "cxxch_count",
    "cxxch_positions_mature",
    "job_name",
    "fasta_file",
    "json_file",
]


with open(
    MANIFEST_OUT,
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore"
    )

    writer.writeheader()

    for t in targets:

        row = dict(t)

        row["full_length"] = len(
            t["full_sequence"]
        )

        writer.writerow(
            row
        )


## ================================================================
## Console QC
## ================================================================

print()
print("=" * 110)
print("MATURE-SEQUENCE / HEME-C PREPARATION")
print("=" * 110)

header = (
    f"{'Family':10s} "
    f"{'Label':18s} "
    f"{'Protein':39s} "
    f"{'SigP':7s} "
    f"{'CS':>4s} "
    f"{'Full':>5s} "
    f"{'Mature':>6s} "
    f"{'CXXCH':>5s}"
)

print(header)
print("-" * len(header))


for t in targets:

    print(
        f"{t['family']:10s} "
        f"{t['old_label']:18s} "
        f"{t['protein_id']:39s} "
        f"{t['signalp_prediction']:7s} "
        f"{t['cleavage_after']:4d} "
        f"{len(t['full_sequence']):5d} "
        f"{t['mature_length']:6d} "
        f"{t['cxxch_count']:5d}"
    )


print()
print("Outputs:")
print(" ", MANIFEST_OUT)
print(" ", ALL_FASTA_OUT)
print(" ", BATCH_JSON_OUT)
print()
print("Individual JSON jobs:")
for t in targets:
    print(" ", t["json_file"])

print()
print("Done.")
