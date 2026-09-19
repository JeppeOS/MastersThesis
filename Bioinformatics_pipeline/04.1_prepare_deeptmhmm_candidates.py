#!/usr/bin/env python3
from __future__ import annotations
import csv, re
from collections import defaultdict
from pathlib import Path

PROJECT = Path.home()/"methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
FEGENIE_RESULTS = PROJECT/"FeGenie_conda_full_run_20260803"/"fegenie_results"
ORF_DIR = FEGENIE_RESULTS/"ORF_calls"
FEGENIE_CSV = FEGENIE_RESULTS/"FeGenie-geneSummary.csv"
FINDMEHEMES_TSV = PROJECT/"03_findmehemes"/"findmehemes_all_hits.tsv"
SIGNALP_ROOT = PROJECT/"02_signalp"/"results_by_genome"
OUT_ROOT = PROJECT/"04_deeptmhmm"
CANDIDATE_DIR = OUT_ROOT/"candidates_by_genome"
MANIFEST_TSV = OUT_ROOT/"deeptmhmm_candidate_manifest.tsv"
GENOME_MANIFEST_TSV = OUT_ROOT/"deeptmhmm_genome_manifest.tsv"
QC_TSV = OUT_ROOT/"deeptmhmm_candidate_qc.tsv"

def normalise_genome(value: str) -> str:
    """Return the bare genome ID from FeGenie/FASTA-derived names.

    Examples:
        BCRBG_00329.fa-proteins.faa -> BCRBG_00329
        BCRBG_00329.fa              -> BCRBG_00329
        BCRBG_00329                 -> BCRBG_00329

    FeGenie's protein FASTAs use a compound suffix, so suffix removal must
    continue until no recognized suffix remains.
    """
    name = Path(value.strip()).name
    suffixes = ("-proteins.faa", ".fasta", ".fna", ".faa", ".fa")

    changed = True
    while changed:
        changed = False
        for suffix in suffixes:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                changed = True
                break

    return name

def fasta_records(path: Path):
    header, seq_parts = None, []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header, seq_parts = line[1:], []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)

def parse_signalp(path: Path):
    """Parse SignalP 6 prediction_results.txt from Prodigal-derived FASTAs.

    SignalP preserves the complete Prodigal FASTA header before its prediction
    columns. Because that header contains spaces, Prediction is not simply the
    second whitespace-delimited field.
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

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")

            if not line or line.startswith("#"):
                continue

            protein_id = line.split()[0]
            match = prediction_pattern.search(line)

            if match is None:
                raise ValueError(
                    f"Could not parse SignalP line in {path}: {line[:200]!r}"
                )

            prediction = match.group(1)
            trailing = (match.group(8) or "").strip()

            cs_match = re.search(
                r"CS pos:\s*(.*?)(?=(?:\s+Pr:)|$)",
                trailing,
            )

            rows[protein_id] = {
                "signalp_prediction": prediction,
                "signalp_cs_position": (
                    cs_match.group(1).strip() if cs_match else ""
                ),
            }

    if not rows:
        raise ValueError(f"No SignalP predictions parsed from {path}")

    return rows

for p in (FEGENIE_CSV, FINDMEHEMES_TSV):
    if not p.is_file():
        raise SystemExit(f"ERROR: missing required file: {p}")
if not ORF_DIR.is_dir():
    raise SystemExit(f"ERROR: missing ORF directory: {ORF_DIR}")

OUT_ROOT.mkdir(parents=True, exist_ok=True)
CANDIDATE_DIR.mkdir(parents=True, exist_ok=True)

# Fail immediately if compound FeGenie filenames are not normalized correctly.
_normalization_tests = {
    "BCRBG_00329.fa-proteins.faa": "BCRBG_00329",
    "BCRBG_00329.fa": "BCRBG_00329",
    "BCRBG_00329": "BCRBG_00329",
}
for raw_name, expected_name in _normalization_tests.items():
    observed_name = normalise_genome(raw_name)
    if observed_name != expected_name:
        raise SystemExit(
            "ERROR: genome-name normalization self-test failed: "
            f"{raw_name!r} -> {observed_name!r}, expected {expected_name!r}"
        )

fegenie = defaultdict(lambda: defaultdict(list))
with FEGENIE_CSV.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    req = {"category","genome/assembly","orf","HMM","bitscore","bitscore_cutoff"}
    missing = req.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"ERROR: FeGenie CSV missing columns: {sorted(missing)}")
    for row in reader:
        genome_raw = row["genome/assembly"].strip()
        pid = row["orf"].strip()

        # FeGenie output can contain a repeated CSV header row.
        # Ignore that metadata line rather than treating it as a protein.
        if genome_raw in {"genome/assembly", "assembly"} and pid == "orf":
            continue

        if not genome_raw or not pid:
            continue

        g = normalise_genome(genome_raw)
        fegenie[g][pid].append(row)

hemes = defaultdict(dict)
with FINDMEHEMES_TSV.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    req = {"genome","protein_id","number_of_hemes"}
    missing = req.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"ERROR: FindMeHemes TSV missing columns: {sorted(missing)}")
    for row in reader:
        g = normalise_genome(row["genome"])
        hemes[g][row["protein_id"].strip()] = int(row["number_of_hemes"])

proteomes = sorted(p for p in ORF_DIR.glob("*-proteins.faa") if p.is_file() and p.stat().st_size > 0)
if len(proteomes) != 631:
    print(f"WARNING: expected 631 proteomes, detected {len(proteomes)}")

manifest_rows, genome_rows, qc_rows = [], [], []
seen_fg, seen_heme = set(), set()

for proteome in proteomes:
    genome = normalise_genome(proteome.name)
    sp_file = SIGNALP_ROOT/genome/"prediction_results.txt"
    sp = parse_signalp(sp_file) if sp_file.is_file() else {}
    out_faa = CANDIDATE_DIR/f"{genome}.faa"
    n_input = n_candidates = n_fg = n_heme = n_both = n_sp = 0
    with out_faa.open("w", encoding="utf-8") as out:
        for full_header, seq in fasta_records(proteome):
            n_input += 1
            pid = full_header.split()[0]
            fg_hits = fegenie.get(genome, {}).get(pid, [])
            heme_count = hemes.get(genome, {}).get(pid, 0)
            is_fg, is_heme = bool(fg_hits), heme_count > 0
            if not (is_fg or is_heme):
                continue
            n_candidates += 1; n_fg += int(is_fg); n_heme += int(is_heme); n_both += int(is_fg and is_heme)
            if is_fg: seen_fg.add((genome,pid))
            if is_heme: seen_heme.add((genome,pid))
            sp_row = sp.get(pid, {"signalp_prediction":"", "signalp_cs_position":""})
            n_sp += int(bool(sp_row["signalp_prediction"]))
            hmms = sorted({x["HMM"].strip() for x in fg_hits if x["HMM"].strip()})
            cats = sorted({x["category"].strip() for x in fg_hits if x["category"].strip()})
            bits = []
            cuts = []
            for x in fg_hits:
                try: bits.append(float(x["bitscore"]))
                except Exception: pass
                try: cuts.append(float(x["bitscore_cutoff"]))
                except Exception: pass
            source = "FeGenie+FindMeHemes" if is_fg and is_heme else ("FeGenie" if is_fg else "FindMeHemes")
            out.write(f">{full_header}\n")
            for i in range(0, len(seq), 80): out.write(seq[i:i+80]+"\n")
            manifest_rows.append({
                "genome": genome,
                "protein_id": pid,
                "length": len(seq),
                "candidate_source": source,
                "fegenie_positive": int(is_fg),
                "fegenie_HMMs": ";".join(hmms),
                "fegenie_categories": ";".join(cats),
                "fegenie_max_bitscore": max(bits) if bits else "",
                "fegenie_max_cutoff": max(cuts) if cuts else "",
                "findmehemes_positive": int(is_heme),
                "number_of_hemes": heme_count,
                "signalp_prediction": sp_row["signalp_prediction"],
                "signalp_cs_position": sp_row["signalp_cs_position"],
                "candidate_fasta": str(out_faa),
                "source_proteome": str(proteome),
            })
    genome_rows.append({
        "genome": genome, "candidate_fasta": str(out_faa), "input_proteins": n_input,
        "candidate_proteins": n_candidates, "fegenie_candidates": n_fg,
        "findmehemes_candidates": n_heme, "both": n_both
    })
    qc_rows.append({
        "genome": genome, "signalp_file_found": int(sp_file.is_file()),
        "candidate_proteins": n_candidates, "candidates_with_signalp_result": n_sp
    })

    if not sp_file.is_file():
        raise SystemExit(
            f"ERROR: SignalP result file missing for {genome}: {sp_file}"
        )

    if n_sp != n_candidates:
        raise SystemExit(
            f"ERROR: SignalP join incomplete for {genome}: "
            f"{n_sp}/{n_candidates} candidates matched."
        )

expected_fg = {(g,p) for g,d in fegenie.items() for p in d}
expected_heme = {(g,p) for g,d in hemes.items() for p in d}
missing_fg = sorted(expected_fg-seen_fg)
missing_heme = sorted(expected_heme-seen_heme)
if missing_fg:
    raise SystemExit(f"ERROR: {len(missing_fg)} FeGenie IDs were not found in authoritative FASTAs. First: {missing_fg[:10]}")
if missing_heme:
    raise SystemExit(f"ERROR: {len(missing_heme)} FindMeHemes IDs were not found in authoritative FASTAs. First: {missing_heme[:10]}")

fields = ["genome","protein_id","length","candidate_source","fegenie_positive","fegenie_HMMs","fegenie_categories","fegenie_max_bitscore","fegenie_max_cutoff","findmehemes_positive","number_of_hemes","signalp_prediction","signalp_cs_position","candidate_fasta","source_proteome"]
with MANIFEST_TSV.open("w", encoding="utf-8", newline="") as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(manifest_rows)
fields2=["genome","candidate_fasta","input_proteins","candidate_proteins","fegenie_candidates","findmehemes_candidates","both"]
with GENOME_MANIFEST_TSV.open("w", encoding="utf-8", newline="") as h:
    w=csv.DictWriter(h,fieldnames=fields2,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(genome_rows)
fields3=["genome","signalp_file_found","candidate_proteins","candidates_with_signalp_result"]
with QC_TSV.open("w", encoding="utf-8", newline="") as h:
    w=csv.DictWriter(h,fieldnames=fields3,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(qc_rows)

print("DeepTMHMM candidate preparation complete.")
print(f"Proteomes:                    {len(proteomes)}")
print(f"Total candidate proteins:     {len(manifest_rows)}")
print(f"FeGenie only:                 {sum(r['candidate_source']=='FeGenie' for r in manifest_rows)}")
print(f"FindMeHemes only:             {sum(r['candidate_source']=='FindMeHemes' for r in manifest_rows)}")
print(f"Both:                         {sum(r['candidate_source']=='FeGenie+FindMeHemes' for r in manifest_rows)}")
print(f"Candidates with SignalP data: {sum(bool(r['signalp_prediction']) for r in manifest_rows)}")
print(f"Candidate FASTAs:             {CANDIDATE_DIR}")
print(f"Protein manifest:             {MANIFEST_TSV}")
print(f"Genome manifest:              {GENOME_MANIFEST_TSV}")
print(f"QC table:                     {QC_TSV}")
