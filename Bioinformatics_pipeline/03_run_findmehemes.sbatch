#!/usr/bin/env bash
#SBATCH --job-name=FindMeHemes
#SBATCH --partition=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/03_findmehemes/logs/findmehemes_%j.out
#SBATCH --error=/home/jeppeos/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow/03_findmehemes/logs/findmehemes_%j.err

set -Eeo pipefail

###############################################################################
## User settings
###############################################################################

PROJECT_DIR="$HOME/methanotrophs/methanotroph_project/jeppe/genome_analysis_workflow"
INPUT_DIR="${PROJECT_DIR}/FeGenie_conda_full_run_20260803/fegenie_results/ORF_calls"
WORK_DIR="${PROJECT_DIR}/03_findmehemes"
FINDMEHEMES_SCRIPT="${WORK_DIR}/FindMeHemes.py"
OUTPUT_ROOT="${WORK_DIR}/results_by_genome"
LOG_DIR="${WORK_DIR}/logs"
TMP_ROOT="${SLURM_TMPDIR:-${WORK_DIR}/tmp/findmehemes_${SLURM_JOB_ID}}"

## MAX_FILES=1 processes one proteome for testing.
## MAX_FILES=0 processes every detected proteome.
MAX_FILES="${MAX_FILES:-0}"

###############################################################################
## Helper functions
###############################################################################

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

die() {
    echo "[$(timestamp)] ERROR: $*" >&2
    exit 1
}

normalise_genome_id() {
    local name="$1"
    name="${name%-proteins.faa}"
    name="${name%.fa}"
    name="${name%.fna}"
    name="${name%.fasta}"
    printf '%s\n' "$name"
}

###############################################################################
## Initial checks
###############################################################################

echo "[$(timestamp)] Starting FindMeHemes workflow."
echo "SLURM job ID:       ${SLURM_JOB_ID:-NA}"
echo "Host:               $(hostname)"
echo "Requested CPUs:     ${SLURM_CPUS_PER_TASK:-NA}"
echo "Maximum files:      ${MAX_FILES} (0 means all)"
echo

[[ -d "$INPUT_DIR" ]] || die "Input directory does not exist: $INPUT_DIR"
[[ -s "$FINDMEHEMES_SCRIPT" ]] || die "FindMeHemes.py is missing or empty: $FINDMEHEMES_SCRIPT"
command -v python3 >/dev/null 2>&1 || die "python3 is not available."
[[ "$MAX_FILES" =~ ^[0-9]+$ ]] || die "MAX_FILES must be a non-negative integer."

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR" "$TMP_ROOT"

echo "Python executable:  $(command -v python3)"
echo "FindMeHemes script: $FINDMEHEMES_SCRIPT"
echo

###############################################################################
## Discover FeGenie proteomes
###############################################################################

mapfile -t ALL_INPUT_FILES < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*-proteins.faa' \
        -size +0c \
        -print |
    LC_ALL=C sort
)

N_DETECTED="${#ALL_INPUT_FILES[@]}"
[[ "$N_DETECTED" -gt 0 ]] || die "No non-empty *-proteins.faa files found in: $INPUT_DIR"

if (( MAX_FILES > 0 )); then
    INPUT_FILES=("${ALL_INPUT_FILES[@]:0:MAX_FILES}")
else
    INPUT_FILES=("${ALL_INPUT_FILES[@]}")
fi

N_SELECTED="${#INPUT_FILES[@]}"

echo "Detected proteomes: $N_DETECTED"
echo "Selected proteomes: $N_SELECTED"
echo

###############################################################################
## Process each proteome
###############################################################################

N_COMPLETED_THIS_RUN=0
N_SKIPPED=0
N_FAILED=0
FAILED_LIST="${WORK_DIR}/failed_genomes_${SLURM_JOB_ID}.txt"
: > "$FAILED_LIST"

for input_file in "${INPUT_FILES[@]}"; do

    input_basename="$(basename "$input_file")"
    genome_id="$(normalise_genome_id "$input_basename")"

    final_dir="${OUTPUT_ROOT}/${genome_id}"
    complete_marker="${final_dir}/FINDMEHEMES_COMPLETE.txt"

    if [[ -s "$complete_marker" ]]; then
        echo "[$(timestamp)] SKIP ${genome_id}: completion marker exists."
        N_SKIPPED=$((N_SKIPPED + 1))
        continue
    fi

    echo "[$(timestamp)] Processing ${genome_id}"

    genome_tmp="${TMP_ROOT}/${genome_id}"
    rm -rf "$genome_tmp"
    mkdir -p "$genome_tmp"

    ## FindMeHemes writes beside the input file. Run on a temporary copy so the
    ## authoritative FeGenie ORF_calls directory remains untouched.
    tmp_input="${genome_tmp}/${genome_id}.faa"

    ## Write each protein sequence on one line. The original FindMeHemes parser
    ## retains FASTA line breaks, which can otherwise hide a motif spanning a
    ## wrapped-line boundary. Headers and protein IDs are preserved unchanged.
    python3 - "$input_file" "$tmp_input" <<'PY_UNWRAP'
import sys
from pathlib import Path

source = Path(sys.argv[1])
target = Path(sys.argv[2])

header = None
sequence_parts = []

with source.open("r", encoding="utf-8") as src, target.open("w", encoding="utf-8") as dst:
    for raw_line in src:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                dst.write(header + "\n")
                dst.write("".join(sequence_parts) + "\n")
            header = line
            sequence_parts = []
        else:
            if header is None:
                raise ValueError(f"Sequence encountered before a FASTA header in {source}")
            sequence_parts.append("".join(line.split()))

    if header is not None:
        dst.write(header + "\n")
        dst.write("".join(sequence_parts) + "\n")

if header is None:
    raise ValueError(f"No FASTA records found in {source}")
PY_UNWRAP

    raw_csv="${tmp_input}.hemeProts.csv"
    raw_fasta="${tmp_input}.hemeProts.fasta"
    tool_log="${genome_tmp}/FindMeHemes.log"

    if ! python3 "$FINDMEHEMES_SCRIPT" \
        -proteins "$tmp_input" \
        > "$tool_log" 2>&1
    then
        echo "[$(timestamp)] FAILED ${genome_id}: FindMeHemes returned a non-zero status." >&2
        echo "$genome_id" >> "$FAILED_LIST"
        N_FAILED=$((N_FAILED + 1))
        continue
    fi

    if [[ ! -s "$raw_csv" ]]; then
        echo "[$(timestamp)] FAILED ${genome_id}: raw CSV was not created." >&2
        echo "$genome_id" >> "$FAILED_LIST"
        N_FAILED=$((N_FAILED + 1))
        continue
    fi

    if [[ ! -f "$raw_fasta" ]]; then
        echo "[$(timestamp)] FAILED ${genome_id}: heme-protein FASTA was not created." >&2
        echo "$genome_id" >> "$FAILED_LIST"
        N_FAILED=$((N_FAILED + 1))
        continue
    fi

    standard_tsv="${genome_tmp}/${genome_id}.findmehemes.tsv"
    summary_tsv="${genome_tmp}/${genome_id}.summary.tsv"

    ## The original tool writes the complete Prodigal FASTA header as ORF.
    ## Create a join-ready protein_id from the first whitespace-delimited token
    ## while retaining the full original header.
    if ! python3 - \
        "$raw_csv" \
        "$tmp_input" \
        "$standard_tsv" \
        "$summary_tsv" \
        "$genome_id" \
        "$input_file" <<'PY'
import csv
import sys
from pathlib import Path

raw_csv = Path(sys.argv[1])
input_fasta = Path(sys.argv[2])
output_tsv = Path(sys.argv[3])
summary_tsv = Path(sys.argv[4])
genome_id = sys.argv[5]
source_proteome = sys.argv[6]

protein_ids = []
with input_fasta.open("r", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith(">"):
            full_header = line[1:].rstrip("\n\r")
            protein_id = full_header.split()[0]
            if not protein_id:
                raise ValueError(f"Empty FASTA identifier in {input_fasta}")
            protein_ids.append(protein_id)

if len(protein_ids) != len(set(protein_ids)):
    raise ValueError(f"Duplicate first-token protein IDs in {input_fasta}")

known_ids = set(protein_ids)
rows = []

with raw_csv.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    expected = {"ORF", "number_of_hemes"}
    if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
        raise ValueError(
            f"Unexpected FindMeHemes CSV header in {raw_csv}: {reader.fieldnames}"
        )

    for row in reader:
        full_header = row["ORF"].strip()
        protein_id = full_header.split()[0]
        heme_count = int(row["number_of_hemes"])

        if protein_id not in known_ids:
            raise ValueError(
                f"FindMeHemes protein ID {protein_id!r} not found in {input_fasta}"
            )
        if heme_count < 1:
            raise ValueError(
                f"Non-positive motif count for {protein_id}: {heme_count}"
            )

        rows.append(
            {
                "genome": genome_id,
                "protein_id": protein_id,
                "full_header": full_header,
                "number_of_hemes": heme_count,
                "source_proteome": source_proteome,
            }
        )

with output_tsv.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = [
        "genome",
        "protein_id",
        "full_header",
        "number_of_hemes",
        "source_proteome",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "genome",
            "input_proteins",
            "heme_positive_proteins",
            "total_predicted_hemes",
            "source_proteome",
        ]
    )
    writer.writerow(
        [
            genome_id,
            len(protein_ids),
            len(rows),
            sum(row["number_of_hemes"] for row in rows),
            source_proteome,
        ]
    )
PY
    then
        echo "[$(timestamp)] FAILED ${genome_id}: output normalization failed." >&2
        echo "$genome_id" >> "$FAILED_LIST"
        N_FAILED=$((N_FAILED + 1))
        continue
    fi

    input_proteins="$(grep -c '^>' "$input_file" || true)"
    heme_positive="$(awk 'NR > 1 { n++ } END { print n + 0 }' "$standard_tsv")"
    total_hemes="$(awk -F '\t' 'NR > 1 { n += $4 } END { print n + 0 }' "$standard_tsv")"

    publish_tmp="${OUTPUT_ROOT}/.${genome_id}.publish_${SLURM_JOB_ID}"
    rm -rf "$publish_tmp"
    mkdir -p "$publish_tmp"

    cp "$raw_csv" "${publish_tmp}/${genome_id}.hemeProts.raw.csv"
    cp "$raw_fasta" "${publish_tmp}/${genome_id}.hemeProts.faa"
    cp "$standard_tsv" "${publish_tmp}/${genome_id}.findmehemes.tsv"
    cp "$summary_tsv" "${publish_tmp}/${genome_id}.summary.tsv"
    cp "$tool_log" "${publish_tmp}/FindMeHemes.log"

    cat > "${publish_tmp}/FINDMEHEMES_COMPLETE.txt" <<MARKER
FindMeHemes completed successfully.
completion_time=$(timestamp)
job_id=${SLURM_JOB_ID}
host=$(hostname)
genome=${genome_id}
input_file=${input_file}
input_proteins=${input_proteins}
heme_positive_proteins=${heme_positive}
total_predicted_hemes=${total_hemes}
MARKER

    rm -rf "$final_dir"
    mv "$publish_tmp" "$final_dir"
    rm -rf "$genome_tmp"

    N_COMPLETED_THIS_RUN=$((N_COMPLETED_THIS_RUN + 1))
    echo "[$(timestamp)] COMPLETE ${genome_id}: ${heme_positive}/${input_proteins} proteins, ${total_hemes} motifs."
done

###############################################################################
## Build combined output tables from all completed genomes
###############################################################################

COMBINED_HITS_TMP="${WORK_DIR}/.findmehemes_all_hits_${SLURM_JOB_ID}.tsv"
COMBINED_SUMMARY_TMP="${WORK_DIR}/.findmehemes_genome_summary_${SLURM_JOB_ID}.tsv"

printf 'genome\tprotein_id\tfull_header\tnumber_of_hemes\tsource_proteome\n' > "$COMBINED_HITS_TMP"

while IFS= read -r result_file; do
    tail -n +2 "$result_file" >> "$COMBINED_HITS_TMP"
done < <(
    find "$OUTPUT_ROOT" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '*.findmehemes.tsv' \
        -print |
    LC_ALL=C sort
)

printf 'genome\tinput_proteins\theme_positive_proteins\ttotal_predicted_hemes\tsource_proteome\n' > "$COMBINED_SUMMARY_TMP"

while IFS= read -r summary_file; do
    tail -n +2 "$summary_file" >> "$COMBINED_SUMMARY_TMP"
done < <(
    find "$OUTPUT_ROOT" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '*.summary.tsv' \
        -print |
    LC_ALL=C sort
)

mv "$COMBINED_HITS_TMP" "${WORK_DIR}/findmehemes_all_hits.tsv"
mv "$COMBINED_SUMMARY_TMP" "${WORK_DIR}/findmehemes_genome_summary.tsv"

N_TOTAL_COMPLETE="$(
    find "$OUTPUT_ROOT" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name 'FINDMEHEMES_COMPLETE.txt' |
    wc -l
)"

N_TOTAL_HITS="$(
    awk 'NR > 1 { n++ } END { print n + 0 }' \
        "${WORK_DIR}/findmehemes_all_hits.tsv"
)"

###############################################################################
## Final status
###############################################################################

echo
echo "===================================================================="
echo "FindMeHemes run summary"
echo "===================================================================="
echo "Detected proteomes:       $N_DETECTED"
echo "Selected for this run:    $N_SELECTED"
echo "Completed this run:       $N_COMPLETED_THIS_RUN"
echo "Skipped as complete:      $N_SKIPPED"
echo "Failed this run:          $N_FAILED"
echo "Total completed genomes:  $N_TOTAL_COMPLETE"
echo "Combined positive hits:   $N_TOTAL_HITS"
echo "Combined hit table:       ${WORK_DIR}/findmehemes_all_hits.tsv"
echo "Genome summary table:     ${WORK_DIR}/findmehemes_genome_summary.tsv"
echo "===================================================================="
echo

if (( N_FAILED > 0 )); then
    echo "[$(timestamp)] ERROR: ${N_FAILED} genome(s) failed. See: $FAILED_LIST" >&2
    exit 1
fi

rm -f "$FAILED_LIST"

if (( N_TOTAL_COMPLETE == N_DETECTED )); then
    cat > "${WORK_DIR}/FINDMEHEMES_ALL_COMPLETE.txt" <<MARKER
All detected FeGenie proteomes were processed successfully.
completion_time=$(timestamp)
job_id=${SLURM_JOB_ID}
detected_proteomes=${N_DETECTED}
completed_proteomes=${N_TOTAL_COMPLETE}
combined_positive_hits=${N_TOTAL_HITS}
MARKER

    echo "[$(timestamp)] All ${N_DETECTED} proteomes are complete."
else
    echo "[$(timestamp)] Partial/test run complete: ${N_TOTAL_COMPLETE}/${N_DETECTED} genomes currently complete."
fi
