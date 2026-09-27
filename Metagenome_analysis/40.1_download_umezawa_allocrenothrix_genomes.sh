#!/bin/bash
set -euo pipefail

## ================================================================== ##
## STAGE 83
##
## DOWNLOAD UMEZAWA Allocrenothrix / Crenothrix GENOMES FROM DDBJ
##
## AF98:
##   AP041007-AP041015  = 9 replicons
##
## MI19235:
##   AP041016-AP041021  = 6 replicons
##
## DDBJ getentry supports accession ranges and:
##
##   format=fasta
##   filetype=text
##
## One FASTA is therefore downloaded per strain.
## ================================================================== ##


PROJECT="$HOME/methanotrophs/methanotroph_project/jeppe/metagenome"

OUTDIR="$PROJECT/comparative_analysis/final_species_phylogeny/00_umezawa_genomes"

mkdir -p "$OUTDIR"


AF98_OUT="$OUTDIR/UME_AF98_Allocrenothrix_methanica.fna"

MI_OUT="$OUTDIR/UME_MI19235_Allocrenothrix_methanica.fna"


download_range () {

    RANGE="$1"
    OUT="$2"

    URL="https://getentry.ddbj.nig.ac.jp/getentry/na/${RANGE}?format=fasta&filetype=text&limit=20"

    echo
    echo "Downloading:"
    echo "  $RANGE"
    echo

    curl \
        -fL \
        --connect-timeout 20 \
        --max-time 300 \
        --retry 4 \
        --retry-delay 3 \
        "$URL" \
        -o "$OUT"


    if [[ ! -s "$OUT" ]]
    then
        echo "ERROR: empty download:"
        echo "$OUT"
        exit 1
    fi


    if ! grep -q '^>' "$OUT"
    then
        echo "ERROR: downloaded file is not FASTA:"
        echo "$OUT"
        echo
        head -n 20 "$OUT"
        exit 1
    fi
}


## ================================================================== ##
## Download AF98
## ================================================================== ##

download_range \
    "AP041007-AP041015" \
    "$AF98_OUT"


## ================================================================== ##
## Download MI19235
## ================================================================== ##

download_range \
    "AP041016-AP041021" \
    "$MI_OUT"


## ================================================================== ##
## Give headers strain-specific prefixes
## ================================================================== ##

sed -i 's/^>/>AF98_/' "$AF98_OUT"

sed -i 's/^>/>MI19235_/' "$MI_OUT"


## ================================================================== ##
## QC
## ================================================================== ##

N_AF98=$(grep -c '^>' "$AF98_OUT")

N_MI=$(grep -c '^>' "$MI_OUT")


echo
echo "============================================================"
echo "STAGE 83 DOWNLOAD QC"
echo "============================================================"
echo

echo "AF98:"
echo "  FASTA records: $N_AF98"
echo "  file: $AF98_OUT"

echo

echo "MI19235:"
echo "  FASTA records: $N_MI"
echo "  file: $MI_OUT"

echo


if [[ "$N_AF98" -ne 9 ]]
then
    echo "ERROR: expected 9 AF98 replicons; found $N_AF98."
    exit 1
fi


if [[ "$N_MI" -ne 6 ]]
then
    echo "ERROR: expected 6 MI19235 replicons; found $N_MI."
    exit 1
fi


## ================================================================== ##
## Total sequence lengths
## ================================================================== ##

python - \
    "$AF98_OUT" \
    "$MI_OUT" \
<<'PY'
from pathlib import Path
import sys


def fasta_stats(path):

    records = []
    name = None
    seq = []


    with open(path) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if name is not None:
                    records.append(
                        (
                            name,
                            len(
                                "".join(seq)
                            )
                        )
                    )

                name = line[1:].split()[0]
                seq = []

            else:
                seq.append(line)


        if name is not None:

            records.append(
                (
                    name,
                    len(
                        "".join(seq)
                    )
                )
            )


    return records


for filename in sys.argv[1:]:

    path = Path(filename)

    records = fasta_stats(path)

    total = sum(
        length
        for _,
        length
        in records
    )


    print(path.name)

    print(
        f"  replicons: {len(records)}"
    )

    print(
        f"  total bp:  {total:,}"
    )


    for name, length in records:

        print(
            f"    {name:<35} {length:>12,}"
        )


    print()
PY


echo "Stage 83 complete."
