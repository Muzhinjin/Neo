#!/bin/bash

set -euo pipefail

mkdir -p 02_qc/BUSCO

for fasta in 01_genomes/*.fasta; do

    [ -f "$fasta" ] || continue

    name=$(basename "$fasta" .fasta)

    echo "Submitting BUSCO: $name"

    ern jobs submit \
        --threads=16 \
        --memory=128gb \
        --hours=100 \
        --input="$fasta" \
        --module="busco/2.0_c6dec63" \
        --command=busco -- \
        "$fasta" \
        -l sordariomycetes_odb10 \
        -m genome \
        -o "${name}_busco" \
        -c 16

done
