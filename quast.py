#!/bin/bash

for fasta in *.fasta
do
    name=$(basename "$fasta" .fasta)

    ern jobs submit \
        --name="QUAST_${name}" \
        --threads=32 \
        --memory=128gb \
        --hours=24 \
        --input="$fasta" \
        --module="quast/5.2.0" \
        --command=quast.py -- \
        "$fasta" \
        -o "${name}_quast" \
        -t 32
done
