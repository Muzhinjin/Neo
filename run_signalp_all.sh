#!/bin/bash

set -euo pipefail

for fasta in *.aa; do

    [ -f "$fasta" ] || continue

    isolate="${fasta%.aa}"
    prefix="${isolate}_signalp"

    echo "Processing $fasta ..."

    # 1. Run SignalP
    signalp \
        -fasta "$fasta" \
        -format short \
        -org euk \
        -prefix "$prefix"

    # 2. Extract signal peptide IDs
    awk '$2 ~ /^SP/ {print $1}' \
        "${prefix}_summary.signalp5" \
        > "${isolate}_sp_ids.txt"

    # 3. Extract signal peptide protein sequences
    seqkit grep \
        -f "${isolate}_sp_ids.txt" \
        "$fasta" \
        > "${isolate}_signalp_only.fa"

    echo "$isolate: $(grep -c '^>' "${isolate}_signalp_only.fa" || true) signal peptide proteins identified."

done

echo "All SignalP analyses completed!"
