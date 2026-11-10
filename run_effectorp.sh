#!/bin/bash

set -u

mkdir -p EffectorP_results

for f in *_proteins_signalp_only.fa; do

    [ -s "$f" ] || continue

    base="${f%_proteins_signalp_only.fa}"

    echo "Processing: $f"

    EffectorP.py -f \
        -i "$f" \
        -o "EffectorP_results/${base}_predictions.tsv" \
        -E "EffectorP_results/${base}_effectors.fa" \
        -N "EffectorP_results/${base}_non_effectors.fa"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: $f"
    else
        echo "FAILED: $f" >> EffectorP_results/failed_genomes.txt
    fi

done

echo "EffectorP processing completed."
