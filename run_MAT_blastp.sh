#!/bin/bash
set -euo pipefail

INPUT_DIR="Protein"
REF_DIR="$HOME/MAT_analysis/references"
OUT_DIR="$HOME/MAT_analysis/BLASTP_results"

DB="$REF_DIR/MAT_reference_db"

mkdir -p "$OUT_DIR"

for file in "$INPUT_DIR"/*_proteins.faa
do
    [ -f "$file" ] || continue

    sample=$(basename "$file" _proteins.faa)

    echo "Screening $sample"

    blastp \
        -query "$file" \
        -db "$DB" \
        -evalue 1e-5 \
        -seg yes \
        -max_target_seqs 20 \
        -outfmt "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore" \
        -num_threads 4 \
        -out "$OUT_DIR/${sample}_MAT_hits.tsv"

done

echo "BLASTP screening completed."
