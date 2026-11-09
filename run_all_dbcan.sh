
#!/bin/bash

set -euo pipefail

module load dbcan/1.0_894582a

# Replace with the actual dbCAN database directory
DB="/PATH/TO/dbCAN_DATABASE"

INPUT="Protein"
OUT="dbCAN_results"

mkdir -p "$OUT"

for file in "$INPUT"/*_proteins.faa
do
    [ -f "$file" ] || continue

    sample=$(basename "$file" _proteins.faa)

    echo "===================================="
    echo "Processing: $sample"
    echo "===================================="

    run_dbcan CAZyme_annotation \
        --input_raw_data "$file" \
        --mode protein \
        --output_dir "$OUT/$sample" \
        --db_dir "$DB" \
        --threads 8

    echo "Completed: $sample"

done

echo "All Neopestalotiopsis proteomes processed."

