#!/bin/bash
set -euo pipefail

echo "=========================================="
echo " AUGUSTUS Gene Prediction"
echo " Species model: fusarium_graminearum"
echo "=========================================="

# Check AUGUSTUS installation
if ! command -v augustus &> /dev/null; then
    echo "ERROR: AUGUSTUS executable not found."
    exit 1
fi

# Species model
SPECIES="fusarium_graminearum"

# Output directory
OUTDIR="AUGUSTUS_results"
mkdir -p "$OUTDIR"

# Check species model availability
echo "Checking AUGUSTUS species model..."

augustus --species="$SPECIES" --help > /dev/null 2>&1 || {
    echo "ERROR: Species model $SPECIES is unavailable."
    exit 1
}

# Process all FASTA files
COUNT=0

for GENOME in *.fasta; do

    [ -f "$GENOME" ] || continue

    BASENAME="${GENOME%.fasta}"

    echo ""
    echo "Processing genome: $GENOME"

    augustus \
        --species="$SPECIES" \
        --gff3=on \
        --protein=on \
        --codingseq=on \
        "$GENOME" \
        > "$OUTDIR/${BASENAME}.gff3"

    echo "Completed: $BASENAME"

    COUNT=$((COUNT + 1))

done

echo ""
echo "=========================================="
echo "AUGUSTUS analysis completed"
echo "Genomes processed: $COUNT"
echo "Results directory: $OUTDIR"
echo "=========================================="

