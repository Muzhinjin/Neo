
cat > run_iqtree.sh <<'EOF'
#!/bin/bash

set -euo pipefail

echo "IQ-TREE version:"
iqtree --version

echo "Starting Neopestalotiopsis BUSCO phylogenomic analysis"

iqtree \
    -s SUPERMATRIX.fasta \
    -p SUPERMATRIX.partitions_AA.txt \
    -m MFP+MERGE \
    -B 1000 \
    --alrt 1000 \
    --bnni \
    -T 16 \
    --prefix Neopestalotiopsis_BUSCO

echo "IQ-TREE analysis completed"
EOF
