# Neo
#Genome statistics

assembly-stats *.fasta > genome_summary.txt
ACC=GCA_023628395.1 && \
datasets download genome accession $ACC --include genome,gff3,rna,cds,protein,seq-report --filename ${ACC}.zip && \
unzip -o ${ACC}.zip -d ${ACC}_data && \
cd ${ACC}_data/Neopestalotiopsis37Mncbi_dataset/data/${ACC}/ && \
for f in *; do mv "$f" "${ACC}_$f"; done && \
mv * ../../../../ && \
cd ../../../../ && \
rm -rf ${ACC}_data && \
echo "✅ Files downloaded, extracted, and renamed for $ACC"
# THMMM
biolib run DTU/DeepTMHMM --fasta 104contigs_effectors.fa
