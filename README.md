module load cluster/hpc
ern cluster create --name=hpc --head-node=login.hpc.ufs.ac.za --node-select=ern --email-address=MuzhinjiN@ufs.ac.za --shared-fs
# Quality


#Loop so many samples
ern jobs submit   --name=Neosamtoolsview   --threads=32   --memory=128gb   --hours=48   --input="*.sam"   --module="samtools/2.1_4854cc4"   --command='
    for f in *.sam; do
        base=${f%.sam}
        samtools view -@ 8 -bS "$f" > "${base}_unsorted.bam"
    done

# change fq.gz t fastq
for file in *_good_2.fq.gz; do
    gzip -dc "$file" > "${file%.fq.gz}.fastq"
done

for file in *_reads.fastq.gz; do
    gzip -dc "$file" > "${file%.fastq.gz}.fastq"
done

ern jobs submit \
  --name=NC100_fastqc \
  --threads=4 \
  --memory=16gb \
  --hours=4 \
  --input="NC100_S423_R1_001.fastq"  --input="NC100_S423_R2_001.fastq" \
  --module='fastqc/1.0_4fa3c56' \
  --command=fastqc \
  -- \
  --threads 4 \
  NC100_S423_R1_001.fastq \
  NC100_S423_R2_001.fastq
    ern jobs submit \
  --name=pilonrun \
  --threads=32 \
  --memory=128gb \
  --hours=48 \
  --input="*_contigs.fasta,*_sorted.bam" \
  --module="java,pilon" \
  --command='
    for ref in *_contigs.fasta; do
        base=${ref%_contigs.fasta}
        bam="${base}_sorted.bam"
        if [[ -f "$bam" ]]; then
            java -Xmx120g -jar pilon.jar \
              --genome "$ref" \
              --bam "$bam" \
              --output "${base}_pilon" \
              --threads 32
        fi
    done
  '
# Spdes
ls *_R1_filtered.fastq.gz | sed 's/_R1_filtered.fastq.gz//' > samples_all.txt
wc -l samples_all.txt
bash submit_spades.sh samples_all.txt
# Remove those wth less than 500bp
for file in *_contigs.fasta
do
    output="${file%_contigs.fasta}.fa"
    seqkit seq -m 500 "$file" > "$output"
    echo "Processed $file -> $output"
done
# Busco
for fasta in *.fasta
do
    name=$(basename "$fasta" .fasta)

    ern jobs submit \
        --name="BUSCO_${name}" \
        --threads=16 \
        --memory=128gb \
        --hours=100 \
        --input="$fasta" \
        --module="busco/1.0_88de6b8" \
        --command=busco -- \
        "$fasta" \
        -l sordariomycetes_odb10 \
        -m genome \
        -o "${name}_busco" \
        -c 32
done

# Augustus
ern jobs submit   --name=AugustusNeo   --threads=16  --memory=128gb   --hours=100   --input="*.fasta"   --input="run_augustus.sh"   --module="augustus/1.0_e86bc41"   --command="bash run_augustus.sh"
# run_augustus.sh


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

# Phylogenomics

ern jobs submit \
    --name=Neopestalotiopsis_IQTree \
    --threads=16 \
    --memory=64gb \
    --hours=48 \
    --input="SUPERMATRIX.fasta" \
    --input="SUPERMATRIX.partitions_AA.txt" \
    --input="run_iqtree.sh" \
    --command="bash run_iqtree.sh"

    # Orthofinder

    orthofinder \
    -f 03_annotation/proteins \
    -S diamond \
    -M msa \
    -A mafft \
    -T iqtree \
    -t 16 \
    -a 8

ern jobs submit
--name=NeoOrthoFinder
--threads=32
--memory=128gb
--hours=100
--input="OrthoFinder_input/*.faa"
--module="orthofinder/1.1_76a008b"
--command=orthofinder --
-f OrthoFinder_input
-t 32
-a 32


 ern jobs submit   --name=NEoOrthoFinder   --threads=32   --memory=128gb   --hours=72   --input="Orthofinder_clean/*_proteins.faa"   --module="orthofinder/1.1_76a008b"   --command=orthofinder -- -f Orthofinder_clean -t 32 -a 32
# Extract sequences

ern jobs submit --name="auspicious-ocelot" --threads=8 --memory=16gb --hours=24 --module='blast/3.0_6bf5a81 python=3.12' -- makeblastdb -in PX112669.fasta -dbtype nucl -out PX112669_db

# MAT gene


ern jobs submit \
  --name="Neo_MAT_BLASTP" \
  --threads=16 \
  --memory=64gb \
  --hours=100 \
  --input="config.sh" \
  --input="01_blastp.sh" \
  --input="MAT_analysis/genomes/" \
  --input="MAT_analysis/references/MAT_refs.faa" \
  --module="blast/3.0_6bf5a81" \
  -- bash 01_blastp.sh


ern jobs submit \
  --name=Neo_MAT_TBLASTN \
  --threads=16 \
  --memory=64gb \
  --hours=100 \
  --input="config.sh" \
  --input="03_tblastn.sh" \
  --input="MAT_analysis/genomes/" \
  --input="MAT_analysis/references/MAT_refs.faa" \
  --module="blast/3.0_6bf5a81" \
  -- bash 03_tblastn.sh

  #FASTANI

  /home/muzhinjin/.local/bin/fastANI
$HOME/.local/bin/fastANI \
--ql genome_list.txt \
--rl genome_list.txt \
--threads 16 \
-o FastANI_all_vs_all.tsv


#fungcde restnce 

bash cyp51_scan.sh Cyp51Areference.fasta *_proteins.faa

# Blast
 ern jobs submit --name="diligent-tiger" --threads=8 --memory=16gb --hours=24 --module='blast/3.0_6bf5a81 python=3.12' -- makeblastdb -in MAT_Neurospora.faa -dbtype prot -parse_seqids -out MAT_reference_db

gatk --java-options "-Xmx4g" HaplotypeCaller -R Neopestalotiopsis_rosae_1902.fasta -I 68_S416_dedup.bam -O 68_S416_g.vcf.gz -ERC BP_RESOLUTION --sample-ploidy 1
gatk CombineGVCFs -R Neopestalotiopsis_rosae_1902.fasta $(for f in *_g.vcf.gz; do echo -V $f; done) -O combined.g.vcf.gz
gatk --java-options "-Xmx4g" HaplotypeCaller -R Neopestalotiopsis_rosae_1902.fasta -I 68_S416_dedup.bam -O 68_S416_g.vcf.gz -ERC BP_RESOLUTION --sample-ploidy 1
gatk CombineGVCFs -R Neopestalotiopsis_rosae_1902.fasta $(for f in *_g.vcf.gz; do echo -V $f; done) -O combined.g.vcf.gz
gatk GenotypeGVCFs  -R Neopestalotiopsis_rosae_1902.fasta -V combined.g.vcf.gz  -O raw_joint.vcf.gz
gatk VariantFiltration -R Neopestalotiopsis_rosae_1902.fasta -V raw_joint.vcf.gz --filter-name "QD_lt_10" --filter-expression "QD < 10.0  --filter-name "FS_gt_55" --filter-expression "FS > 55.0" --filter-name "MQ_lt_45" --filter-expression "MQ < 45.0" --filter-name "ReadPosRankSum_outside" --filter-expression "ReadPosRankSum < -4.0 || ReadPosRankSum > 4.0"

gatk VariantFiltration   -R Neopestalotiopsis_rosae_1902.fasta   -V raw_joint.vcf.gz   --filter-name "QD_lt_10" --filter-expression "QD < 10.0"   --filter-name "FS_gt_55" --filter-expression "FS > 55.0"   --filter-name "MQ_lt_45" --filter-expression "MQ < 45.0"   --filter-name "ReadPosRankSum_outside" --filter-expression "ReadPosRankSum < -4.0 || ReadPosRankSum > 4.0"   -O raw_joint.filtered.vcf.gz

bcftools view -m2 -M2 -v snps filtered_PASS.vcf.gz -O z -o filtered_PASS_biallelic_snps.vcf.gz


java -Xmx64G -jar pilon-1.24.jar --genome flye_assembly/assembly.fasta --frags illumina_mapped.sorted.bam --output polished
java -Xmx64G -jar pilon-1.24.jar --genome /home/muzhinjin/Tika_neosporadis/neosporadis_project/Neopcontigs/genomes/NC65contigs.fasta  --frags NC65_S32_sorted.bam --output pilonresults/NC65_assembly_pilon  

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
/home/muzhinjin/Tika_neosporadis/neosporadis_project/Neopcontigs/Proteins/signalp_only/effector_only/effector_only/./effectornotm.py



#gatk
for bam in *_sorted.bam; do   base=$(basename "$bam" _sorted.bam);   ern jobs submit --name="GATKrD_${base}" --threads=8 --memory=32gb --hours=72 --input "$bam" --module='gatk/1.0_5ea5a38' --command= gatk -- AddOrReplaceReadGroups -I $bam -O ${base}_sorted_RG.bam -RGID $base -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM $base; done

# Fungicides resistance 

samtools view -@ 8 -Sb sample.sam | samtools sort -@ 8 -o sample.sorted.bam
samtools index sample.sorted.bam
bcftools mpileup --threads 8 -f reference_genome.fna *.sorted.bam I bcftools call --threads 8 -mv -Ob -o all_samples_variants.bcf
bcftools filter -e 'QUAL<30 || DP<10' all_samplevariants.bcf -Ov -o all_samples_variants.filtered.vcf

snpEff build -c snpEff.config -gff3 -v Your_Neo_Genome
snpEff Your_Neo_Genome all_samples_variants.filtered.vcf > all_samples_variants.annotated.vcf

# Create a pattern file with your resistance gene names
echo -e "CYP51\nCYP51A\nCYP51B\ntubulin\nTUB2\nSdhB\nSdhC\nSdhD\nCYTB" > resistance_genes.txt

# Filter VCF for these genes
grep -w -f resistance_genes.txt all_samples_variants.annotated.vcf > resistance_candidates.vcf

# Or if you want case-insensitive search
grep -i -w -f resistance_genes.txt all_samples_variants.annotated.vcf > resistance_candidates.vcf


# Filter for Specific Amino Acid Changes
# Search for specific amino acid changes
grep -E "p\.(Tyr137Phe|Y137F|Glu198Ala|E198A|Gly143Ala|G143A)" all_samples_variants.annotated.vcf > exact_mutations.vcf

# Case insensitive and more flexible
grep -i -E "p\.[yY]137[fF]|p\.[eE]198[aA]|p\.[gG]143[aA]" all_samples_variants.annotated.vcf

 # Extract All Missense Mutations in Resistance Genes

# First get variants in resistance genes, then filter for missense
grep -w -f resistance_genes.txt all_samples_variants.annotated.vcf | \
grep "missense_variant" > resistance_missense.vcf

# filter_resistance.awk
#!/bin/awk -f
# Filter for resistance gene mutations with quality control

BEGIN {
    # Define your resistance genes (case insensitive)
    resistance_genes["cyp51"] = 1;
    resistance_genes["tubulin"] = 1;
    resistance_genes["tub2"] = 1;
    resistance_genes["sdhb"] = 1;
    resistance_genes["sdhc"] = 1;
    resistance_genes["sdhd"] = 1;
    resistance_genes["cytb"] = 1;
    
    # Define known resistance mutations
    known_mutations["Y137F"] = 1;
    known_mutations["E198A"] = 1;
    known_mutations["G143A"] = 1;
    # Add more from your Gold List
}

# Skip header lines but print them to output
/^#/ { print; next }

# Process variant lines
{
    # Extract INFO field
    split($8, info, ";");


    # Neopestalotiopsis

# Create a working directory and move all FASTQ files there.
mkdir fungal_project
cd fungal_project
mkdir raw_data trimmed_data qc_reports
mv *_R1_001.fastq.gz *_R2_001.fastq.gz raw_data/

/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query sdhB_Zymoseptoria_tritici_resistance_DNA.fasta  -db all_contigs_db -out cytb_hits.txt -evalue 1e-10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

awk '{
  contig=$2; start=$9; end=$10;
  if (start < end) { strand="plus"; s=start; e=end }
  else { strand="minus"; s=end; e=start }
  print contig, s, e, strand
}' cytb_hits.txt > sdhBb_coords.txt


  while read contig start end strand; do
  /home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastdbcmd \
    -db all_contigs_db \
    -entry "$contig" \
    -range "$start-$end" \
    -strand "$strand" \
    -out ${contig}_${start}-${end}.fasta
done < sdhBb_coords.txt


cat sdhB_Zymoseptoria_tritici_resistance_DNA.fast sdhBb_from_contigs.fasta > sdhBb_all.fasta
sdhB_Zymoseptoria_tritici_resistance_DNA.fasta

awk '
/^>/{name=$0; getline seq; seq=toupper(seq);
    # extract codons
    P225F=substr(seq,673,3);
    N225H=substr(seq,673,3);
    H272R=substr(seq,814,3);
    H272Y=substr(seq,814,3);
    I279V=substr(seq,835,3);
    
    # check resistance
    P225F_status=(P225F ~ /TTT|TTC/) ? "resistant":"sensitive";
    N225H_status=(N225H ~ /CAC/) ? "resistant":"sensitive";
    H272R_status=(H272R ~ /CGT|CGC|CGA|CGG/) ? "resistant":"sensitive";
    H272Y_status=(H272Y ~ /TAT|TAC/) ? "resistant":"sensitive";
    I279V_status=(I279V ~ /GTA|GTC|GTG|GTT/) ? "resistant":"sensitive";
    
    # print report
    printf "%s\nP225F: %s -> %s\nN225H: %s -> %s\nH272R: %s -> %s\nH272Y: %s -> %s\nI279V: %s -> %s\n\n", \
    name, P225F, P225F_status, N225H, N225H_status, H272R, H272R_status, H272Y, H272Y_status, I279V, I279V_status
}'  SdhBb_aligned.fasta > sdhb_resistance_report.txt


# Quality Control and Trimming (fastp)
conda install -c bioconda fastp
cd raw_data
for sample in *_R1_001.fastq.gz
do
  base=$(basename ${sample} _R1_001.fastq.gz)
  fastp -i ${base}_R1_001.fastq.gz -I ${base}_R2_001.fastq.gz \
        -o ../trimmed_data/${base}_R1_trimmed.fastq.gz \
        -O ../trimmed_data/${base}_R2_trimmed.fastq.gz \
        -h ../qc_reports/${base}_fastp.html -j ../qc_reports/${base}_fastp.json
done
fastp -i *_R1.fastq.gz -I *_R2.fastq.gz -o clean_R1.fastq.gz -O clean_R2.fastq.gz -h report.html -j report.json

# Mapping to Reference Genome (BWA + Samtools)
# Index Reference
bwa index ref_genome.fasta

# Align Reads

cd ../trimmed_data
mkdir ../aligned_bam
for sample in 73_S419 80_S420 86_S421 93_S422
do
  echo "Aligning $sample"

  bwa mem -t 4 GCNrosiarefgenome.fna \
    ${sample}_clean_R1_001.fastq.gz \
    ${sample}_clean_R2_001.fastq.gz | \
  samtools view -bS - | \
  samtools sort -o ${sample}_aligned_sorted.bam

  samtools index ${sample}_aligned_sorted.bam
done


bwa mem -t 4 ../ref_genome.fasta \
    ${base}_R1_trimmed.fastq.gz ${base}_R2_trimmed.fastq.gz | \
    samtools view -Sb - > ../aligned_bam/${base}.bam

# For singles 
bwa mem GCNrosiarefgenome.fna 73_S419_clean_R1_001.fastq.gz 73_S419_clean_R2_001.fastq.gz > 73_S419_alignment.sam
# Convert SAM to BAM
samtools view -bS 73_S419_alignment.sam > 73_S419_alignment.bam

# Sort the BAM file
samtools sort 73_S419_alignment.bam -o 73_S419_alignment_sorted.bam

# Index bam files
samtools index 73_S419_alignment_sorted.bam


    
# Sort and Index

cd ../aligned_bam
samtools sort -o sorted_${bam} $bam
samtools index sorted_${bam}

# Variant Calling (bcftools)
conda install -c bioconda bcftools

# Create combined pileup
samtools mpileup -uf reference/ref_genome.fasta aligned_bam/*_sorted.bam | \
    bcftools call -mv -Oz -o variants.vcf.gz

# Index VCF
bcftools index variants.vcf.gz

# Call variants (from previous alignment)
bcftools mpileup -f ref_genome.fasta aligned_bam/*.bam | bcftools call -mv -Oz -o all_samples.vcf.gz

# Filter variants
bcftools filter -e 'QUAL<30 || DP<10' all_samples.vcf.gz -Oz -o filtered.vcf.gz

# Convert to PHYLIP
vcftools --vcf filtered.vcf.gz --out snp_matrix --phylip

# Build SNP tree
fasttree -nt -gtr snp_matrix.phylip > snp_tree.nwk

# Detect Housekeeping Genes & Resistance Genes
# Create BLAST database
makeblastdb -in reference/ref_genome.fasta -dbtype nucl -out blast_results/ref_db

# Create query files
cat > blast_resistance_genes.fa << 'EOF'
>CYP51A
ATGGCTGTATCAGTTCGTTTCTGCTGCTCG...
>SdhB
ATGGACGTGAAGTGGAGGTCAGGATTTTGAG...
EOF
>
# Extract sequences

/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/makeblastdb -in 100contigs.fasta -dbtype nucl
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query NeopestalotroseITS.fasta -db 100contigs.fasta -outfmt 6 -out 100contigITS.tsv
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools faidx 100contigs.fasta
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools faidx 100contigs.fasta NODE_499_length_3025_cov_4728.818860:450-995 > 100cintigITS_raw.fasta
cat 

# Run BLAST searches
blastn -query blast_resistance_genes.fa -db blast_results/ref_db \
       -out blast_results/resistance_blast.txt -evalue 1e-5 -outfmt 6

# Use blastn to search genes
makeblastdb -in ref_genome.fasta -dbtype nucl -out ref_db
blastn -query housekeeping.fasta -db ref_db -out hk_blast.txt -evalue 1e-5 -outfmt 6
blastn -query fungicide_genes.fasta -db ref_db -out resistance_blast.txt -evalue 1e-5 -outfmt 6
# House keeping genes
>ACT1_Fungal_Actin
ATGGTCGGTATGGGTCAAAAAGAAATGGTAAGGTTATGTCACCAACTGGGACGATATGGAGAAGATCTGGCATCACACCTTCTACAACGAGCTGGAA

>EF1A_Translation_Elongation_Factor_1_Alpha
ATGGGTAAGGAGGACAAGACTGGAAGTTGAAGAGGAGCTTGATTTCAAGAACATGACGCTGAAGGTTGACGATATTGATGGAATCGTCGGAGAGGACA

>GAPDH_Glyceraldehyde-3-phosphate_dehydrogenase
ATGGTTGTTGAGGTTGTTGTTGACGTTGGTGTTGGTGATGCTGTTGTTGCTGATGTTGTTGTTGTTGGTTCTGATGCTGTTGGTGTTGGTGGTCGAG

>TUB2_Beta_Tubulin
ATGCGTGAGTGCATCAAGAAGGACTACCTGGGCACCCACCTCGGAGCGTGTGTACAGTCCACCGCCGATGTCTACGTCGTACAGTTGGGAGCCG

>RPB2_RNA_polymerase_II_largest_subunit
ATGGGCAAGGGTCAGATTGATGTCTGGGTCGATTGTTGAGGCGAGTTCGGTGAGTTCGAGGCTGACTTCGAGGAGTTGATGACGAGCTTGAGG


# Fungicide ressistance genes 

>CYP51A_Sterol_14α-demethylase_Azole_target
ATGGCTGTATCAGTTCGTTTCTGCTGCTCGATGTGTCTCGTGCTGGTGCTGTTGGCGGCTCGTTGTCTCGTTGGAGCTGCTGTTGTTGCTGGT

>SdhB_Succinate_Dehydrogenase_Subunit_B_SDHI_target
ATGGACGTGAAGTGGAGGTCAGGATTTTGAGGTCGTTGAGCTGGTTCTGTTCTGGTGGTGGTGATGTGGTGTCTGCTGATGGTCAGGATGTGGAC

>FKS1_1,3-beta-D-glucan_synthase_Echinocandin_target
ATGGTTGGTGCTGGTTGCTGGTGGTTGCTGATGCTGCTGTTGCTGCTGCTGGTCAGGGTTCTGGTGTTGTTGTTGGTTCTGTTGCTGATGTT

>ABC1_ABC_Transporter_Multidrug_Resistance
ATGGCGCTGGTTGTTGGTTCTGATGTTGTTGATGCTGGTGTTGCTGCTGTTGCTGGTGCTGCTGGTGCTGGTGGTGCTGGTGCTGCTGGTTCT

>MDR1_Major_Facilitator_Superfamily_MDR
ATGGTGCTGGTTGTTGCTGCTGCTGGTGTTGTTGTTGCTGGTCAGGTTGGTGCTGGTCAGGGTGTTGGTGGTCAGGGTCAGGGTGGTCAG

# 1. FKS1 – β-1,3-glucan synthase (echinocandin resistance)
wget -O FKS1.fasta "https://www.uniprot.org/uniprot/P38631.fasta"

# 2. ERG27 – 3-keto-steroid reductase (fenhexamid resistance)
wget -O ERG27.fasta "https://www.uniprot.org/uniprot/Q6CJC2.fasta"

# 3. sdhB – Succinate dehydrogenase iron-sulfur subunit (SDHI resistance)
wget -O sdhB.fasta "https://www.uniprot.org/uniprot/T1QD65.fasta"

# 4. β‑tubulin (TUB2) – benzimidazole (MBC) resistance
wget -O TUB2.fasta "https://www.uniprot.org/uniprot/P40904.fasta"

# 5. Cyt b – Mitochondrial cytochrome b (QoI / azoxystrobin resistance)
wget -O cytb.fasta "https://www.uniprot.org/uniprot/Q7YF16.fasta"

# 6. CYP51 (also known as ERG11) – Azole/DMI fungicide resistance
wget -O CYP51.fasta "https://www.uniprot.org/uniprot/Q4WNT5.fasta"

# 7. HOG1 – Protein kinase (MAPK pathway, stress signaling)
wget -O HOG1.fasta "https://www.uniprot.org/uniprot/A0A8H6BVP0.fasta"



  

    
    # Look for EFF field (SnpEff annotation)
    for (i in info) {
        if (info[i] ~ /^EFF=/) {
            # Check if it's in a resistance gene
            gene_found = 0;
            mutation_found = 0;
            
            for (gene in resistance_genes) {
                if (tolower(info[i]) ~ tolower(gene)) {
                    gene_found = 1;
                    break;
                }
            }
            
            # Check for known mutations
            for (mutation in known_mutations) {
                if (info[i] ~ mutation) {
                    mutation_found = 1;
                    break;
                }
            }
            
            # Print if in resistance gene OR is a known mutation
            if (gene_found || mutation_found) {
                print $0;
            }
        }
    }
}


# Population Genomics

Population structure and genetic diversity: step-by-step

Strawberry crown isolates (NC series), Illumina PE reads, SPAdes assemblies, sample metadata sample_metadatafinal.tsv.

Approach

You have raw reads for every isolate, so the standard route is reference-based SNP calling: map all isolates to one reference, call SNPs jointly with GATK, filter hard, then analyse the SNPs. This is the same logic as the method you pasted (GATK HaplotypeCaller GVCF, joint genotyping, biallelic SNPs, VCFtools filters, repeat masking, adegenet/poppr analyses, index of association, GST/G'ST, clone correction).

The SPAdes assemblies serve only to supply the reference. The Augustus proteins are not needed for SNP-based structure. You could use them later for BUSCO or orthology trees, or to annotate which genes carry the SNPs.

Decision to check: ploidy. 00_config.sh sets PLOIDY=1 (haploid, typical for ascomycete fungi). If your organism is diploid (e.g. an oomycete like Phytophthora, as in the pasted methods), set PLOIDY=2. Everything else adapts automatically.

Sample naming

Read files are named NC100_S423_R1_filtered.fastq.gz, so the VCF sample name is NC100_S423. The metadata Strain_ID is NC100. The R scripts strip the _S### suffix to join the two. Isolates without a metadata row are reported with a warning.

Pipeline
module load bwa/1.0_bec6a80 iqtree/2.1_fb8f580 vcf2phylip/1.1_bf013e8 vcftools/1.0_9326bcf blast/3.0_6bf5a81 bbtools/2.0_48bf7d9 bcftools/1.1_a8ebeaf bamtools/1.0_a3b280e gatk/1.0_5ea5a38 samtools/1.0_5128b90
module load  plink/1.0_9cff8a0 picard/1.0_90fd432 bedtools/1.0_52b61cf snippy/1.1_acdd356
module load conda
conda create -n vcf_tools -c conda-forge -c bioconda htslib bcftools -y
conda activate vcf_tools
source ~/miniconda3/etc/profile.d/conda.sh
conda activate
conda info --base
source /home/muzhinjin/miniconda3/etc/profile.d/conda.sh
conda activate vcf_tools
conda env list
conda run -n vcf_tools bcftools --version
bcftools plugin -l | grep -i fixploidy

#admxture
conda install -c bioconda admixture
export PATH="$PATH:/home/muzhinjin/.conda/envs/admixture/bin"
ls /home/muzhinjin/.conda/envs/admixture/bin/admixture
/home/muzhinjin/.conda/envs/admixture/bin/admixture --help | head -3
env -i PATH=/usr/bin:/bin bash -c 'source ./00_config_v2.sh >/dev/null 2>&1; command -v admixture'

Step	Script	Env	What it does	Check before moving on
0	edit scripts/00_config.sh	–	Paths, ploidy, thresholds	META path is right; PLOIDY correct
1	01_setup_env.sh	–	Creates conda envs popgen, popgen-r	– conda activate popgen-r 
2	02_reference.sh	popgen	Assembly stats, picks reference (highest N50, contigs ≥1 kb), indexes it	Open ref/assembly_stats.tsv. Reference should be a typical genome size for the species, not an outlier. Override with REF_ASM=...
ern jobs submit \
  --name="NeoGATK" \
  --threads=8 \
  --memory=64gb \
  --hours=100 \
  --input="ref/reference.fasta" \
  --module='gatk/1.0_5ea5a38 python=3.6' \
  --command=gatk -- \
  CreateSequenceDictionary \
  -R ref/reference.fasta \
  -O ref/reference.dict
  
2b	02b_repeats_optional.sh	popgen-rep	RepeatModeler/RepeatMasker, writes ref/repeats.bed	Optional but recommended. Step 6 uses the BED automatically if it exists
ern jobs submit \
  --name="NEOreapeatmodeller" \
  --threads=16 \
  --memory=64gb \
  --hours=48 \
  --input="ref/reference.fasta" \
  --module='repeatmodeler/1.0_881fa7a python=3.12' \
  --command=BuildDatabase -- \
  -name refdb \
  reference.fasta
3	03_map.sh	popgen	bwa mem, fixmate, sort, markdup per isolate	logs/bwa_*.log mapping rates
ern jobs submit \
  --name="NeoBWA_Index" \
  --threads=8 \
  --memory=16gb \
  --hours=4 \
  --input="ref/reference.fasta" \
  --module='bwa/1.0' \
  --command=bwa -- \
  index ref/reference.fasta

4	04_bam_qc.sh	popgen	Mapped %, mean depth, breadth; writes samples_pass.txt	Drop low-depth or low-mapping isolates (possible contamination or wrong species)
ern jobs submit  --name=Neo_BAM_QC --threads=4  --memory=32gb  --hours=12  --input="04_bam_qc.sh" --input="00_config.sh"  --input="popgen/samples.txt" --input="popgen/bam"  --module="samtools/2.1_4854cc4 python=3.12" --command="bash 04_bam_qc.sh"

5	05_call_variants.sh	popgen	HaplotypeCaller (-ERC GVCF), CombineGVCFs, GenotypeGVCFs	vcf/raw.vcf.gz exists

ern jobs submit \
  --name="Neo_Variant_Calling" \
  --threads=16 \
  --memory=64gb \
  --hours=100 \
  --input="05_call_variants_v2.sh" \
  --input="00_config.sh" \
  --input="popgen/samples_pass.txt" \
  --input="popgen/bam" \
  --input="popgen/ref" \
  --module="gatk/1.0_5ea5a38" \
  --command="bash 05_call_variants_v2.sh"

6	06_filter_vcf.sh	popgen	Biallelic SNPs, GATK hard filters, DP/GQ/QUAL, repeats, missingness, LD pruning	Read the printed SNP counts after each filter
7	07_admixture.sh	popgen	ADMIXTURE, K = 1–10, 10-fold CV	admix/cv.tsv
8	08_tree_iqtree.sh	popgen	ML SNP tree, +ASC, 1000 ultrafast bootstraps	tree/snp_tree.treefile
10	10_structure.R	popgen-r	PCA, DAPC, UPGMA + 1000 bootstraps, MSN	See outputs below
11	11_diversity.R	popgen-r	MLGs and clone correction, π, θW, Tajima's D, Ia/r̄d, GST, G'ST, FST, AMOVA	See outputs below
12	12_plot_admixture.R	popgen-r	CV curve and ancestry barplot	–
Commands
bash
# 0. edit paths/ploidy
nano scripts/00_config.sh

# 1. once
bash scripts/01_setup_env.sh
conda activate popgen

# 2-8
bash scripts/02_reference.sh        # review ref/assembly_stats.tsv; rerun with REF_ASM=... if needed
# optional: conda activate popgen-rep && bash scripts/02b_repeats_optional.sh && conda activate popgen
bash scripts/03_map.sh
bash scripts/04_bam_qc.sh
bash scripts/05_call_variants.sh
bash scripts/06_filter_vcf.sh

ern jobs submit   --name="Neo_VCF_Filtering_v2"   --threads=16   --memory=64gb   --hours=100   --input="06_filter_vcf_v3.sh"   --input="00_config_v2.sh"   --input="popgen/vcf"   --input="popgen/ref"   --module="gatk/1.0_5ea5a38"   --module="vcftools/1.0_9326bcf"   --module="plink/1.0_9cff8a0"   --module="samtools/1.0_5128b90"   --command="bash 06_filter_vcf_v3.sh"

bash scripts/07_admixture.sh
ern jobs submit   --name="Neo_Admixture"   --threads=16   --memory=32gb   --hours=24   --input="08_admixture_v2.sh"   --input="00_config_v2.sh"   --input="popgen/vcf"   --module="plink/1.0_9cff8a0"   --command="bash 08_admixture_v2.sh"

bash scripts/08_tree_iqtree.sh
 ern jobs submit   --name="Neo_IQTreeSNP"   --threads=16   --memory=32gb   --hours=48   --input="09_iqtree_v3.sh"   --input="00_config_v2.sh"   --input="popgen/vcf"   --module="vcf2phylip/1.1_bf013e8"   --module="iqtree/2.1_fb8f580"   --command="bash 09_iqtree_v3.sh"


# R steps
conda activate popgen-r
source scripts/00_config.sh         # sets WORK, META, PLOIDY (exported)
R_OUT=$WORK/results

Rscript scripts/10_structure.R   $WORK/vcf/structure.vcf.gz   $META $R_OUT/structure
Rscript scripts/11_diversity.R   $WORK/vcf/diversity.vcf.gz   $META $R_OUT/diversity $WORK/ref/callable_bp.txt
Rscript scripts/12_plot_admixture.R $WORK/admix $META $R_OUT/admixture      # optional args: K  facet_variable

For long runs on a server, wrap steps 3 and 5 in nohup ... & or submit through SLURM. Both are resumable, because finished BAMs and GVCFs are skipped.

What the results mean
PCA, DAPC, ADMIXTURE and the trees should agree. If they do, the structure is real. DAPC_cluster_vs_metadata.txt cross-tabulates clusters against County, Region, Cultivar, Host, Plant part and collection month. Clusters that track a metadata variable point to geographic, host or temporal structure. Clusters spread across all of them suggest a single, mixed population.
ADMIXTURE K choice: take the lowest CV error, or the elbow if the curve is flat. Plot neighbouring K values (12_plot_admixture.R ... <K>). K = 1 being best means no structure.
MLGs (MLG_assignments.csv): many isolates collapsing into few MLGs means a clonal population, as expected for a strawberry pathogen spread on transplants. Check pairwise_distance_hist.pdf. Set CLONE_THR at the gap between near-identical pairs and distinct pairs. The default 0.01 follows your pasted methods. Rerun step 11 with CLONE_THR=0.005 Rscript ... to test sensitivity.
Diversity (diversity_by_group.csv): π and θW are per site, with one row for all isolates and one per group, computed on both the full and the clone-corrected data. Tajima's D < 0 suggests expansion or a recent sweep (typical of a clonal lineage), and D > 0 suggests balancing selection or a structured or admixed population.
Index of association (index_of_association.csv): r̄d near 0 with a non-significant permutation test (ia_permutation_test.txt) is consistent with recombination. r̄d clearly > 0 means linkage, i.e. clonal reproduction.
Differentiation (differentiation.txt): pairwise GST, G'ST and Weir & Cockerham FST, plus AMOVA per metadata variable. Interpret groups with few isolates cautiously, since MIN_N defaults to 5.
Caveats you should know about
π and θ are approximations. They are computed from SNPs only and divided by the "callable" length (reference size minus repeats). This ignores invariant sites that failed depth or quality filters, so per-site values are somewhat underestimated. For publication-grade π, dxy and D, call all sites (GenotypeGVCFs --all-sites) and use pixy. The environment recipe is in 01_setup_env.sh.
Metadata balance. Most isolates share Host = Strawberry and Plant part = Crown. Those columns, and any group with fewer than MIN_N isolates, are skipped automatically. County, Cultivar, and collection date are the informative variables. If many cultivars are NA, they form an "Unknown" group that is excluded.
Reference bias. Isolates far from the reference map poorly. Check mapping rates in bam_qc.tsv. If one isolate or clade has a much lower rate, repeat with a reference from that clade.
Not tested on your data. I could not run these scripts in my environment (no R, no network), so I only checked the shell scripts' syntax. Expect to fix small package-version issues on first run, especially in the R steps. Send me the error messages and I'll patch them.

#Tribase

mkdir -p extracted_loci/refs
for p in ITS:KY271740 TEF:KY271093 TUB2:KY271094; do
  curl -sSf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${p#*:}&rettype=fasta&retmode=text" > extracted_loci/refs/${p%%:*}.fasta
done
head -c 200 extracted_loci/refs/*.fasta   # check they look like FASTA


mafft --auto --thread 8 TUB2_trimmed.fasta > TUB2.aln.fasta
 python3 concat_loci.py Neopestalotiopsis  ITS.aln.fasta TEF.aln.fasta TUB2.aln.fasta

