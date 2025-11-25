module load cluster/hpc
ern cluster create --name=hpc --head-node=login.hpc.ufs.ac.za --node-select=ern --email-address=MuzhinjiN@ufs.ac.za --shared-fs
#Loop so many samples
ern jobs submit   --name=Neosamtoolsview   --threads=32   --memory=128gb   --hours=48   --input="*.sam"   --module="samtools/2.1_4854cc4"   --command='
    for f in *.sam; do
        base=${f%.sam}
        samtools view -@ 8 -bS "$f" > "${base}_unsorted.bam"
    done

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
