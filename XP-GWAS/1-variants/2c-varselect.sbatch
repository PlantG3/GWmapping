#!/bin/bash -l
#SBATCH --mem-per-cpu=60G
#SBATCH --time=0-23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=ksu-gen-highmem.q,batch.q,ksu-biol-ari.q,ksu-plantpath-liu3zhen.q
module load Java/1.8.0_192

# generate a bam list
vcf=./CmnRS/CmnRS.vcf
ref=/homes/liu3zhen/references/B73Ref3/GATK/Zea_mays.AGPv3.23.dna.genome.fa
out=CmnRS.2.select

gatk SelectVariants \
	-R ${ref} \
	-V ${vcf} \
	-select 'DP >= 150.0' \
	--restrict-alleles-to BIALLELIC \
	-select-type SNP \
	-O ${out}.vcf &>${out}.log

