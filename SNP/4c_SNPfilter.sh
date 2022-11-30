#!/bin/bash
## select bi-allelic SNPs ##
vcf=sample.raw.0.vcf
ref=reference.fasta
gatk SelectVariants \
        -R $ref \
        -V $vcf \
        --restrict-alleles-to BIALLELIC \
        -select-type SNP \
        -O sample.bi.1.vcf.gz

## Hard-filtering germline short variants ##
gatk VariantFiltration --java-options '-Xmx24g' -R $ref -V sample.bi.1.vcf.gz \
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filter-name "hard_filter"  \
-O sample.HF.2.vcf.gz

## Extract PASS SNPs ##
gatk SelectVariants --java-options '-Xmx24g' -R $ref -V sample.HF.2.vcf.gz \
--exclude-filtered \
-O sample.PASS.3.vcf.gz

## (Optional) Convert heterozygous SNPs to Missing ##
gatk VariantFiltration --java-options '-Xmx24g' \
-V sample.PASS.3.vcf.gz \
-O sample.mark.hetero.4.vcf.gz \
--genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHetFilter"

gatk SelectVariants --java-options '-Xmx24g' SelectVariants \
-V sample.mark.hetero.4.vcf.gz \
--set-filtered-gt-to-nocall \
-O sample.het2miss.4.vcf.gz

## (Optional) Filter SNPs by Minor Allele Frequency (MAF) & Missing Rate (MR) ##
vcftools .sample.het2miss.4.vcf.gz \
--maf 0.02 --max-missing 0.2 --recode --recode-INFO-all \
--out sample.MAF.MR.5.vcf.gz