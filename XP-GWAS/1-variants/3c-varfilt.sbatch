#!/bin/bash -l
#SBATCH --mem-per-cpu=60G
#SBATCH --time=0-23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=ksu-gen-highmem.q
invcf=CmnRS.2.select.vcf
outvcf=CmnRS.3.filt
outtxt=CmnRS.4.filt.AC.txt
#vcftools --vcf $invcf --maf 0.15 --max-missing 0.5 --recode --out $outvcf
perl ~/scripts/vcf/vcf2AC.pl $outvcf.recode.vcf > $outtxt

