#!/bin/bash
rbam=/bulk/liu3zhen/research/projects/CMN/RSwgs/1-variants/0-bam/1-R
sbam=/bulk/liu3zhen/research/projects/CMN/RSwgs/1-variants/0-bam/2-S
bamdir=$rbam,$sbam
#bamdir=/bulk/hecheng/cheng_project_beocat/maize_cmn_WGS/R_S_set/R_S_bam/R_bam,/bulk/hecheng/cheng_project_beocat/maize_cmn_WGS/R_S_set/R_S_bam/S_bam
ref=/homes/liu3zhen/references/B73Ref3/GATK/Zea_mays.AGPv3.23.dna.genome.fa
perl ~/local/slurm/snp/gatk.sbatch.pl \
  --outbase CmnRS \
  --bampaths $bamdir \
  --ref $ref \
  --mem 24G --maxlen 2000000
  #--checkscript

# change output name
#mv CmnRS_FriJun212244302019/ CmnRSvcf

