#!/bin/bash
module load Java
/bulk/hecheng/software/gwas/tassel5/run_pipeline.pl -Xmx18g -fork1 -importGuess ../maize269_mis10_maf5_10000snp_sort.hmp.txt -PrincipalComponentsPlugin -covariance true -endPlugin -export pca_tassel -runfork1
/bulk/hecheng/software/gwas/tassel5/run_pipeline.pl -Xmx18g -importGuess ../maize269_mis10_maf5_10000snp_sort.hmp.txt -KinshipPlugin -method Centered_IBS -endPlugin -export kins_tassel.txt -exportType SqrMatrix
