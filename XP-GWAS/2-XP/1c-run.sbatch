#!/bin/bash
#SBATCH --array=1-45
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-23:00:00
#SBATCH --partition=ksu-gen-highmem.q
module load R
id=$SLURM_ARRAY_TASK_ID
Rscript 1c-association.split.R /bulk/liu3zhen/research/projects/CMN/RSwgs/1-variants/5-split/ac$id

