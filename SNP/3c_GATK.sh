#!/bin/bash
## MarkDuplicates ##
gatk MarkDuplicates --java-options '-Xmx24g' --REMOVE_DUPLICATES true -I sample.parse.sort.bam -M sample.metrics.txt -O sample.RD.bam

## gatk index for reference genome ##
samtool faidx reference.fasta
gatk CreateSequenceDictionary -R reference.fasta -O reference.dict

## GATK SNP Calling  ##
gatk HaplotypeCaller --java-options '-Xmx24g' -R reference.fasta -I sample.RD.bam -O sample.raw.0.vcf

