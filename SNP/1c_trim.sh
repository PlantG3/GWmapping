#!/bin/bash
## Reads trimming ##
java -Xmx16g -jar trimmomatic-0.38.jar PE \
sample.R1.fq.gz sample.R2.fq.gz \
sample_trim.R1.fq.gz sample_unpaired.R1.fq.gz \
sample_trim.R2.fq.gz sample_unpaired.R2.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:3:20:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:13 MINLEN:40 2