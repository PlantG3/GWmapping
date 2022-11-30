
## BWA alignment ##
bwadb=reference.fa
reads_R1=sample_trim.R1.fq.gz
reads_R2=sample_trim.R2.fq.gz
bwa mem -t 8 -R '@RG\tID:B73\tSM:B73' $bwadb $reads_R1 $reads_R2 > sample_align.sam

## SAM filter ##
perl samparser.bwa.pl -i sample_align.sam -e 100 -m 3 100 --tail 5 100 --gap 0 --insert 100 800 1>sample.parse.sam 2>sample.parse.log
samtools view -u sample.parse.sam | samtools sort -o sample.parse.sort.bam
samtools index sample.parse.sort.bam