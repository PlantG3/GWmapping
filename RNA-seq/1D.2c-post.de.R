setwd("/bulk/liu3zhen/research/projects/GW/main/1_RNA-Seq/1D_v3DE/")

### annotation
annot <- read.delim("~/references/B73Ref3/genemodels/gene_v3_annotation_Cheng.txt")

### DE result:
rr <- read.delim("./R.CMN_CK/CMN_CK.DESeq2.txt", stringsAsFactors = F)
sr <- read.delim("./S.CMN_CK/CMN_CK.DESeq2.txt", stringsAsFactors = F)
rr2 <- rr[, c("GeneID", "CMN_CK.log2FC", "CMN_CK.qval")]
colnames(rr2) <- c("GeneID", "R.log2FC", "R.qval")
sr2 <- sr[, c("GeneID", "CMN_CK.log2FC", "CMN_CK.qval")]
colnames(sr2) <- c("GeneID", "S.log2FC", "S.qval")
rs <- merge(rr2, sr2, by="GeneID", all = T)
head(rs)
### output
rs2 <- merge(rs, annot, by = "GeneID")
write.table(rs2, "1D.2o-CMN_DE_result.txt", quote = F, row.names = F, sep = "\t")

