################################################################################
### Genotype * treatment interaction analysis
### Sanzhen Liu
### 3/5/2018
################################################################################
### required packages
library(DESeq2)
#source("~/scripts/star/starRCmerger.R")

### setup the working directory:
setwd("/bulk/liu3zhen/research/projects/GW/main/1_RNA-Seq/1G_interaction")

### annotation
annot <- read.delim("/homes/liu3zhen/references/B73Ref3/genemodels/gene_v3_annotation_Cheng.txt", stringsAsFactors = F)
head(annot)

### data
### import read counts:
indata <- read.delim("../1C_readcounts/1C.1o-raw.read.counts", stringsAsFactors=F)

### experimental design
expdesign <- read.delim("../../../data/RvsS_RNASeq/expdesign.original.txt", stringsAsFactors = F)
expdesign <- expdesign[order(expdesign$Sample), ]

### subset (8 time points in each condition)
# minimum total reads
minCounts <- 100
d2 <- indata[, c("Gene", expdesign$Sample)]
head(d2)
### require at least "minCounts" total counts from all the samples
d2 <- d2[rowSums(d2[, -1])>=minCounts, ]
nrow(d2)
### sample design information in the DESeq format
sample.info <- data.frame(row.names=expdesign$Sample,
                          geno=as.factor(expdesign$Geno),
                   			  trt=as.factor(expdesign$Inoculation),
                  			  batch=as.factor(expdesign$Batch))

### DESeq analysis
gene <- d2[, 1] ### gene sequence
dds <- DESeqDataSetFromMatrix(countData=as.matrix(d2[, -1]),
                              colData=sample.info,
                              formula(~ geno + trt + batch + geno*trt))

### leave effect, including interaction with "day"
dds <- DESeq(dds, "LRT",
             full=formula(~ geno + trt + batch + geno*trt),
             reduced=formula(~ geno + trt + batch))

### test result
res <- results(dds, cooksCutoff=F,independentFiltering=T, alpha=0.1)

res$gene <- gene
out.coef <- coef(dds, SE=F) ### model coefficients
out.coef <- data.frame(out.coef)
out.coef$gene <- gene
out.se <- coef(dds, SE=T) ### standard errors
out.se <- data.frame(out.se)
out.se$gene <- gene
out.disp <- data.frame(gene=gene, Dispersion=dispersions(dds))

#### filter
sum(!is.na(res$padj))
res2 <- res[!is.na(res$padj), ]
res2 <- data.frame(res2)
nrow(res2)

fdr.cutoff <- 0.05
res2$Sig <- "no"
res2$Sig[res2$padj < fdr.cutoff] <- "yes"

sum(!is.na(res$padj) & res$padj<fdr.cutoff)

sig.genes <- res[!is.na(res$padj) & res$padj<fdr.cutoff, "gene"]

res2annot <- merge(res2[, c("pvalue", "padj", "Sig", "gene")], annot, by.x = "gene", by.y = "GeneID")

sig.out.coef <- out.coef[out.coef$gene %in% sig.genes, ]
sig.out.se <- out.se[out.se$gene %in% sig.genes, ]
sig.out.dispersion <- out.disp[out.disp$gene %in% sig.genes, ]

### the "alpha" value has been adjusted to obtain more a uniform p-value
### histogram shape at the high p-value side
hist(res2$pvalue, xlab="p-values", ylab="Number of genes", main="interaction effect")

### output
write.table(res2annot, "1G.1o-gene.interaction.raw.output.txt", quote=F, row.names=F, sep="\t")
write.table(sig.out.coef, "1G.1o-sig.interaction.genes.coefficients.output.txt",
            quote=F, row.names=F, sep="\t")
write.table(sig.out.se, "1G.1o-sig.interaction.genes.se.output.txt",
            quote=F, row.names=F, sep="\t")
write.table(sig.out.dispersion, "1G.1o-sig.interaction.genes.dispersion.output.txt",
            quote=F, row.names=F, sep="\t")

