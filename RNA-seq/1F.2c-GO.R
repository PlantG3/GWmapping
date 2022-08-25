setwd("/bulk/liu3zhen/research/projects/GW/main/1_RNA-Seq/1F_pathway")
source("~/scripts/RNA-Seq/DE.summary.R")
source("~/scripts/RNA-Seq/RNA_seq.plot.R")
source("~/scripts/go/goseq2019.R")
source("~/scripts/go/goplot.R")

godb <- read.delim("/homes/liu3zhen/references/B73Ref3/genemodels/B73Ref3.gene2go.txt")

# check and create the directory
if (dir.exists("GOenrichments")) {
	system("rm -rf GOenrichments")
}
system("mkdir GOenrichments")

###############################################################################
# parameters
###############################################################################
fdr.cutoff <- 0.05
log2fc.cutoff <- 0.6

###############################################################################
# load data
###############################################################################
readcounts <- read.delim("../1C_readcounts/1C.1o-raw.read.counts")
de <- read.delim("../1D_v3DE/1D.2o-CMN_DE_result.txt")
rde <- !is.na(de$R.qval) & de$R.qval<fdr.cutoff & abs(de$R.log2FC)>log2fc.cutoff
sde <- !is.na(de$S.qval) & de$S.qval<fdr.cutoff & abs(de$S.log2FC)>log2fc.cutoff
rsde <- rde & sde
sum(rsde)

nonde <- !rde & !sde

input.de <- de[rsde | nonde, ]

# sig col
sigcolname <- "cm_ck.sig"
de[, sigcolname] <- "no"
de[rsde, sigcolname] <- "yes"

# log2fc
fc_colname <- "cm_ck.log2fc"
de[, fc_colname] <- apply(de[, c(2,4)], 1, mean, na.rm=T)
nrow(de)
# subset
input.de <- de[rsde | nonde, c("GeneID", sigcolname, fc_colname)]
nrow(input.de)
input.de2 <- merge(input.de, readcounts, by.x="GeneID", by.y="Gene")
head(input.de2)

# GO-Seq
gores1 <- goseq.auto(data=input.de2, godb=godb, geneheader="GeneID", sigcolname=sigcolname,
                     nsampling=100000, rawdatacol = 3:ncol(input.de2), 
                     log2colname=fc_colname, padjmethod="BY", qvalcutoff=1,
                     up.down="up", pvalcutoff=0.01, outpath = "GOenrichments")

gores2 <- goseq.auto(data=input.de2, godb=godb, geneheader="GeneID", sigcolname=sigcolname,
                     nsampling=100000, rawdatacol=3:ncol(input.de2),
                     log2colname=fc_colname,  padjmethod="BY", qvalcutoff=11,
                     up.down="down", pvalcutoff=0.01, outpath = "GOenrichments")

# GO plot
gores1$Term <- gsub(",.*", "", gores1$Term)
goplot(go.out = gores1, main.space = 2, order.by = "pvals", 
       term.space=16, pval.cutoff = 0.001, xlim = c(0, 800),
       outsave=T, outfile="1F.2o-cm_ck.up.GO.pdf",
       pdf.width=6, pdf.height=5,
       main="GO enrichment of UP genes")

gores2$Term <- gsub(",.*", "", gores2$Term)
goplot(go.out = gores2, main.space = 2, order.by = "pvals", 
       term.space=14, pval.cutoff = 0.0001, xlim = c(0, 220),
       outsave=T, outfile="1F.2o-cm_ck.dn.GO.pdf",
       pdf.width=6, pdf.height=5,
       main="GO enrichment of DOWN genes")
#############################################






