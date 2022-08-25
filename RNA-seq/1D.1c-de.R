#############################################
### Project: CMN RNASEQ 2018 BATCH1&2
### RNA-Seq
### Kansas State University
### Sanzhen Liu
### date: 11/7/2020
#############################################

#############################################
### load modules
#############################################
source("~/scripts/RNA-Seq/DESeq/DESeq2.single.trt.R")
source("~/scripts/RNA-Seq/DESeq/DESeq2.block.R")
source("~/scripts/RNA-Seq/DE.summary.R")
source("~/scripts/RNA-Seq/normalization.R")
#source("~/scripts/RNA-Seq/genecount.merger.R")
source("~/scripts/star/starRCmerger.R")
source("~/scripts/RNA-Seq/RNA_seq.plot.R")
#############################################

#############################################
### Parameters - Subject to change
#############################################
library("DESeq2")
working.path <- "/bulk/liu3zhen/research/projects/GW/main/1_RNA-Seq/1D_v3DE/"
setwd(working.path)

### import read counts:
inpath <- "../1B_star"
allfiles <- list.files(path = inpath, all.files = T)
allfiles <- allfiles[grep("ReadsPerGene", allfiles)]
#allfiles <- allfiles[-grep("[RS][12]CMN|[RS][12]CK", allfiles)]  ### remove two reps
allfiles

groups <- c("R", "S")

### experimental design
expdesign <- read.delim("../../../data/RvsS_RNASeq/expdesign.original.txt")
expdesign

###
### perform pair-wise comparison
###
for (cg in groups) {
  g1 <- "CK"
  g2 <- "CMN"
  common.out.name <- paste0(g2, "_", g1)
  feature <- paste0("^", cg)
  selfiles <- allfiles[grep(feature, allfiles)]
  print(selfiles)

  ### readcounts:
  indata <- starRCmerger(datapath = inpath, selection = feature)
  print(head(indata, 5))
  #colnames(indata) <- gsub("CK", "CMN", colnames(indata))
  print(head(indata, 1))


  ### import library size
  #lib <- read.delim("lib.size.txt", header=F)
  #lib[, 1] <- gsub("-", ".", lib[, 1])
  libsize <- colSums(indata[, -1])
  #libsize <- lib[,2]
  #names(libsize) <- lib[,1]

  ### comparison information
  samples <- colnames(indata)[2:ncol(indata)]
  comparison <- c(g1, g2)

  ### DE parameters
  fdr.cutoff <- 0.05

  ### output directory:
  outpath <- paste0(working.path, "/", cg, ".", common.out.name)
  system(paste("mkdir", outpath))

  de.out.file <- paste(common.out.name, ".DESeq2.txt", sep="")
  de.summary.file <- paste(common.out.name, ".DESeq2.summary.txt", sep="")

  ### report
  logoutfile <- paste(common.out.name, ".log.md", sep="")
  #############################################

  #############################################
  ### testing using diff parameters (DE)
  #############################################
  # paie-wise comparison:
  input <- indata[,c(grep(g1, colnames(indata)),
                  grep(g2, colnames(indata)))]
  rownames(input) <- indata[, 1]

  # DE:
  dcolnames <- colnames(input)
  DE.out <- DESeq2.block(input.matrix = input, block = substr(dcolnames, 0, 2),
                         min.mean.reads = 1, min.positive.samples = 3,
                         group1.col = grep(g1, dcolnames),
                         group2.col = grep(g2, dcolnames),
                         comparison = comparison, geneID = rownames(input),
                         fdr = fdr.cutoff, logpath = outpath, logfile = logoutfile)

  indata.norm <- normalization(counts=indata, methods="RPM", normcolname=samples, libsize=libsize)
  DE.out <- data.frame(DE.out)
  head(DE.out)
  final.out <- merge(indata.norm, DE.out, by=intersect(colnames(indata.norm), colnames(DE.out)))
  nrow(final.out)
  head(final.out)

  # output DE:
  final.out.path.file <- paste(outpath, "/", de.out.file, sep="")
  write.table(final.out, final.out.path.file, sep="\t", quote=F, row.names=F )
  cat("\n**Differential expression result was saved in the file:**  \n",
  file=paste(outpath, "/", logoutfile, sep=""), append=T)
  cat(final.out.path.file, "  \n",
  file=paste(outpath, "/", logoutfile, sep=""), append=T)

  # generate figures
  #comparison_info <- paste(comparison[2], "_", comparison[1], sep="")
  contrast.plot(dedata=final.out, normtext="RPM", contrast=common.out.name,
                outpath=outpath, fdrco=fdr.cutoff, report=logoutfile)
  #############################################

  #############################################
  ### summary in a table:
  #############################################
  de.summary <- DE.summary(DE.path=outpath, DE.files=de.out.file, qval.feature=".qval",
                           log2FC.feature=".log2FC", fdr=fdr.cutoff,
                           out.path=outpath, out.file=de.summary.file)
  #############################################
}

