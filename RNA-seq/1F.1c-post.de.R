setwd("/bulk/liu3zhen/research/projects/GW/main/1_RNA-Seq/1F_pathway")
options(stringsAsFactors = F)

###############################################################################
# parameters
###############################################################################
fdr.cutoff <- 0.05
log2fc.cutoff <- 1

###############################################################################
# load data
###############################################################################
meanDE <- read.delim("../1E_postDE/1E.1o-DE.RSmean.txt")
load("~/scripts2/mapmanplus/data/B73_5b.rda")

###############################################################################
# mapman
###############################################################################
source("~/scripts2/mapmanplus/R/mapman.format.from.DE.R")
mapman.format(data = meanDE,
             geneid.colname = "GeneID",
             log2fc.namefeature = "log2FC",
             min.abs.log2fc = log2fc.cutoff,
             fdrcol.namefeature = "qval",
             fdr=fdr.cutoff, output.path = ".",
             output.file = "meanRS.DE.result.mapman.txt",
             primary.transcript=b73_5b.primary.transcript)

### data
source("~/scripts2/mapmanplus/R/mapman.plot.R")
de <- read.delim("meanRS.DE.result.mapman.txt")
head(de)
aliasdb <- read.delim("~/references/maizeGDB/B73v3.cc.genes.v0.1.txt")
aliasdb <- aliasdb[, 1:2]
nrow(aliasdb)
aliasdb <- aliasdb[!duplicated(aliasdb$Gene.Model.ID), ]

# plot select hormone families - I:
hormonekeys <- c("hormone metabolism.abscisic acid",
                 "hormone metabolism.auxin",
                 "hormone metabolism.ethylene",
                 "hormone metabolism.jasmonate",
                 "hormone metabolism.gibberelin")
hormonetitles <- c("ABA", "auxin", "ethylene", "jasmonate", "gibberelin")
mapman.plot(mapmankey=hormonekeys, db=b73_5b.mapman.db, de=de, pdfout=T,
            alias=T, pwidth = 8, outpath = ".", outfile = "1F.1o-Cmn.hormone.1.pdf",
            aliasDB=aliasdb, title=hormonetitles, log2fc.max=3, text.font=0.7)

# plot select hormone families - II:
hormonekeys2 <- c("hormone metabolism.brassinosteroid",
                 "hormone metabolism.cytokinin",
                 "hormone metabolism.salicylic acid.synthesis-degradation")
hormonetitles2 <- c("brassinosteroid", "cytokinin", "SA")
mapman.plot(mapmankey=hormonekeys2, db=b73_5b.mapman.db, de=de, pdfout=T,
            alias=T, pwidth = 7.5, outpath = ".", outfile = "1F.1o-Cmn.hormone.2.pdf",
            aliasDB=aliasdb, title=hormonetitles2, log2fc.max=3, text.font=1)


##################################################################################
# TF
##################################################################################
# plot select TF families:
tfkeys <- c("RNA.regulation of transcription.AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family",
          "RNA.regulation of transcription.MYB domain transcription factor family",
          "RNA.regulation of transcription.WRKY domain transcription factor family",
          "RNA.regulation of transcription.bHLH,Basic Helix-Loop-Helix family",
          "RNA.regulation of transcription.C2H2 zinc finger family",
          "RNA.regulation of transcription.NAC domain transcription factor family",
          "RNA.regulation of transcription.HB,Homeobox transcription factor family")
titles <- c("EREBP", "MYB", "WRKY", "bHLH", "C2H2", "NAC", "Homeobox")
mapman.plot(mapmankey=tfkeys, db=b73_5b.mapman.db, de=de, pdfout=T, alias=T,
            pwidth = 12, outpath = ".", outfile = "1F.1o-Cmn.TF.pdf",
            aliasDB=aliasdb, title=titles, log2fc.max=3)


