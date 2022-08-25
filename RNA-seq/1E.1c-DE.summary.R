setwd("/bulk/liu3zhen/research/projects/GW/main/1_RNA-Seq/1E_postDE")
options(stringsAsFactors = F)

###############################################################################
# parameters
###############################################################################
fdr.cutoff <- 0.05
log2fc.cutoff <- 0.6

###############################################################################
# load data
###############################################################################
de <- read.delim("../1D_v3DE/1D.2o-CMN_DE_result.txt")
head(de)
###############################################################################
# DE
###############################################################################
rde <- !is.na(de$R.qval) & de$R.qval<fdr.cutoff & abs(de$R.log2FC)>log2fc.cutoff
sum(rde)

sde <- !is.na(de$S.qval) & de$S.qval<fdr.cutoff & abs(de$S.log2FC)>log2fc.cutoff
sum(sde)

rsde <- rde & sde
sum(rsde)

rsde_same_direction <- !is.na(de$S.qval) & !is.na(de$R.qval) &
   ((de$S.qval<fdr.cutoff & de$S.log2FC>log2fc.cutoff & de$R.qval<fdr.cutoff & de$R.log2FC>log2fc.cutoff) |
   (de$S.qval<fdr.cutoff & de$S.log2FC< -1*log2fc.cutoff & de$R.qval<fdr.cutoff & de$R.log2FC< -1*log2fc.cutoff))
sum(rsde_same_direction)

rsde_diff_direction <- !is.na(de$S.qval) & !is.na(de$R.qval) &
        ((de$S.qval<fdr.cutoff & de$S.log2FC>log2fc.cutoff & de$R.qval<fdr.cutoff & de$R.log2FC< -1*log2fc.cutoff) |
        (de$S.qval<fdr.cutoff & de$S.log2FC< -1*log2fc.cutoff & de$R.qval<fdr.cutoff & de$R.log2FC>log2fc.cutoff))

de[rsde_diff_direction, ]

###############################################################################
# assign log2FC with NA to random number between -0.6 to 0.6
###############################################################################
de2 <- de[, 1:5]
Rna <- de2$S.qval<fdr.cutoff & is.na(de2$R.log2FC)
nRna <- sum(Rna)
de2$R.log2FC[Rna] <- runif(nRna, min=-1*log2fc.cutoff, max=log2fc.cutoff)
Sna <- de2$R.qval<fdr.cutoff & is.na(de2$S.log2FC)
nSna <- sum(Sna)
de2$S.log2FC[Sna] <- runif(nSna, min=-1*log2fc.cutoff, max=log2fc.cutoff)


###############################################################################
# scatter plot
###############################################################################
pdf("1E.1o-DE.scatter.plot.pdf", width=4.5, height=4.5)
par(mfrow=c(1,1), mgp=c(2, 0.7, 0))
plot(de2$R.log2FC, de2$S.log2FC, xlim=c(-4, 10), ylim=c(-4, 10),
     xlab="log2FC in R lines", ylab="log2FC in S lines",
     col="gray75", cex=0.5, las=1,
     main="DE upon bacterial infection")

points(de2$R.log2FC[rde & ! sde], de2$S.log2FC[rde & ! sde], col="tomato3", cex=0.5)
points(de2$R.log2FC[sde & ! rde], de2$S.log2FC[sde & ! rde], col="dodgerblue3", cex=0.5)
points(de2$R.log2FC[rsde], de2$S.log2FC[rsde], col="olivedrab4", cex=0.5)

abline(a=0, b=1, col="red", lty=2)
abline(h=0, v=0, col="orange", lty=2)

legend("bottomright", col = c("tomato3", "dodgerblue3", "olivedrab4"),
       legend = c(paste0("DE only in R, N=", sum(rde & ! sde)),
                  paste0("DE only in S, N=", sum(sde & ! rde)),
                  paste0("DE in both, N=", sum(rsde))),
       bty = "n",  pch = 19, cex=0.8)
dev.off()

###############################################################################
# R DE only output
###############################################################################
Ronly <- de[rde & ! sde, ]
Ronly <- Ronly[order(Ronly$R.qval), ]
write.table(Ronly, "1E.1o-Ronly.DE.txt", quote=F, row.names=F, sep="\t")

Ronly.large.diff <- Ronly[(!is.na(Ronly$R.log2FC) & is.na(Ronly$S.log2FC)) |
                        (Ronly$R.log2FC>0 & Ronly$R.log2FC - Ronly$S.log2FC>1) |
                        (Ronly$R.log2FC<0 & Ronly$R.log2FC - Ronly$S.log2FC< -1), ]
nrow(Ronly.large.diff)
write.table(Ronly.large.diff, "1E.1o-Ronly.log2FC.diff.gt1.DE.txt", quote=F, row.names=F, sep="\t")

###############################################################################
# average log2FC of R and S 
###############################################################################
mean.de <- de[1]
mean.de$log2FC <- apply(de[, c(2,4)], 1, mean, na.rm=T)
mean.de$qval <- apply(de[, c(3,5)], 1, mean, na.rm=T)

write.table(mean.de, "1E.1o-DE.RSmean.txt", quote=F, row.names=F, sep="\t")
