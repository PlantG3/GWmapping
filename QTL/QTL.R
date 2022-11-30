library(qtl)
indata <- read.cross("csvs", ".", genfile = "IBMRIL.geno.map.csv",
                     phefile = "IBMRIL13D.pheno.point2.csv",na.strings = "-")
indata <- jittermap(indata,amount = 1e-6)
###scanno###
cal.data <- calc.genoprob(indata)
out.em <- scanone(cal.data,model = "normal",method = "em")
operm.em <- scanone(cal.data,method = "em",n.perm = 1000)
cutoff.em <- summary(operm.em,alpha=c(0.05,0.1,0.2))
write.table(cutoff.em,"permutation.cutoff.txt",quote = F)
write.table(out.em,"20200813_qtl_IBM.txt",row.names = T,col.names = T,sep = "\t",quote = FALSE)
qtl.em <- summary(out.em, perms = operm.em, alpha = 0.2,pvalues = TRUE)
write.table(qtl.em,"20200813_sig_qtl_IBM.txt",row.names = T,col.names = T,sep = "\t")
###make qtl###
out.em <- scanone(cal.data,model = "normal",method = "em")
mqtl <- makeqtl(cal.data,chr=qtl.em[,1],pos=qtl.em[,2],what = "prob")

### interaction between QTLs ##
intqtl <- addint(cal.data,pheno.col = 1,qtl = mqtl,formula = y~Q1+Q2,method = "hk",model = "normal",pvalues = T)
summary(intqtl) #### the pvalue of these three QTL were all > 0.05, means there is no interaction between them.

## estimate the pve (percent of variance ) per QTL
fqtl <- fitqtl(cal.data,dropone = T,get.ests = T,model = "normal",qtl= mqtl,method = "hk",formula = y~Q1+Q2)
summary(fqtl)
