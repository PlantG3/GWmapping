infile <- commandArgs(trailingOnly=TRUE)
source("0-binomial.od.Dan.R")
#infile <- "../1-variants/test.txt"
outfile <- gsub(".*\\/", "", infile)
outfile <- paste0("binomial.", outfile)

### input
d <- read.delim(infile)
ncols <- ncol(d)

# convert NA to 0
for (i in 6:ncols) {
	d[is.na(d[, i]), i] <- 0
}

nsamples <- (ncols - 5) / 2 # number of samples
pheno <- substr(grep("_REF", colnames(d), value = T), 0 ,1)  # phenotype data

### test
dx2 <- apply(d[, 6:ncols], 1, binomial.od,
             cat1cols = 2*(1:nsamples) - 1,
             cat2cols = 2*(1:nsamples),
             geno = pheno,
             lambda = 1)

dx2t <- t(dx2)
colnames(dx2t)[(ncol(dx2t)-2):ncol(dx2t)] <- c("log2OR", "overdispersion", "pvalue")

dm <- cbind(d[, 1:4], dx2t)
head(dm)
### output
write.table(dm, outfile, quote = F, row.names = F, sep = "\t")

which.min(dm$pvalue)
dm[21, ]
