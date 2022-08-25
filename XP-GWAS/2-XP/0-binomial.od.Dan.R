#### function to do binomial test ####
binomial.od <- function(data, cat1cols, cat2cols, geno, lambda = 1) {
  #library(dispmod)
  cat1counts <- data[cat1cols]
  cat1counts <- as.numeric(cat1counts)
  cat2counts <- data[cat2cols]
  cat2counts <- as.numeric(cat2counts)
  geno <- as.character(geno)
  # check the levels of geno
  geno.level <- unique(geno)
  if (length(geno.level) == 2) {
    geno <- as.factor(geno)
    full <- glm(cbind(cat1counts, cat2counts) ~ geno,
             family = binomial(link=logit)) ## Binomial test
    logm2b <- full[[1]][2] # log odd ratio
    log2m2b <- log2(exp(logm2b))
    reduced <- glm(cbind(cat1counts, cat2counts) ~ 1,
             family = binomial(link=logit)) ## Binomial test

	###
	### method1: quasi bionomial test, developed by Dan Nettleton
	###
	lrt <- anova(reduced, full, test="LRT")
	
	### Adjust for overdispersion using a quasilikelihood approach
	dfr <- df.residual(reduced)
	dff <- df.residual(full)
	devr <- deviance(reduced)
	devf <- deviance(full)
	
	### Check for overdispersion
	overdispersion.pval <- 1 - pchisq(devf, dff)
	od <- 1
	od.new <- devf/dff
	if (overdispersion.pval < 0.05) {
		od <- max(1, od.new)
	}

	Fstat <- ((devr-devf)/(dfr-dff))/od/lambda
	cis.pvalue <- 1 - pf(Fstat, dfr-dff, dff)
  } else {
    cis.pvalue <- NA
	log2m2b <- NA
	od <- NA
  }
  # output variables:
  log2m2b <- round(log2m2b, 3)
  #cis.pvalue <- format(cis.pvalue, digits = 4)
  c(data, log2m2b, od, cis.pvalue)
}

