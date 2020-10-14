### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
library(gamlss)
source("F00.00.generic.R")
source("F01.01.base.R")
source("F01.01.goodness_of_fit.R")
args = commandArgs(trailingOnly=TRUE)  # passed from script

if (is.na(args[1])|args[1] == "") {
  zoe = 1
  nrm = "tpm5"
} else {
  zoe = as.numeric(args[1])
  nrm = as.numeric(args[2])
  nrm = switch(nrm, "1" = "tpm5", "2" = "rpk", "3" = "asin")
}
print(sprintf("zoe %s - %s", zoe, nrm))

### 0.2 Data
# Raw data of 118 subjects
if (zoe == 1) {
  gene.marginal.RPK.DRNA <- readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE1.rds")
} else if (zoe == 2) {
  gene.marginal.RPK.DRNA <- readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE2.rds")
}
excluded.subject <- gene.marginal.RPK.DRNA$meta$id %in% c(352, 420, 10083, 11210, 11259, 11790, 12623)
DataMeta = gene.marginal.RPK.DRNA$meta[!excluded.subject,]
RNA     = gene.marginal.RPK.DRNA$otu[,, 2]
DataRPK  = RNA[,colnames(RNA) %in% DataMeta$id]
n.test   = 300


rm(gene.marginal.RPK.DRNA, excluded.subject, RNA)
gc()

ST    = apply(DataRPK, 2, sum)
mean(ST) # 5,551,718 (ZOE1), 21M for ZOE2
const = if (zoe == 1) 5E+6 else 20E+6
DataTPM <- t(t(DataRPK)/ST) * const

# screening and sampling
prev.filter = DataTPM %>% apply(1, function(x) mean(x>0)) %>% {which(. >= 0.1)}
set.seed(1)
i.sample    = sample(prev.filter, n.test)
i.taxa      = rownames(DataTPM)[i.sample]

# Empty shell
tpm.i = DataTPM[i.sample[1], ]
otu.i = DataRPK[i.sample[1], ]
con = gamlss.control(n.cyc = 100, trace = FALSE)
a1 <- est.beta(otu.i[otu.i>0]/const);   nm.beta  = names(a1)
a2 <- est.gamma(otu.i[otu.i>0]);        nm.gamma = names(a2)
a3 <- est.ln(otu.i[otu.i>0]);           nm.ln    = names(a3)
nm.ks <- c("ks.coef", "ks.pval", "lil.pval")

# full.zinb    <- gamlss(as.integer(tpm.i) ~ 1, family = ZINBI(sigma.link = "log"), control = con)
full.zinb    <- zeroinfl(as.integer(tpm.i) ~ 1, dist = "negbin")

gof.beta =
  matrix(NA, nrow = n.test, ncol = length(a1) + 4, 
         dimnames = list(i.taxa, c(nm.beta, "zero.prop", nm.ks)))
gof.gamma =
  matrix(NA, nrow = n.test, ncol = length(a2) + 4, 
         dimnames = list(i.taxa, c(nm.gamma, "zero.prop", nm.ks)))
gof.ln =
  matrix(NA, nrow = n.test, ncol = length(a3) + 4, 
         dimnames = list(i.taxa, c(nm.ln, "zero.prop", nm.ks)))


for (i in seq_along(i.sample)) {
  # if (i < 101) next
  # update dat.reg with the corresponding otu values.
  j = i.sample[i]
  tpm.i[] <- as.numeric(DataTPM[j, ])
  if (nrm == "asin") tpm.i = asn(tpm.i / const) * const
  otu.i[] <- if (nrm == "rpk") {as.numeric(DataRPK[j, ])} else {tpm.i } #tpm5 or asin
  taxa.i  <- rownames(DataTPM)[j] %>% gsub(" \\(TOTAL\\)", "", .) %>% gsub("\\_", " ", .)
# name.i  <- names(DataRPK)[j]
  # if (mean(dat.reg$composition > 0) < 0.1) next
  
  # fitting & KS tests
  args.gof = list(return.est = TRUE, lilliefors = TRUE, n.lilliefors = 300, ks.pval = TRUE)
  gof.beta [i, c(nm.beta, nm.ks)]  = do.call(ks.empirical, c(list(tpm.i[tpm.i > 0]/const, model = "beta"), args.gof))
  gof.gamma[i, c(nm.gamma, nm.ks)] = do.call(ks.empirical, c(list(otu.i[otu.i > 0], model = "gamma"), args.gof))
  gof.ln[i, c(nm.ln, nm.ks)]       = do.call(ks.empirical, c(list(otu.i[otu.i > 0], model = "ln"), args.gof))
  gof.beta [i, "zero.prop"]= gof.gamma[i, "zero.prop"] = gof.ln[i, "zero.prop"] = mean(otu.i == 0)
  if (i <= 5) full.zinb  <- try(zeroinfl(as.integer(otu.i) ~ 1, dist = "negbin"))
  
  
  ## plots
  if (i <= 5) {
  p1  <- 
    qqplot1(values = tpm.i[tpm.i > 0]/const, ks.pval = gof.beta[i, "ks.pval"],
            pbeta, shape1 = gof.beta[i, "alpha"], shape2 = gof.beta[i, "beta"], title = "Beta")
  p2 <-
    qqplot1(values = otu.i[otu.i > 0], ks.pval = gof.gamma[i, "ks.pval"], 
            pGA, mu = gof.gamma[i, "mu"], sigma = gof.gamma[i, "sig"], title = "Gamma")
  p3 <-
    qqplot1(values = otu.i[otu.i > 0], ks.pval = gof.ln[i, "ks.pval"], 
            pLOGNO, mu = gof.ln[i, "mu"], sigma = gof.ln[i, "sig"], title = "Log-normal")
  p4 <-
    qqplot.zinb(values = as.integer(tpm.i), estimates = param.zinb2(full.zinb))
  
  pp <- gridExtra::grid.arrange(p1, p2, p3, p4, 
                                top = grid::textGrob(paste0("Fitted models for ", taxa.i)))
  ggsave(paste0("figure/C0102QQplot_zoe", zoe, "_", nrm, "_",i, ".png"), plot = pp,
         width = 12, height = 9)
  }
  if (! i %% 10) {cat("i = ", i, 
                      ", mean.beta.ks.p = ", sprintf("%.5f", mean(gof.beta[, "ks.pval"], na.rm = T)), 
                      " mean.gamma.ks.p = ", sprintf("%.5f", mean(gof.gamma[, "ks.pval"], na.rm = T)), 
                      " mean.ln.ks.p = ", sprintf("%.5f", mean(gof.ln[, "ks.pval"], na.rm = T)), "\n")}
  saveRDS(list(beta = gof.beta, gamma = gof.gamma, ln = gof.ln),
          paste0("output/gof_zoe", zoe, "_", nrm, ".tmp.rds"))
}

saveRDS(list(beta = gof.beta, gamma = gof.gamma, ln = gof.ln),
        paste0("output/gof_zoe", zoe, "_", nrm, ".rds"))


if (0) {
  for (zoe in 1:2) {
    gof.total <- NULL
    for (nrm in c("rpk", "tpm5", "asin")) {
      nrm2 = switch(nrm, "rpk" = "RPK", "tpm5" = "TPM", "asin" = "arcsin")
      gof <- readRDS(paste0("output/gof_zoe", zoe, "_", nrm, ".rds"))
      gof.total <- 
        rbind(gof.total,
              gof$beta[, "lil.pval"] %>% data.frame(pval = ., distribution = "Beta", nrm = nrm2),
              gof$ln[, "lil.pval"] %>% data.frame(pval = ., distribution = "Log-normal", nrm = nrm2),
              gof$gamma[, "lil.pval"] %>% data.frame(pval = ., distribution = "Gamma", nrm = nrm2))
      
    }
    gof.stat <- 
      aggregate(pval ~ distribution + nrm, data = gof.total, 
                FUN = function(s) mean(s < 0.05, na.rm = TRUE)) %>% 
      mutate(reject = paste0("%(p < 0.05) = ", signif(pval, 2)), 
             pval = 0.5, count = if (zoe == 1) 20 else 100,
             reject = ifelse(reject == "%(p < 0.05) = 0", "%(p < 0.05) < 0.001", reject))
    gof.total %>% 
      mutate(nrm = factor(nrm, levels = c("RPK", "TPM", "arcsin"))) %>% 
      ggplot(aes(pval)) + 
      facet_grid(distribution ~ nrm) + 
      geom_histogram(binwidth = 0.01) + 
      # ggtitle("Kolmogorov-Smirnov test p-value histogram for Beta, Log-normal and Gamma distribution") +
      xlab("KS (Lilliefors) test p-values") + ggtitle(if (zoe == 1) "(A) ZOE 1.0 (n = 116)" else "(A) ZOE 2.0 (n = 297)") + 
      geom_vline(xintercept = 0.05, col = "red") + 
      geom_text(data = gof.stat, aes(pval, count, label = reject)) +
      theme_bw()
      # annotate("text", x = 0.75, y = 25, label = paste0("%(p < 0.05) = ", mean(gof.gamma[, "ks.pval"] < 0.05)))
    ggsave(paste0("figure/C0102KS_zoe",zoe, "-p-Histogram_.png"))
  }
}
