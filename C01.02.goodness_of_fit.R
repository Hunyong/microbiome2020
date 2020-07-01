### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
library(gamlss)
source("F00.00.generic.R")
source("F01.01.base.R")
args = commandArgs(trailingOnly=TRUE)  # passed from script

if (is.na(args[1])|args[1] == "") {
  zoe = 2
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
n.test   = 1000


rm(gene.marginal.RPK.DRNA, excluded.subject, RNA)
gc()

param.beta <- function(beta.obj) {
  mu     = beta.obj$mu.coefficients %>% plogis %>% as.numeric
  sig    = beta.obj$sigma.coefficients %>% exp %>% as.numeric
  aplusb = 1/sig^2 - 1  # see ?BE
  alpha  = mu * aplusb
  beta   = aplusb - alpha
  return(c(alpha = alpha, beta = beta, mu = mu, sig = sig))
}
param.gamma <- function(gamma.obj) {
  mu = gamma.obj$mu.coefficients %>% exp %>% as.numeric
  sig = gamma.obj$sigma.coefficients %>% exp %>% as.numeric
  # alpha = mu / beta
  return(c(mu = mu, sig = sig))
}
param.ln <- function(ln.obj) {
  mu = ln.obj$mu.coefficients %>% as.numeric
  sig = ln.obj$sigma.coefficients %>% exp %>% as.numeric
  return(c(mu = mu, sig = sig))
}
param.zinb <- function(zinb.obj) {
  mu  = zinb.obj$mu.coefficients %>% exp %>% as.numeric
  sig = zinb.obj$sigma.coefficients %>% exp %>% as.numeric
  pi  = zinb.obj$nu.coefficients %>% plogis %>% as.numeric
  return(c(mu = mu, sig = sig, pi = pi))
}

qqplot1 <- function(values, distnFn, ks.pval, ..., title = "Beta", taxa = NULL) {
  values = sort(values)
  pval   = distnFn(values, ...)
  index  = (1:length(values) - 0.5)/length(values)
  kstest = paste0("Kolmogorov-Smirnov p = ", signif(ks.pval, 2))
  taxa   = if (!is.null(taxa)) paste0(" - ", taxa) else NULL
  
  data.frame(index = -log10(index), p = -log10(pval)) %>%
    ggplot(aes(index, p)) + 
    geom_line(col = "slateblue4", size = 1, alpha = 0.5) +
    geom_point() +
    annotate("text", x = 1, y = 0.25, label = kstest) +
    # xlim(c(0,max.q.x)) + ylim(c(0,max.q.y)) +
    xlim(c(0,2)) + ylim(c(0,2)) +
    xlab(expression(paste("Expected (",-log[10], " quantile)"))) + 
    ylab(expression(paste("Observed (",-log[10], " quantile)"))) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_abline(slope = 1, intercept = 0, col = "red", linetype = "dashed", size = 0.7) +
    ggtitle(paste0("QQ-plot under ", title, " distribution", taxa)) + theme_bw()
}

qqplot.zinb <- function(values, rng = seq(min(values), max(values), by = 1), 
                        estimates, taxa = NULL) {
  values = sort(values)
  ecdf   = sapply(rng, function(s) mean(values <= s))
  cdf    = pZINBI(rng, mu = estimates["mu"], sigma = estimates["sig"], nu = estimates["pi"])
  taxa   = if (!is.null(taxa)) paste0(" - ", taxa) else NULL
  
  data.frame(value = rng, ecdf = ecdf, cdf =  cdf) %>%
    ggplot(aes(value, ecdf)) +
    geom_col() +
    geom_line(aes(value, cdf), col = "red", size = 1, alpha = 0.5) +
    # geom_point(aes(value, cdf), col = "red", size = 1) +
    xlab("counts") + 
    ylab("cumulative distribution") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggtitle(paste0("Empirical (bars) v.s. model (red line) CDF", taxa)) + theme_bw()
}


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
full.beta  <- gamlss(tpm.i[tpm.i > 0]/const ~ 1, family = BE(sigma.link = "log"), control = con)
full.gamma <- gamlss(tpm.i[tpm.i > 0] ~ 1, family = GA(sigma.link = "log"), control = con)
full.ln    <- gamlss(tpm.i[tpm.i > 0] ~ 1, family = LOGNO(sigma.link = "log"), control = con)
full.zinb    <- gamlss(as.integer(tpm.i) ~ 1, family = ZINBI(sigma.link = "log"), control = con)

a1 = param.beta(full.beta);   nm.beta  = names(a1)
a2 = param.gamma(full.gamma); nm.gamma = names(a2)
a3 = param.ln(full.ln);       nm.ln    = names(a3)
# a4 = param.zinb(full.zinb)

gof.beta =
  matrix(NA, nrow = n.test, ncol = length(a1) + 3, 
         dimnames = list(i.taxa, c(nm.beta, "zero.prop", "ks.coef", "ks.pval")))
gof.gamma =
  matrix(NA, nrow = n.test, ncol = length(a2) + 3, 
         dimnames = list(i.taxa, c(nm.gamma, "zero.prop", "ks.coef", "ks.pval")))
gof.ln =
  matrix(NA, nrow = n.test, ncol = length(a3) + 3, 
         dimnames = list(i.taxa, c(nm.ln, "zero.prop", "ks.coef", "ks.pval")))


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
  
  # fitting
  full.beta  <- gamlss(tpm.i[tpm.i > 0]/const ~ 1, family = BE(sigma.link = "log"), control = con)
  full.gamma <- gamlss(otu.i[otu.i > 0] ~ 1, family = GA(sigma.link = "log"), control = con)
  full.ln    <- gamlss(otu.i[otu.i > 0] ~ 1, family = LOGNO(sigma.link = "log"), control = con)
  if (i <= 5) full.zinb  <- try(gamlss(as.integer(otu.i) ~ 1, family = ZINBI(sigma.link = "log"), control = con))
  
  gof.beta [i, nm.beta]    = param.beta(full.beta)
  gof.gamma[i, nm.gamma]   = param.gamma(full.gamma)
  gof.ln[i, nm.ln]         = param.ln(full.ln)
  gof.beta [i, "zero.prop"]= gof.gamma[i, "zero.prop"] = gof.ln[i, "zero.prop"] = mean(otu.i == 0)
  
  # KS test for nonzeros / Beta (KS does not allow ties.)
  tmp.ks <- ks.test(tpm.i[tpm.i > 0]/const, "pbeta", shape1 = gof.beta[i, "alpha"], shape2 = gof.beta[i, "beta"])
  gof.beta[i, "ks.coef"] = tmp.ks$statistic %>% as.numeric
  gof.beta[i, "ks.pval"] = tmp.ks$p.value %>% as.numeric
  
  # KS test for nonzeros / Gamma (KS does not allow ties.)
  tmp.ks <- ks.test(otu.i[otu.i > 0], "pGA", mu = gof.gamma[i, "mu"], sigma = gof.gamma[i, "sig"])
  gof.gamma[i, "ks.coef"] = tmp.ks$statistic %>% as.numeric
  gof.gamma[i, "ks.pval"] = tmp.ks$p.value %>% as.numeric
  
  # KS test for nonzeros / log-normal (KS does not allow ties.)
  tmp.ks <- ks.test(otu.i[otu.i > 0], "pLOGNO", 
                    mu = gof.ln[i, "mu"], sigma = gof.ln[i, "sig"])
  gof.ln[i, "ks.coef"] = tmp.ks$statistic %>% as.numeric
  gof.ln[i, "ks.pval"] = tmp.ks$p.value %>% as.numeric
  
  
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
    qqplot.zinb(values = as.integer(tpm.i), estimates = param.zinb(full.zinb))
  
  pp <- gridExtra::grid.arrange(p1, p2, p3, p4, 
                                top = grid::textGrob(paste0("Fitted models for ", taxa.i)))
  ggsave(paste0("figure/C0102QQplot_zoe", zoe, "_", nrm, "_",i, ".png"), plot = pp)
  }
  if (! i %% 10) {cat("i = ", i, 
                      ", mean.beta.ks.p = ", sprintf("%.5f", mean(gof.beta[, "ks.pval"], na.rm = T)), 
                      " mean.gamma.ks.p = ", sprintf("%.5f", mean(gof.gamma[, "ks.pval"], na.rm = T)), 
                      " mean.ln.ks.p = ", sprintf("%.5f", mean(gof.ln[, "ks.pval"], na.rm = T)), "\n")}
}


saveRDS(list(beta = gof.beta, gamma = gof.gamma, ln = gof.ln),
        paste0("output/gof_zoe", zoe, "_", nrm, ".rds"))

if (0) {
  for (zoe in 1:2) {
    gof.total <- NULL
    for (nrm in c("rpk", "tpm5", "asn")) {
      nrm2 = switch(nrm, "rpk" = "RPK", "tpm5" = "TPM", "asn" = "arcsin")
      gof <- readRDS(paste0("output/gof_zoe", zoe, "_", nrm, ".rds"))
      gof.total <- 
        rbind(gof.total,
              gof.beta[, "ks.pval"] %>% data.frame(pval = ., distribution = "Beta", nrm = nrm2),
              gof.ln[, "ks.pval"] %>% data.frame(pval = ., distribution = "Log-normal", nrm = nrm2),
              gof.gamma[, "ks.pval"] %>% data.frame(pval = ., distribution = "Gamma", nrm = nrm2))
    }
    gof.total %>% 
      ggplot(aes(pval)) + 
      facet_grid(distribution ~ nrm) + 
      geom_histogram(binwidth = 0.01) + 
      ggtitle("Kolmogorov-Smirnov test p-value histogram for Beta, Log-normal and Gamma distribution") +
      geom_vline(xintercept = 0.05, col = "red") + 
      annotate("text", x = 0.75, y = 25, label = paste0("%(p < 0.05) = ", mean(gof.gamma[, "ks.pval"] < 0.05)))
    ggsave(paste0("figure/C0102KS_zoe",zoe, "-p-Histogram_.png"))
  }
}

# h1 <-
#   gof.beta[, "ks.pval"] %>% data.frame(x = .) %>% ggplot(aes(x)) + 
#     geom_histogram(binwidth = 0.01) + ggtitle("KS p-value histogram for Beta distribution") +
#     geom_vline(xintercept = 0.05, col = "red") + xlab("p-values of the Kolmogorov-Smirnov tests")
#     annotate("text", x = 0.75, y = 25, label = paste0("%(p < 0.05) = ", mean(gof.beta[, "ks.pval"] < 0.05)))
# h2 <-
#   gof.ln[, "ks.pval"] %>% data.frame(x = .) %>% ggplot(aes(x)) + 
#     geom_histogram(binwidth = 0.01) + ggtitle("KS p-value histogram for Log-normal distribution") +
#     geom_vline(xintercept = 0.05, col = "red") + 
#     annotate("text", x = 0.75, y = 25, label = paste0("%(p < 0.05) = ", mean(gof.ln[, "ks.pval"] < 0.05)))
# h3 <-
#   gof.gamma[, "ks.pval"] %>% data.frame(x = .) %>% ggplot(aes(x)) + 
#     geom_histogram(binwidth = 0.01) + ggtitle("KS p-value histogram for Gamma distribution") +
#     geom_vline(xintercept = 0.05, col = "red") + 
#     annotate("text", x = 0.75, y = 25, label = paste0("%(p < 0.05) = ", mean(gof.gamma[, "ks.pval"] < 0.05)))
# hh <- gridExtra::grid.arrange(h1, h2, h3, top = grid::textGrob("KS test p-value histograms "))
# ggsave(paste0("figure/C0102KS_zoe",zoe, "_", nrm, "-p-Histogram_.png"), plot = hh)
#
#
# g1 <- 
#   gof.gamma %>% 
#   as.data.frame %>% 
#   ggplot(aes(zero.prop, ks.pval)) + 
#   geom_point() +
#   geom_smooth(method = "loess")
# 
# g2 <- 
#   gof.beta %>% 
#   as.data.frame %>% 
#   ggplot(aes(zero.prop, ks.pval)) + 
#   geom_point() +
#   geom_smooth(method = "loess")
# 
# g3 <- 
#   gof.ln %>% 
#   as.data.frame %>% 
#   ggplot(aes(zero.prop, ks.pval)) + 
#   geom_point() +
#   geom_smooth(method = "loess")
# 
# gg <- gridExtra::grid.arrange(g1, g2, g3, top = grid::textGrob("KS test p-value vs zero proportion
#                                                                    \nGamma / Beta / LogNormal"))
# ggsave(paste0("figure/C0102KS_zoe",zoe, "_", nrm, "-p_vs_zero-Histogram_.png"), plot = gg)
