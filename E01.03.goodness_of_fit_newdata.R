### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra); library(cowplot)
library(gamlss)
source("F00.00.generic.R")
source("F01.01.base.R")
source("F01.01.goodness_of_fit.R")
args = commandArgs(trailingOnly=TRUE)  # passed from script

nrm = "tpm"
type = "gene"
for (type in c("genebact", "bact", "gene")) {
  print(type)
  type.full = switch(type, gene = "geneProp.marginal", 
                     genebact = "geneProp.joint", bact = "bactProp.marginal")
  
  fn <- sprintf("Nature2019data/data.ecs_relab.%s.RNA.rds", type.full)
  data <- readRDS(fn)
  DataMeta = data$meta
  DataMeta <-
    DataMeta %>% 
    mutate(group_batch =  ifelse(grepl("MGH", site_name), 1, 2), # dichotomize sites into either related to MGH (MGH and MGH pediatrics)
           group_disease = ifelse(diagnosis == "nonIBD", "H", "D"),
           group = paste0(group_disease, group_batch))
  # DRNA = "RNA"; DR.no = 2
  RNA      = data$otu  # proportion (0~1)
  
  # scale1 <- 5e+5 * if (type == "bact") {1/2500} else 1
  # scale1 <- if (type == "bact") 10 else 1e+3
  scale1 <- if (type == "bact") 500 else 1e+5
  
  DataComp    = RNA  # for LB tests
  DataTPM     = RNA * scale1
  if (nrm == "asin") {
    DataTPM = asn(DataTPM/scale1) * scale1
    DataComp = asn(DataComp)
  }
  
  
  # screening and sampling
  prev.filter = DataTPM %>% apply(1, function(x) mean(x>0)) %>% {which(. >= 0.1)}
  n.test   = min(300, length(prev.filter))
  set.seed(1)
  i.sample    = sample(prev.filter, n.test)
  i.taxa      = rownames(DataTPM)[i.sample]
  
  # Empty shell
  tpm.i = DataTPM[i.sample[1], ]
  comp.i = DataComp[i.sample[1], ]
  con = gamlss.control(n.cyc = 100, trace = FALSE)
  a1 <- est.beta(comp.i[comp.i>0]);   nm.beta  = names(a1)
  a2 <- est.gamma(tpm.i[tpm.i>0]);        nm.gamma = names(a2)
  a3 <- est.ln(tpm.i[tpm.i>0]);           nm.ln    = names(a3)
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
    j = i.sample[i]
    tpm.i[] <- as.numeric(DataTPM[j, ])
    comp.i[] <- as.numeric(DataComp[j, ])
    if (nrm == "asin") tpm.i = asn(tpm.i / scale1) * scale1
    taxa.i  <- rownames(DataTPM)[j] %>% gsub(" \\(TOTAL\\)", "", .) %>% gsub("\\_", " ", .)
    # name.i  <- names(DataRPK)[j]
    # if (mean(dat.reg$composition > 0) < 0.1) next
    
    # fitting & KS tests
    args.gof = list(return.est = TRUE, lilliefors = TRUE, n.lilliefors = 300, ks.pval = TRUE)
    gof.beta [i, c(nm.beta, nm.ks)]  = do.call(ks.empirical, c(list(comp.i[comp.i > 0], model = "beta"), args.gof))
    gof.gamma[i, c(nm.gamma, nm.ks)] = do.call(ks.empirical, c(list(tpm.i[tpm.i > 0], model = "gamma"), args.gof))
    gof.ln[i, c(nm.ln, nm.ks)]       = do.call(ks.empirical, c(list(tpm.i[tpm.i > 0], model = "ln"), args.gof))
    gof.beta [i, "zero.prop"]= gof.gamma[i, "zero.prop"] = gof.ln[i, "zero.prop"] = mean(tpm.i == 0)
    
    
    ## plots
    if (i <= 0) {
      p1  <- 
        qqplot1(values = tpm.i[tpm.i > 0]/scale1, ks.pval = gof.beta[i, "ks.pval"],
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
      ggsave(paste0("figure/C0102QQplot_newdata_", type, "_", nrm, "_",i, ".png"), plot = pp,
             width = 12, height = 9)
    }
    if (! i %% 10) {cat("i = ", i, 
                        ", mean.beta.ks.p = ", sprintf("%.5f", mean(gof.beta[, "ks.pval"], na.rm = T)), 
                        " mean.gamma.ks.p = ", sprintf("%.5f", mean(gof.gamma[, "ks.pval"], na.rm = T)), 
                        " mean.ln.ks.p = ", sprintf("%.5f", mean(gof.ln[, "ks.pval"], na.rm = T)), "\n")}
    saveRDS(list(beta = gof.beta, gamma = gof.gamma, ln = gof.ln),
            paste0("output/gof_newdata_", type, "_", nrm, ".tmp.rds"))
  }
  
  saveRDS(list(beta = gof.beta, gamma = gof.gamma, ln = gof.ln),
          paste0("output/gof_newdata_", type, "_", nrm, ".rds"))

}

if (0) {
  gof.total <- NULL
  for (nrm in c("tpm", "asin")) {
    nrm2 = switch(nrm, "tpm" = "TPM", "asin" = "arcsin")
    gof <- readRDS(paste0("output/gof_newdata_", nrm, ".rds"))
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
           pval = 0.5, count = 100,
           reject = ifelse(reject == "%(p < 0.05) = 0", "%(p < 0.05) < 0.001", reject))
  gof.total %>% 
    mutate(nrm = factor(nrm, levels = c("TPM", "arcsin"))) %>% 
    ggplot(aes(pval)) + 
    facet_grid(distribution ~ nrm) + 
    geom_histogram(binwidth = 0.01) + 
    # ggtitle("Kolmogorov-Smirnov test p-value histogram for Beta, Log-normal and Gamma distribution") +
    xlab("KS (Lilliefors) test p-values") + ggtitle("new data (n = 104)") + 
    geom_vline(xintercept = 0.05, col = "red") + 
    geom_text(data = gof.stat, aes(pval, count, label = reject)) +
    theme_bw()
  # annotate("text", x = 0.75, y = 25, label = paste0("%(p < 0.05) = ", mean(gof.gamma[, "ks.pval"] < 0.05)))
  ggsave(paste0("figure/C0102KS_newdata-p-Histogram_.png"))

}

