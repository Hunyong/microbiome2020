library(tidyverse)
library(latex2exp)
library(ggh4x) # for nested facetting
source("C01.02.simulation.setup.R")


####power plot####

powerplot <- function(model, size,  width = 12,  height = 8, metrics.name = "type1error", delta.base = TRUE, 
                      fn = paste0("figure/",  model, "_power_size", size,  if (!delta.base) "_effectSize(no_batch)",  
                                  if (include.null) "_with_null",  ".pdf"),
                      res.tmp = TRUE, include.null = FALSE, stop.if.absent = TRUE)
{
  require(tidyr)
  parameter = switch(model,  
                     zinb = parameterNB,  
                     zig = parameterLN,  
                     ziln = parameterLN)
  k.index = c(1)
  param.k = apply(parameter[k.index, -1],  1,  function(x) paste0("(",  paste(x,  collapse=",  "),  ")"))
  
  if (delta.base) {
    i.rng = c(2, 3, 4, 6, 8)
    disease.levels = c("Effect_mu(D>H)",  "Effect_theta(D>H)",  "Effect_pi(D<H)", 
                       "Effect_mu(D>H).pi(D<H)",  "Effect_mu(D>H),pi(D>H)")
    disease.labels = c(TeX("D2 ($\\mu_D$>$\\mu_{H}$)"),  
                       TeX("D3 ($\\theta_D$>$\\theta_{H}$)"), 
                       TeX("D4 ($\\pi_D$<$\\pi_{H}$)"), 
                       TeX("D6 ($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"), 
                       TeX("D8 ($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)"))
    j.rng = c(1)
    batch.levels = c("no batch effect",  "large(+,+,-) batch effect", "large(+,-,-) batch effect")
    batch.labels = c("K1 (0, 0, 0)", "K3 (1, 1, -1)", "K5 (1, -1, -1)")
    
  } else {
    i.rng = c(11, 13, 15, 2, 4, 6, 12, 14, 16)
    disease.levels = c("Effect_mu_small", "Effect_mu(D>H)", "Effect_mu_large", 
                       "Effect_pi_small", "Effect_pi(D<H)", "Effect_pi_large", 
                       "Effect_mu.pi_small", "Effect_mu(D>H).pi(D<H)", "Effect_mu.pi_large")
    disease.labels = c(TeX("D11. $\\delta = (0.5, 0, 0)$"),
                       TeX("D2. $\\delta = (1, 0, 0)$"),
                       TeX("D13. $\\delta = (2, 0, 0)$"),
                       TeX("D13. $\\delta = (0, 0, -0.5)$"),
                       TeX("D4. $\\delta = (0, 0, -1)$"),
                       TeX("D14. $\\delta = (0, 0, -2)$"),
                       TeX("D15. $\\delta = (0.5, 0, -0.5)$"),
                       TeX("D6. $\\delta = (1, 0, -1)$"),
                       TeX("D16. $\\delta = (2, 0, -2)$"))
    disease2.labels = c(TeX("$\\mu$ effect (D2 and its variants)"),
                       TeX("$\\mu$ effect (D2 and its variants)"),
                       TeX("$\\mu$ effect (D2 and its variants)"),
                       TeX("$\\pi$ effect (D4 and its variants)"),
                       TeX("$\\pi$ effect (D4 and its variants)"),
                       TeX("$\\pi$ effect (D4 and its variants)"),
                       TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
                       TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
                       TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"))
    j.rng = 1
    batch.levels = c("no batch effect")
    batch.labels = c("K1 (0, 0, 0)")
  }
  if (include.null) {
    i.rng = c(1, i.rng)
    disease.levels = c("Effect_null", disease.levels)
    disease.labels = c("D1 (null)", disease.labels)
  }
  dict = list("sensitivity"=1, "type1error" =2, "FDR"=3, "accuracy"=4, "AUC"=5)

  res <- NULL
  for(i in i.rng)
  {
    
    for(j in j.rng)
    {
      cat("\ni = ",i,"j = ",j, "k = ")
      for(k in k.index)
      {
        cat(k, " ")
        fn.tmp <- paste0("output/stat-n",size,"-pert0.5-signal0.1-",model,"-",i,".",j,".",k,".rds")
        if (file.exists(fn.tmp)) {
          result <- readRDS(fn.tmp)
          result.metrics <- data.frame(result$metrics %>% t)
          
          tmp <- 
            result.metrics %>% 
            mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                    "i" = i,"j" = j,"k" = k,
                    "batch" = as.character(result$setting$kappa[4]),
                    "effect" = as.character(result$setting$delta[4]))%>%
            dplyr::select("LB", "LN", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS", 
                          "ANCOM", "i","j","k","batch","effect")
          
          res <- rbind(res,tmp[dict[[metrics.name]],])
        } else {
          if (stop.if.absent) stop("Not available")
          cat("(Not available) ")
        }
      }
    }
  }
  
  res <- res %>% 
    gather(key = "method", value = "metrics",
           `LB`,`LN`,`MAST`,`KW`,`KW-II`,`DS2`, `DS2ZI`, `MGS`, `ANCOM`)
  res$method_f = factor(res$method,
                         levels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II"),
                         labels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II"))
  res$batch_f = factor(res$batch, levels = batch.levels, labels = batch.labels)
  res$effect_f = factor(res$effect, levels = disease.levels, labels = disease.labels)
  if (!delta.base) res$effect2_f = factor(res$effect, levels = disease.levels, labels = disease2.labels)
  # res[res$method %in% c("MGS", "ANCOM", "ANCOM.sz") 
  #     & res$j != 1, "p.value"] <- NA #NA for MGS, ANCOM, and ANCOM.sz with batch effects
  res$k <- factor(res$k)
  if (res.tmp) res.tmp <<- res  
  
  res %>%
    ggplot(aes(k, metrics, fill = batch_f)) +
    geom_bar(stat="identity", position = position_dodge(width = .8)) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + #ylim(c(0, 1)) + 
    geom_point(aes(col = batch_f), position = position_dodge(width = .8), shape = 15, size = 0.1) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k) +
    {if (delta.base) 
      scale_fill_manual(name = TeX("Batch effects ($\\kappa_\\mu, \\kappa_\\theta, \\kappa_{\\pi}$)"),
                        values=c("K1 (0, 0, 0)"  = "dodgerblue",
                                 "K3 (1, 1, -1)" = "chartreuse3",
                                 "K5 (1, -1, -1)" = "tomato1")) 
      else
      scale_fill_manual(values=c("K1 (0, 0, 0)"  = "dodgerblue"))} +
    {if (delta.base) 
      scale_color_manual(values=c("K1 (0, 0, 0)"  = "dodgerblue",
                                 "K3 (1, 1, -1)" = "chartreuse3",
                                 "K5 (1, -1, -1)" = "tomato1")) 
      else
        scale_color_manual(values=c("K1 (0, 0, 0)"  = "dodgerblue"))} +
    {if (!delta.base) guides(fill = FALSE, color = FALSE)} +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    guides(col = FALSE) +
    {if (!delta.base) facet_nested(method_f ~ effect2_f + effect_f, labeller = label_parsed)} +
    {if (delta.base) facet_grid(rows = vars(method_f), cols = vars(effect_f),scales = "free_y", labeller = label_parsed)} +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  ggsave(file = fn, p, width = width, height= height)
  p
}


# modify to add the head points!!!!

# powerplot(model = "ziln", size = 400, stop.if.absent = FALSE)
powerplot(model = "ziln", size = 80, stop.if.absent = FALSE)
powerplot(model = "ziln", size = 80, include.null = TRUE, stop.if.absent = FALSE)
powerplot(model = "ziln", size = 400, delta.base = FALSE, stop.if.absent = FALSE) # disease effect sensitivity analysis
powerplot(model = "ziln", size = 80, delta.base = FALSE, stop.if.absent = FALSE) # disease effect sensitivity analysis

powerplot(model = "zig", size = 400, stop.if.absent = FALSE)
powerplot(model = "zig", size = 80, stop.if.absent = FALSE)

powerplot(model = "zinb", size = 400, stop.if.absent = FALSE)
powerplot(model = "zinb", size = 80, stop.if.absent = FALSE)


powercurve <- function(model,  width = 12,  height = 9,
                       fn = paste0("figure/", model,"_power_curve", ".pdf"), res.tmp = TRUE,
                       stop.if.absent = TRUE)
{
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
  k.index = k.core # c(7, 9, 10, 12, 25, 27, 28, 30, 43, 45, 46, 48)
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  
  {
    # i.rng = c(11, 13, 15, 2, 4, 6, 12, 14, 16)
    # disease.levels = c("Effect_mu_small", "Effect_mu(D>H)", "Effect_mu_large", 
    #                    "Effect_pi_small", "Effect_pi(D<H)", "Effect_pi_large", 
    #                    "Effect_mu.pi_small", "Effect_mu(D>H).pi(D<H)", "Effect_mu.pi_large")
    # disease.labels = c(TeX("D11. $\\delta = (0.5, 0, 0)$"),
    #                    TeX("D2. $\\delta = (1, 0, 0)$"),
    #                    TeX("D13. $\\delta = (2, 0, 0)$"),
    #                    TeX("D13. $\\delta = (0, 0, -0.5)$"),
    #                    TeX("D4. $\\delta = (0, 0, -1)$"),
    #                    TeX("D14. $\\delta = (0, 0, -2)$"),
    #                    TeX("D15. $\\delta = (0.5, 0, -0.5)$"),
    #                    TeX("D6. $\\delta = (1, 0, -1)$"),
    #                    TeX("D16. $\\delta = (2, 0, -2)$"))
    # disease2.labels = c(TeX("$\\mu$ effect (D2 and its variants)"),
    #                     TeX("$\\mu$ effect (D2 and its variants)"),
    #                     TeX("$\\mu$ effect (D2 and its variants)"),
    #                     TeX("$\\pi$ effect (D4 and its variants)"),
    #                     TeX("$\\pi$ effect (D4 and its variants)"),
    #                     TeX("$\\pi$ effect (D4 and its variants)"),
    #                     TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
    #                     TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
    #                     TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"))
    i.rng = c(1, 2, 4, 6)
    disease.levels = c("Effect_null", "Effect_mu(D>H)", "Effect_pi(D<H)", "Effect_mu(D>H).pi(D<H)")
    disease.labels = c(TeX("D1. $\\delta = (0, 0, 0)$"),
                       TeX("D2. $\\delta = (1, 0, 0)$"),
                       TeX("D4. $\\delta = (0, 0, -1)$"),
                       TeX("D6. $\\delta = (1, 0, -1)$"))
    k.rng = c(7, 25, 43) # c(7,9,10,12,25,27,28,30,43,45,46,48)
    base.labels = c(TeX("$(\\mu, \\theta, \\pi)$ = (1, 0.5, 0.6)"),
                    TeX("$(\\mu, \\theta, \\pi)$ = (1, 0.5, 0.75)"),
                    TeX("$(\\mu, \\theta, \\pi)$ = (1, 0.5, 0.9)"))
    j.rng = 1
    batch.levels = c("no batch effect")
    batch.labels = c("K1 (0, 0, 0)")
    
    n.rng = c(80, 400)
    n.labels = paste0("n==", n.rng)
  }
  
  res <- NULL
  for (size in n.rng) {
    for (i in i.rng) {
      for (j in j.rng) {
        cat("\nn = ", size, "i = ", i, "j = ", j, "k = ")
        for(k in k.rng) {
          cat(k, " ")
          fn.tmp = paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds")
          if (file.exists(fn.tmp)) {
            result <- readRDS(paste0("output/stat-n", size, "-pert0.5-", model, "-", i, ".", j, ".", k, ".rds"))
            # result.stat <- data.frame(result$stat)
            result.cdf <- data.frame(result$cdf %>% t)
            
            tmp <- 
              result.cdf%>% 
              mutate (cutoff = attr(result$cdf, "cutoff"),
                      "LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                      n = size, i = i, j = j, k = k,
                      "batch" = as.character(result$setting$kappa[4]),
                      "effect" = as.character(result$setting$delta[4])) %>%
              dplyr::select(cutoff, "LB", "LN", "MAST", "KW", "KW-II", "DS2", 
                            "DS2ZI", "MGS", "ANCOM",
                            n, i, j, k, batch, effect)
            res <- rbind(res, tmp)
          } else {
            if (stop.if.absent) stop("Not available")
            cat("(Not available) ")
          }
        }
      }
    }
  }
  # res[res$method_f == "MGS" & res$j != 1, "p.value"] <- NA #NA for MGS with batch effects # not needed for this fn.
  
  res <- res %>% 
    gather(key = "method", value = "rejection.rate",
                  `LB`,`LN`,`MAST`,`KW`,`KW-II`,`DS2`, `DS2ZI`, `MGS`, `ANCOM`)
  res$method_f = factor(res$method,
                         levels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II"),
                         labels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II"))
  res$batch_f = factor(res$batch, levels = batch.levels, labels = batch.labels)
  res$effect_f = factor(res$effect, levels = disease.levels, labels = disease.labels)
  # res$effect2_f = factor(res$effect, levels = disease.levels, labels = disease2.labels)
  res$baseline_f = factor(res$k, levels = k.rng, labels = base.labels)
  res$size_f = factor(res$n, levels = n.rng, labels = n.labels)
  
  # grid points for dots, different x-values for each method.
  res.points = 
    res %>% 
    dplyr::filter(((cutoff * 100) %% 9) == {as.numeric(method_f) %% 9} )
tmp.p <<- res.points
  # res[res$method %in% c("MGS", "ANCOM", "ANCOM.sz") & res$j != 1, "p.value"] <- NA #NA for MGS with batch effects
  if (res.tmp) res.tmp <<- res
  # res %>%
  res %>% 
    ggplot(aes(cutoff, rejection.rate, col = method_f, linetype = method_f)) +
    geom_line() +
    geom_point(data = res.points, size = 2,
               aes(cutoff, rejection.rate, col = method_f, shape = method_f)) +
    geom_abline(slope = 1, intercept = 0, col = "gray") +
    geom_vline(xintercept=0.05, col="black", linetype = 2) + 
    ylim(c(0,1)) + 
    # xlim(c(0,1)) + 
    scale_x_continuous(breaks = c(0, 0.03, 0.05, 0.10, 0.15, 0.2), limits = c(0, 0.2)) +
    scale_shape_manual(values = c(LN = 16, LB = 3, MAST = 1, KW = 18, `KW-II` = 3, DS2 = 8, `DS2ZI` = 3, MGS = 3, ANCOM = 16)) +
    scale_linetype_manual(values = c(LN = 1, LB = 3, MAST = 1, KW = 1, `KW-II` = 3, DS2 = 1, `DS2ZI` = 3, MGS = 1, ANCOM = 3)) +
    scale_color_manual(values = c(LN = "firebrick", LB = "firebrick", MAST = "darkseagreen4", KW = "darkorchid3", `KW-II` = "darkorchid3", 
                                  DS2 = "dodgerblue3", `DS2ZI` = "dodgerblue3", MGS = "goldenrod3", ANCOM = "goldenrod3" )) +
    guides(fill = FALSE, col = guide_legend(nrow = 1, title = NULL), 
           shape = guide_legend(nrow = 1, title = NULL), linetype = guide_legend(nrow = 1, title = NULL)) +
    xlab("cut-off values") +
    ylab("power (rejection rate)") +
    facet_nested(baseline_f ~ size_f + effect_f, labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle=90), plot.title = element_text(hjust = 0.5))  -> p
  ggsave(file = fn, p, width = width, height= height)
  p
}
powercurve(model = "ziln", stop.if.absent = FALSE)


#### Full plot####

fullplot <- function(size, model, res.tmp = TRUE, stop.if.absent = TRUE) { # res.tmp being globally assigned
  
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
  k.index = dim(parameter)[1]
  param.k = apply(parameter[,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  
  
  res <- NULL
  for(i in 1:10)
  {
    
    for(j in 1:5)
    {
      cat("\ni = ",i,"j = ",j, "k = ")
      
      for(k in 1:k.index)
      {
        cat(k , " ")
        fn.tmp <- paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds")
        if (file.exists(fn.tmp)) {
          result <- readRDS(fn.tmp)
          result.stat <- data.frame(result$stat)
          tmp <- 
            result.stat %>% 
            mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                    "i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),
                    "effect" = as.character(result$setting$delta[4])) %>%
            dplyr::select(LB, LN, MAST, KW, `KW-II`, DS2, `DS2ZI`, MGS, 
                          ANCOM.sz, ANCOM, i, j, k, batch_f, effect)
          
          res <- rbind(res,tmp[1,])
          
        } else {
          if (stop.if.absent) stop("Not available")
          cat("(Not available) ")
        }
      }
    }
  }
  res <- res %>% gather(key = "method", value = "p.value",`LB`,`LN`,`MAST`,`KW`,
                        `KW-II`,`DS2`, `DS2ZI`, `MGS`, `ANCOM.sz`, `ANCOM`)
  
  res$method_f = factor(res$method,
                        levels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM.sz", "ANCOM", "KW", "KW-II"),
                        labels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOMsz", "ANCOM", "KW", "KW-II"))
  res$effect_f = factor(res$effect,
                        levels = c("Effect_null", "Effect_mu(D>H)", 
                                   "Effect_theta(D>H)", "Effect_pi(D<H)",
                                   "Effect_mu(D>H).theta(D>H)", "Effect_mu(D>H).pi(D<H)", 
                                   "Effect_theta(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)",
                                   "Effect_mu(D>H).theta(D<H)", "Effect_theta(D<H).pi(D<H)"), 
                        labels = c(TeX("D1 (null)"),
                                   TeX("D2 ($\\mu_D$>$\\mu_{H}$)"), 
                                   TeX("D3 ($\\theta_D$>$\\theta_{H}$)"),
                                   TeX("D4 ($\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D5 ($\\mu_D$>$\\mu_H$, $\\theta_D$>$\\theta_{H}$)"),
                                   TeX("D6 ($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D7 ($\\theta_D$>$\\theta_H$, $\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D8 ($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)"),
                                   TeX("D9 ($\\mu_D$>$\\mu_H$, $\\theta_D$<$\\theta_{H}$)"),
                                   TeX("D10 ($\\theta_D$<$\\theta_H$, $\\pi_D$<$\\pi_{H}$)")))
  res$batch_f = factor(res$batch_f,
                       levels = c("no batch effect", 
                                  "small(+,+,-) batch effect",
                                  "large(+,+,-) batch effect",
                                  "large(+,-,-) batch effect",
                                  "small(+,-,-) batch effect"),
                       labels = c("K1 (0, 0, 0)",
                                  "K2 (0.5, 0.5, -0.5)",
                                  "K3 (1, 1, -1)",
                                  "K4 (0.5, -0.5, -0.5)",
                                  "K5 (1, -1, -1)"))
  # res.tmp <<- res
  # res[res$method %in% c("MGS", "ANCOM", "ANCOM.sz") 
  #     & res$j != 1, "p.value"] <- NA #NA for MGS, ANCOM, and ANCOM.sz with batch effects
  res$k <- factor(res$k)
  if (res.tmp) res.tmp <<- res
  res %>%
    ggplot(aes(k, p.value, fill = batch_f)) +
    geom_bar(stat="identity", position = position_dodge(width = .8)) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(c(0, 1)) + 
    geom_point(aes(col = batch_f), position = position_dodge(width = .8), shape = 15, size = 0.1) +
    # geom_text(data = res %>% filter(p.value==0), aes(label = p.value, col = batch_f), 
    #           position = position_dodge(width = .8), size = 0.5) +
    # geom_point(data = res %>% filter(p.value==0), aes(k, p.value, fill = batch_f, col = batch_f),
    #            shape = 21, size = 0.2,
    #            position = position_dodge()) +
    guides(col = FALSE) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k,guide = guide_axis(n.dodge = 2)) +
    scale_color_manual(values=c("K1 (0, 0, 0)"  = "dodgerblue",
                                "K2 (0.5, 0.5, -0.5)"  = "darkgreen",
                                "K3 (1, 1, -1)" = "chartreuse3",
                                "K4 (0.5, -0.5, -0.5)"  = "red3",
                                "K5 (1, -1, -1)" = "tomato1")) +
    scale_fill_manual(name = TeX("Batch effects ($\\kappa_\\mu, \\kappa_\\theta, \\kappa_{\\pi}$)"),
                      values=c("K1 (0, 0, 0)"  = "dodgerblue",
                               "K2 (0.5, 0.5, -0.5)"  = "darkgreen",
                               "K3 (1, 1, -1)" = "chartreuse3",
                               "K4 (0.5, -0.5, -0.5)"  = "red3",
                               "K5 (1, -1, -1)" = "tomato1")) +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    facet_grid(rows = vars(method_f), cols = vars(effect_f), labeller = label_parsed) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  
  p
}

#width = 32, height = 16
pn1 <- fullplot(80, model="ziln", stop.if.absent = FALSE)
ggsave(file="figure/ziln_full_size80.pdf", pn1, width = 30, height = 24)
pn2 <- fullplot(400, model="ziln", stop.if.absent = FALSE)
ggsave(file="figure/ziln_full_size400.pdf", pn2, width = 30, height = 24)


p1 <- fullplot(80, model="zinb", stop.if.absent = FALSE)
ggsave(file="figure/zinb_full_size80.pdf", p1, width = 30, height = 24)
p2 <- fullplot(400, model="zinb", stop.if.absent = FALSE)
ggsave(file="figure/zinb_full_size400.pdf", p2, width = 30, height = 24)

pp1 <- fullplot(80, model="zig", stop.if.absent = FALSE)
ggsave(file="figure/zig_full_size80.pdf", pp1, width = 30, height = 24)
pp2 <- fullplot(400, model="zig", stop.if.absent = FALSE)
ggsave(file="figure/zig_full_size400.pdf", pp2, width = 30, height = 24)
