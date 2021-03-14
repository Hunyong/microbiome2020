library(tidyverse)
library(latex2exp)
library(ggh4x) # for nested facetting
source("C01.02.simulation.setup.R")

fullplot <- function(size,model)
{
  
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
  k.index = dim(parameter)[1]
  param.k = apply(parameter[,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  ylim = c(0,1)
  
  
  res <- NULL
  for(i in 1:10)
  {
    
    for(j in 1:5)
    {
      cat("\ni = ",i,"j = ",j, "k = ")
      
      for(k in 1:k.index)
      {
        cat(k , " ")
        result <- readRDS(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
        result.stat <- data.frame(result$stat)
        
        tmp <- 
          result.stat %>% 
          mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                  "i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),
                  "effect" = as.character(result$setting$delta[4])) %>%
          dplyr::select(LB, LN, MAST, KW, `KW-II`, DS2, `DS2ZI`, MGS, i, j, k, batch_f, effect)
        
        res <- rbind(res,tmp[1,])
      }
    }
  }
  res <- res %>% gather(key = "method", value = "p.value",`LB`,`LN`,`MAST`,`KW`,`KW-II`,`DS2`, `DS2ZI`, `MGS`)
  
  res$method_f = factor(res$method,
                             levels = c("LN", "LB", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS"),
                             labels = c("LN", "LB", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS"))
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
  res[res$method == "MGS" & res$j != 1, "p.value"] <- NA #NA for MGS with batch effects
# res.tmp2 <<- res
  res %>%
    ggplot(aes(factor(k), p.value,fill = batch_f)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(ylim) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k,guide =guide_axis(n.dodge = 2)) +
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
pn1 <- fullplot(80,model="ziln")
ggsave(file="figure/ziln_full_size80.png",pn1,width = 20,height = 16)
pn2 <- fullplot(400,model="ziln")
ggsave(file="figure/ziln_full_size400.png",pn2,width = 20,height = 16)


p1 <- fullplot(80,model="zinb")
ggsave(file="figure/zinb_full_size80.png",p1,width = 20,height = 16)
p2 <- fullplot(400,model="zinb")
ggsave(file="figure/zinb_full_size400.png",p2,width = 20,height = 16)

pp1 <- fullplot(80,model="zig")
ggsave(file="figure/zig_full_size80.png",pp1,width = 20,height = 16)
pp2 <- fullplot(400,model="zig")
ggsave(file="figure/zig_full_size400.png",pp2,width = 20,height = 16)


####power plot####

powerplot <- function(model,size, width = 20, height=12, delta.base = TRUE,
                      fn = paste0("figure/", model,"_power_size",size, if (!delta.base) "_effectSize(no_batch)", ".png"))
{
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
  k.index = k.core # c(7,9,10,12,25,27,28,30,43,45,46,48)
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  
  if (delta.base) {
    i.rng = c(2,3,4,6,8)
    disease.levels = c("Effect_mu(D>H)", "Effect_theta(D>H)", "Effect_pi(D<H)",
                       "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)")
    disease.labels = c(TeX("D2 ($\\mu_D$>$\\mu_{H}$)"), 
                       TeX("D3 ($\\theta_D$>$\\theta_{H}$)"),
                       TeX("D4 ($\\pi_D$<$\\pi_{H}$)"),
                       TeX("D6 ($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                       TeX("D8 ($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)"))
    j.rng = c(1,5,3)
    batch.levels = c("no batch effect", "large(+,+,-) batch effect", "large(+,-,-) batch effect")
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
  
  res3 <- NULL
  for(i in i.rng)
  {
    
    for(j in j.rng)
    {
      cat("\ni = ",i,"j = ",j, "k = ")
      for(k in k.index)
      {
        cat(k, " ")
        result <- readRDS(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
        result.stat <- data.frame(result$stat)
        
        tmp <- 
          result.stat%>% 
          mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                  "i" = i,"j" = j,"k" = k,
                  "batch" = as.character(result$setting$kappa[4]),
                  "effect" = as.character(result$setting$delta[4]))%>%
          dplyr::select("LB", "LN", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS", "i","j","k","batch","effect")
        
        res3 <- rbind(res3,tmp[1,])
      }
    }
  }
  res3[res3$method_f == "MGS" & res3$j != 1, "p.value"] <- NA #NA for MGS with batch effects
  
  res3 <- res3 %>% 
    gather(key = "method", value = "p.value",
           `LB`,`LN`,`MAST`,`KW`,`KW-II`,`DS2`, `DS2ZI`, `MGS`)
  res3$method_f = factor(res3$method,
                         levels = c("LN", "LB", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS"),
                         labels = c("LN", "LB", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS"))
  res3$batch_f = factor(res3$batch, levels = batch.levels, labels = batch.labels)
  res3$effect_f = factor(res3$effect, levels = disease.levels, labels = disease.labels)
  if (!delta.base) res3$effect2_f = factor(res3$effect, levels = disease.levels, labels = disease2.labels)
  res3[res3$method == "MGS" & res3$j != 1, "p.value"] <- NA #NA for MGS with batch effects
# tmp.res <<- res3  
  res3 %>%
    ggplot(aes(factor(k), p.value, fill = batch_f)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(c(0,1)) + 
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
    {if (!delta.base) guides(fill = FALSE)} +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    {if (!delta.base) facet_nested(method_f ~ effect2_f + effect_f, labeller = label_parsed)} +
    {if (delta.base) facet_grid(rows = vars(method_f), cols = vars(effect_f), labeller = label_parsed)} +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  ggsave(file = fn, p, width = width, height= height)
  p
}



powerplot(model = "ziln",size =400)
powerplot(model = "ziln",size =80)
powerplot(model = "ziln",size =400, delta.base = FALSE)
powerplot(model = "ziln",size =80, delta.base = FALSE)

powerplot(model = "ziln",size =400, width = 10, height = 7, 
          fn = paste0("figure/ziln_power_size400_2.png"))
powerplot(model = "ziln",size =80, width = 10, height = 7, 
          fn = paste0("figure/ziln_power_size80_2.png"))

powerplot(model = "zig",size =400)
powerplot(model = "zig",size =80)
powerplot(model = "zig",size =400, width = 10, height = 7, 
          fn = paste0("figure/zig_power_size400_2.png"))
powerplot(model = "zig",size =80, width = 10, height = 7, 
          fn = paste0("figure/zig_power_size80_2.png"))
powerplot(model = "zinb",size =400)
powerplot(model = "zinb",size =80)
powerplot(model = "zinb",size =400, width = 10, height = 7, 
          fn = paste0("figure/zinb_power_size400_2.png"))
powerplot(model = "zinb",size =80, width = 10, height = 7, 
          fn = paste0("figure/zinb_power_size80_2.png"))


powercurve <- function(model, width = 20, height=12,
                       fn = paste0("figure/", model,"_power_size",size, "_curve", ".png"))
{
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
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
  
  res3 <- NULL
  for (size in n.rng) {
    for (i in i.rng) {
      for (j in j.rng) {
        cat("\nn = ", size, "i = ", i, "j = ", j, "k = ")
        for(k in k.rng) {
          cat(k, " ")
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
            dplyr::select(cutoff, "LB", "LN", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS", 
                          n, i, j, k, batch, effect)
          res3 <- rbind(res3, tmp)
        }
      }
    }
  }
  # res3[res3$method_f == "MGS" & res3$j != 1, "p.value"] <- NA #NA for MGS with batch effects # not needed for this fn.
  
  res3 <- res3 %>% 
    gather(key = "method", value = "rejection.rate",
           `LB`,`LN`,`MAST`,`KW`,`KW-II`,`DS2`, `DS2ZI`, `MGS`)
  res3$method_f = factor(res3$method,
                         levels = c("LN", "LB", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS"),
                         labels = c("LN", "LB", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS"))
  res3$batch_f = factor(res3$batch, levels = batch.levels, labels = batch.labels)
  res3$effect_f = factor(res3$effect, levels = disease.levels, labels = disease.labels)
  # res3$effect2_f = factor(res3$effect, levels = disease.levels, labels = disease2.labels)
  res3$baseline_f = factor(res3$k, levels = k.rng, labels = base.labels)
  res3$size_f = factor(res3$n, levels = n.rng, labels = n.labels)
  
  # res3[res3$method == "MGS" & res3$j != 1, "p.value"] <- NA #NA for MGS with batch effects
  # tmp.res <<- res3  
  # res3 %>%
  res3 %>% 
    ggplot(aes(cutoff, rejection.rate, col = method_f)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, col = "gray") +
    geom_vline(xintercept=0.05, col="black", linetype = 2) + 
    ylim(c(0,1)) + xlim(c(0,1)) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    guides(fill = FALSE) +
    xlab("cut-off values") +
    ylab("power (rejection rate)") +
    facet_nested(baseline_f ~ size_f + effect_f, labeller = label_parsed) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  ggsave(file = fn, p, width = width, height= height)
  p
}
powercurve(model = "ziln")
