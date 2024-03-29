library(tidyverse)
library(latex2exp)
source("C01.02.simulation.setup.R")

pertplot <- function(model= "ziln", size, res.tmp = TRUE, stop.if.absent = TRUE)
{
  j.index <- c(1,5,3)
  k.index = c(7,9,10,12,25,27,28,30,43,45,46,48)
  ylim =c(0,1)
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  res <- NULL
  for(pe in c(0,0.25,0.5))
  {
    for(i in c(1,2,3,4,6,8))
    {
      
      for(j in j.index)
      {
        for(k in k.index)
        {
          fn.tmp <- paste0("output/stat-n",size,"-pert",pe,"-ziln-",i,".",j,".",k,".rds")
          if (file.exists(fn.tmp)) {
            result <- readRDS(fn.tmp)
            result.stat <- data.frame(result$stat)
            
            tmp <- result.stat%>% mutate ("i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),"effect" = as.character(result$setting$delta[4]),"pert"=pe)%>%dplyr::select("LB.glob","i","j","k","batch_f","effect","pert")
            
            res <- rbind(res,tmp[1,])
          } else {
            if (stop.if.absent) stop("Not available")
            cat("(Not available) ")
          }
        }
      }
    }
  }
  
  res$perturbation_f = factor(res$pert,
                              levels = c( 0,0.25, 0.5),
                              labels = c(TeX("perturbation=0"),
                                         TeX("perturbation=0.25"),
                                         TeX("perturbation=0.5")))
  res$effect_f = factor(res$effect,
                        levels = c("Effect_null", "Effect_mu(D>H)", 
                                   "Effect_theta(D>H)", "Effect_pi(D<H)",
                                   "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"), 
                        labels = c(TeX("D1 (null)"),
                                   TeX("D2 ($\\mu_D$>$\\mu_{H}$)"), 
                                   TeX("D3 ($\\theta_D$>$\\theta_{H}$)"),
                                   TeX("D4 ($\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D6 ($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D8 ($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)")))
  res$batch_f = factor(res$batch_f,
                       levels = c("no batch effect", 
                                  "large(+,+,-) batch effect",
                                  "large(+,-,-) batch effect"),
                       labels = c("K1 (0, 0, 0)",
                                  "K3 (1, 1, -1)",
                                  "K5 (1, -1, -1)"))
  res$k <- factor(res$k)
  if (res.tmp) res.tmp <<- res
    
  res %>%
    ggplot(aes(k, LB.glob, fill = batch_f)) +
    geom_bar(stat="identity", position = position_dodge(width = .8)) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(ylim) + 
    geom_point(aes(col = batch_f), position = position_dodge(width = .8), shape = 15, size = 0.5) +
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k) +
    scale_fill_manual(name = TeX("Batch effects ($\\kappa_\\mu, \\kappa_\\theta, \\kappa_{\\pi}$)"),
                      values=c("K1 (0, 0, 0)"  = "dodgerblue",
                               "K3 (1, 1, -1)" = "chartreuse3",
                               "K5 (1, -1, -1)" = "tomato1")) +
    scale_color_manual(values=c("K1 (0, 0, 0)"  = "dodgerblue",
                               "K3 (1, 1, -1)" = "chartreuse3",
                               "K5 (1, -1, -1)" = "tomato1")) +
    guides(col = FALSE) +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    facet_grid(cols = vars(perturbation_f), rows = vars(effect_f), labeller = label_parsed) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> plb
  
  ggsave(file = paste0("figure/LB_pert_size",size,".pdf"), plb, width = 20, height = 12)
  plb
}

pertplot(model= "ziln", size = 80, stop.if.absent = FALSE)
pertplot(model= "ziln", size = 400, stop.if.absent = FALSE)
