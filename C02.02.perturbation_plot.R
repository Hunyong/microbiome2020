library(tidyverse)
library(latex2exp)
source("C01.02.simulation.setup.R")

pertplot <- function(model= "ZILN", size)
{
  j.index <- c(1,5,3)
  k.index = c(7,9,10,12,25,27,28,30,43,45,46,48)
  ylim =c(0,1)
  parameter = switch(model, 
                     zinb = parameterNB2, 
                     zig = parameterLN2, 
                     ziln = parameterLN2)
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
          result <- readRDS(paste0("stat-n",size,"-pert",pe,"-ziln-",i,".",j,".",k,".rds"))
          result.stat <- data.frame(result$stat)
          
          tmp <- result.stat%>% mutate ("i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),"effect" = as.character(result$setting$delta[4]),"pert"=pe)%>%dplyr::select("LB.glob","i","j","k","batch_f","effect","pert")
          
          res <- rbind(res,tmp[1,])
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
  res %>%
    ggplot(aes(factor(k),LB.glob,fill = batch_f)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(ylim) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k) +
    scale_fill_manual(name = TeX("Batch effects ($\\kappa_\\mu, \\kappa_\\theta, \\kappa_{\\pi}$)"),
                      values=c("K1 (0, 0, 0)"  = "dodgerblue",
                               "K3 (1, 1, -1)" = "chartreuse3",
                               "K5 (1, -1, -1)" = "tomato1")) +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    facet_grid(cols = vars(perturbation_f), rows = vars(effect_f), labeller = label_parsed) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> plb
  
  ggsave(file = "figure/LB_pert_size",size,".png", plb, width = 20, height=12)
  plb
}

pertplot(80)
pertplot(400)
