library(tidyverse)
library(latex2exp)
source("C01.02.simulation.setup.R")

LBplot <- function(model,size)
{
  j.index <- c(1,5,3)
  k.index = c(7,9,10,12,25,27,28,30,43,45,46,48)
  parameter = switch(model, 
                     zinb = parameterNB, 
                     zig = parameterLN, 
                     ziln = parameterLN)
  
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  i=1
  res <- NULL
  for(j in j.index)
  {
    for(k in k.index)
    {
      result <- readRDS(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
      result.stat <- data.frame(result$stat)
      
      tmp <- result.stat%>% mutate ("i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),"effect" = as.character(result$setting$delta[4]))%>%dplyr::select("LB.nonz","LB.zero","LB.glob","LB.min","i","j","k","batch_f","effect")
      
      res <- rbind(res,tmp[1,])
    }
  }
  
  res <- res %>% gather(key = "method", value = "p.value",`LB.nonz`,`LB.zero`,`LB.glob`,`LB.min`)
  res$method_f = factor(res$method,
                        levels = c("LB.nonz", "LB.zero", "LB.glob", "LB.min"),
                        labels = c("LB - continuous","LB - discrete", "LB - global", "LB - minimum"))
  res$batch_f = factor(res$batch_f,
                       levels = c("no batch effect", 
                                  "large(+,+,-) batch effect",
                                  "large(+,-,-) batch effect"),
                       labels = c("K1 (0, 0, 0)",
                                  "K3 (1, 1, -1)",
                                  "K5 (1, -1, -1)"))
  res$effect_f = factor(res$effect,
                        levels = c("Effect_mu(D>H)", "Effect_theta(D>H)", "Effect_pi(D<H)",
                                   "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"), 
                        labels = c(TeX("D2 ($\\mu_D$>$\\mu_{H}$)"), 
                                   TeX("D3 ($\\theta_D$>$\\theta_{H}$)"),
                                   TeX("D4 ($\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D6 ($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                   TeX("D8 ($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)")))
  
  
  res %>%
    ggplot(aes(factor(k), p.value,fill = batch_f)) +
    #ggplot(aes(factor(k), p.value, col=batch, group=batch, fill = batch)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(c(0,0.2)) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k) +
    scale_fill_manual(name = TeX("Batch effects ($\\kappa_\\mu, \\kappa_\\theta, \\kappa_{\\pi}$)"),
                      values=c("K1 (0, 0, 0)"  = "dodgerblue",
                               "K3 (1, 1, -1)" = "chartreuse3",
                               "K5 (1, -1, -1)" = "tomato1")) +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    facet_grid(rows = vars(method_f)) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")-> p
  
  ggsave(file = paste0(model,"_LB_size",size,".pdf"), p, width = 10, height=12)
  
  p
}
LBplot(model = "ziln",size =400)

LBplot(model = "ziln",size =80)

LBplot(model = "zig",size =400)
LBplot(model = "zig",size =80)
LBplot(model = "zinb",size =400)
LBplot(model = "zinb",size =80)
  
