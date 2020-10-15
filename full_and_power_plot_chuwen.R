library(tidyverse)
library(latex2exp)
source("C01.02.simulation.setup.R")

fullplot <- function(size,model)
{
  
  parameter = switch(model, 
                     zinb = parameterNB2, 
                     zig = parameterLN2, 
                     ziln = parameterLN2)
  k.index = dim(parameter)[1]
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  ylim = c(0,1)
  
  
  res <- NULL
  for(i in 1:10)
  {
    
    for(j in 1:5)
    {
      for(k in 1:k.index)
      {
        result <- readRDS(paste0("stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
        result.stat <- data.frame(result$stat)
        
        tmp <- result.stat%>% mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,"i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),"effect" = as.character(result$setting$delta[4]))%>%dplyr::select("LB","LN","MAST","KW","KW-II","DESeq2","MGS","i","j","k","batch_f","effect")
        
        res <- rbind(res,tmp[1,])
      }
    }
  }
  res[which(res$batch !="no batch effect"),"MGS"] = NA
  
  res <- res %>% gather(key = "method", value = "p.value",`LB`,`LN`,`MAST`,`KW`,`KW-II`,`DESeq2`,`MGS`)
  
  res$method_f = factor(res$method,
                             levels = c("LN", "LB", "MAST", "KW", "KW-II", "DESeq2","MGS"),
                             labels = c("LN", "LB", "MAST", "KW", "KW-II", "DESeq2","MGS"))
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
  res %>%
    ggplot(aes(factor(k), p.value,fill = batch_f)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(ylim) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
    scale_x_discrete(labels=param.k) +
    scale_fill_manual(name = TeX("Batch effects ($\\kappa_\\mu, \\kappa_\\theta, \\kappa_{\\pi}$)"),
                      values=c("K1 (0, 0, 0)"  = "dodgerblue",
                               "K2 (0.5, 0.5, -0.5)"  = "darkgreen",
                               "K3 (1, 1, -1)" = "chartreuse3",
                               "K4 (0.5, -0.5, -0.5)"  = "red3",
                               "K5 (1, -1, -1)" = "tomato1")) +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("rejection rate") +
    facet_grid(cols = vars(method_f), rows = vars(effect_f), labeller = label_parsed) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  ggsave(file = paste0("../",model,"_full_size",size,".png"), p, width = 32, height=16)
  
  p
}

fullplot(80,model="zinb")
 fullplot(400,model="zinb")

 fullplot(80,model="zig")
 fullplot(400,model="zig")

fullplot(80,model="ziln")

fullplot(400,model="ziln")


####power plot####

powerplot <- function(model,size)
{
  parameter = switch(model, 
                     zinb = parameterNB2, 
                     zig = parameterLN2, 
                     ziln = parameterLN2)
  j.index <- c(1,5,3)
  k.index = c(7,9,10,12,25,27,28,30,43,45,46,48)
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  ylim = c(0,1)
  
  res3 <- NULL
  for(i in c(2,3,4,6,8))
  {
    
    for(j in j.index)
    {
      for(k in k.index)
      {
        result <- readRDS(paste0("stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
        result.stat <- data.frame(result$stat)
        
        tmp <- result.stat%>% mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,"i" = i,"j" = j,"k" = k,"batch" = as.character(result$setting$kappa[4]),"effect" = as.character(result$setting$delta[4]))%>%dplyr::select("LB","LN","MAST","KW","KW-II","DESeq2","MGS","i","j","k","batch","effect")
        
        res3 <- rbind(res3,tmp[1,])
      }
    }
  }
  res3[which(res3$batch !="no batch effect"),"MGS"] = NA
  
  res3 <- res3 %>% gather(key = "method", value = "p.value",`LB`,`LN`,`MAST`,`KW`,`KW-II`,`DESeq2`,`MGS`)
  res3$method_f = factor(res3$method,
                         levels = c("LN", "LB", "MAST", "KW", "KW-II","DESeq2","MGS"),
                         labels = c("LN", "LB", "MAST", "KW", "KW-II","DESeq2","MGS"))
  res3$batch_f = factor(res3$batch,
                        levels = c("no batch effect", 
                                   "large(+,+,-) batch effect",
                                   "large(+,-,-) batch effect"),
                        labels = c("K1 (0, 0, 0)",
                                   "K3 (1, 1, -1)",
                                   "K5 (1, -1, -1)"))
  res3$effect_f = factor(res3$effect,
                         levels = c("Effect_mu(D>H)", "Effect_theta(D>H)", "Effect_pi(D<H)",
                                    "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"), 
                         labels = c(TeX("D2 ($\\mu_D$>$\\mu_{H}$)"), 
                                    TeX("D3 ($\\theta_D$>$\\theta_{H}$)"),
                                    TeX("D4 ($\\pi_D$<$\\pi_{H}$)"),
                                    TeX("D6 ($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                    TeX("D8 ($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)")))
  res3 %>%
    ggplot(aes(factor(k), p.value,fill = batch_f)) +
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
    facet_grid(cols = vars(method_f), rows = vars(effect_f), labeller = label_parsed) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  ggsave(file = paste0("../",model,"_power_size",size,".png"), p, width = 22, height=12)
  p
  
}

powerplot(model = "ziln",size =400)

powerplot(model = "ziln",size =80)

powerplot(model = "zig",size =400)
powerplot(model = "zig",size =80)
powerplot(model = "zinb",size =400)
powerplot(model = "zinb",size =80)