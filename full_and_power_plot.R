library(tidyverse)
library(latex2exp)
load("micro.RData")

fullplot <- function(size,model)
{
  
  parameter = switch(model, 
                     zinb = parameterNB2, 
                     zig = parameterLN2, 
                     ziln = parameterLN2)
  k.index = dim(parameter)[1]
  param.k = apply(parameter[,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  ylim = c(0,1)
  
  
  res <- NULL
  for(i in 1:10)
  {
    
    for(j in 1:5)
    {
      cat("i= ",i,"j= ",j, "\n")
      
      for(k in 1:k.index)
      {
        result <- readRDS(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
        result.stat <- data.frame(result$stat)
        
        tmp <- 
          result.stat %>% 
          mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                  "i" = i,"j" = j,"k" = k,"batch_f" = as.character(result$setting$kappa[4]),
                  "effect" = as.character(result$setting$delta[4])) %>%
          dplyr::select(LB, LN, MAST, KW, `KW-II`, DESeq2, MGS, i, j, k, batch_f, effect)
        
        res <- rbind(res,tmp[1,])
      }
    }
  }
  res <- res %>% gather(key = "method", value = "p.value",`LB`,`LN`,`MAST`,`KW`,`KW-II`,`DESeq2`, `MGS`)
  
  res$method_f = factor(res$method,
                             levels = c("LN", "LB", "MAST", "KW", "KW-II", "DESeq2", "MGS"),
                             labels = c("LN", "LB", "MAST", "KW", "KW-II", "DESeq2", "MGS"))
  res$effect_f = factor(res$effect,
                             levels = c("Effect_null", "Effect_mu(D>H)", 
                                        "Effect_theta(D>H)", "Effect_pi(D<H)",
                                        "Effect_mu(D>H).theta(D>H)", "Effect_mu(D>H).pi(D<H)", 
                                        "Effect_theta(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)",
                                        "Effect_mu(D>H).theta(D<H)", "Effect_theta(D<H).pi(D<H)"), 
                             labels = c(TeX("D1(null)"),
                                        TeX("D2($\\mu_D$>$\\mu_{H}$)"), 
                                        TeX("D3($\\theta_D$>$\\theta_{H}$)"),
                                        TeX("D4($\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D5($\\mu_D$>$\\mu_H$, $\\theta_D$>$\\theta_{H}$)"),
                                        TeX("D6($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D7($\\theta_D$>$\\theta_H$, $\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D8($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)"),
                                        TeX("D9($\\mu_D$>$\\mu_H$, $\\theta_D$<$\\theta_{H}$)"),
                                        TeX("D10($\\theta_D$<$\\theta_H$, $\\pi_D$<$\\pi_{H}$)")))
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
  
  p
}

#width = 32, height = 16
pn1 <- fullplot(80,model="ziln")
ggsave(file="figure/ziln_full_80.png",pn1,width = 10,height = 7)
pn2 <- fullplot(400,model="ziln")
ggsave(file="figure/ziln_full_400.png",pn2,width = 10,height = 7)


p1 <- fullplot(80,model="zinb")
ggsave(file="figure/zinb_full_80.png",p1,width = 10,height = 7)
p2 <- fullplot(400,model="zinb")
ggsave(file="figure/zinb_full_400.png",p2,width = 10,height = 7)

pp1 <- fullplot(80,model="zig")
ggsave(file="figure/zig_full_80.png",pp1,width = 10,height = 7)
pp2 <- fullplot(400,model="zig")
ggsave(file="figure/zig_full_400.png",pp2,width = 10,height = 7)


####power plot####

powerplot <- function(model,size, width = 20, height=12, 
                      fn = paste0("figure/", model,"_power_size",size,".png"))
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
        cat("i= ",i,"j= ",j,"k= ",k, "\n")
        result <- readRDS(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
        result.stat <- data.frame(result$stat)
        
        tmp <- 
          result.stat%>% 
          mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                  "i" = i,"j" = j,"k" = k,
                  "batch" = as.character(result$setting$kappa[4]),
                  "effect" = as.character(result$setting$delta[4]))%>%
          dplyr::select("LB","LN","MAST","KW","KW-II","DESeq2","MGS", "i","j","k","batch","effect")
        
        res3 <- rbind(res3,tmp[1,])
      }
    }
  }
  res3[res3$method_f == "MGS" & res3$j != 1, "p.value"] <- NA #NA for MGS with batch effects
  
  res3 <- res3 %>% 
    gather(key = "method", value = "p.value",
           `LB`,`LN`,`MAST`,`KW`,`KW-II`,`DESeq2`, `MGS`)
  res3$method_f = factor(res3$method,
                         levels = c("LN", "LB", "MAST", "KW", "KW-II","DESeq2", "MGS"),
                         labels = c("LN", "LB", "MAST", "KW", "KW-II","DESeq2", "MGS"))
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
                         labels = c(TeX("D2($\\mu_D$>$\\mu_{H}$)"), 
                                    TeX("D3($\\theta_D$>$\\theta_{H}$)"),
                                    TeX("D4($\\pi_D$<$\\pi_{H}$)"),
                                    TeX("D6($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                    TeX("D8($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)")))
  
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
  ggsave(file = fn, p, width = width, height= height)
  p
  
}

powerplot(model = "ziln",size =400)
powerplot(model = "ziln",size =80)

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
