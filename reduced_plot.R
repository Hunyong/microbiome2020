###reduced plot###

library(latex2exp)
library(tidyverse)
load("micro.RData")

reducedplot <- function(model,size) {
  i = 1
  j.index <- c(1,5,3)
  k.index = c(7,9,10,12,25,27,28,30,43,45,46,48)
  parameter = switch(model, 
                     zinb = parameterNB2, 
                     zig = parameterLN2, 
                     ziln = parameterLN2)
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  res3 <- NULL
  for(j in j.index)
  {
    for(k in k.index)
    {
      result <- readRDS(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds"))
      result.stat <- data.frame(result$stat)
      
      tmp <- result.stat%>% mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
                                    "i" = i,"j" = j,"k" = k,
                                    "batch" = as.character(result$setting$kappa[4]),
                                    "effect" = as.character(result$setting$delta[4])) %>%
        dplyr::select("LB","LN","MAST","KW","KW-II","DESeq2", "MGS", "i","j","k","batch","effect")
      
      res3 <- rbind(res3,tmp[1,])
    }
  }
  # 
  res3 <- res3 %>% gather(key = "method", value = "p.value",`LB`,`LN`,`MAST`,`KW`,`KW-II`,`DESeq2`, `MGS`)
res.tmp <<- res3
  res3[res3$method == "MGS" & res3$j != 1, "p.value"] <- NA #NA for MGS with batch effects
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
  
  res3 %>%
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
  
  ggsave(file = paste0("figure/", model,"_null_size",size,".png"), p, width = 10, height=12)
  p
}


reducedplot(model = "ziln",size =400)

reducedplot(model = "ziln",size =80)

reducedplot(model = "zig",size =400)
reducedplot(model = "zig",size =80)
reducedplot(model = "zinb",size =400)
reducedplot(model = "zinb",size =80)



