source("C01.02.simulation.setup.R")
# rounding numbers
if (FALSE) {result.stat %<>% mutate(LB.nonz = round(LB.nonz, 3), LB.zero = round(LB.zero, 3),
                                    LB.glob = round(LB.glob, 3),
                                    LN = round(LN, 3), KW = round(KW, 3),
                                    MAST.nonz = round(MAST.nonz, 3), MAST.zero = round(MAST.zero, 3),
                                    MAST.glob = round(MAST.glob, 3),
                                    Wg.nonz = round(Wg.nonz, 3), Wg.zero = round(Wg.zero, 3),
                                    Wg.glob = round(Wg.glob, 3))}

test.name = data.frame(abbr = c("LB.nonz", "LB.zero", "LB.glob", "LN", "KW",
                                "MAST.nonz", "MAST.zero", "MAST.glob",
                                "Wg.nonz", "Wg.zero", "Wg.glob", "(Reserved)"),
                       full = c("Logistic Beta - nonzero", "Logistic Beta - zero", "Logistic Beta - global", 
                                "Log-Normal", "Kruskal-Wallis",
                                "MAST - nonzero", "MAST - zero", "MAST - global",
                                "Wagner - nonzero", "Wagner - zero", "Wagner - global", "(Reserved)"))

pval.plot <- function(result.stat, parameter_in_use = parameter1,
                      i = 1, test = "LB.glob", k.index=1:dim(parameter_in_use)[1], 
                      j.index = 1:dim(kappa)[1], ylim=0:1, title = NULL
) {
  
  test.full = test.name[test == test.name[,1],2]  # full name
  delta.description = delta[i,]
  param.k = apply(parameter_in_use[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  
  if (i == (dim(delta)[1]+1)) { # legend page
    #table <- parameter1
    #p <- tableGrob(table, theme = ttheme_default(base_size=4))
    
    grid:::textGrob(paste0("Test = ", test.full, "\n",
                           ##   "\nColor = batch effect (red: large, blue: small, green: no)",
                           "\nColor = batch effect (red: no, brown: small_1, green: large_1, 
                                                    blue:small_2 , purple:large_2)",
                           "\n\nParam = 1~8 (pi=30%), 9~16 (pi=50%), 17~24 (pi=60%)",
                           "\n          25~32 & 41~43 (pi=90%),  33~40 & 44~46 (pi=95%)",
                           ##      "\n        9n + 1~3 (theta=1), 4~6 (theta=3), 7~9 (theta=10)",
                           ##      "\n        3n + 1 (mu=0.5), 2 (mu=1), 3(mu=2)",
                           ##      "\n  e.g. setting 22 = pi 90%, theta 3, mu 0.5.",
                           "\n\nEffect = mu(exp(+/-2)), theta (exp(+/-0.5), pi (expit +/- 0.5)"
    ), just = "left", x = unit(0.1, "npc"), y = unit(0.5, "npc"),
    gp=grid:::gpar(fontsize=10, col="black")) -> p
    
  } else {
    result.stat %>% 
      dplyr::filter_(.dots=paste0("i==",i)) %>%
      dplyr::filter(k %in% k.index) %>%
      dplyr::filter(j %in% j.index) %>%
      mutate_(.dots=setNames(list(as.name(test)),"p.value")) %>%
      ggplot(aes(factor(k), p.value, col=batch, group=batch, fill = batch)) +
      # geom_point() + geom_line() +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_hline(yintercept=0.05, col="black", linetype = 2) + ylim(ylim) + 
      # theme(legend.position="bottom") + 
      theme(legend.position = "none", axis.text.x = element_text(angle=90)) +
      scale_x_discrete(labels=param.k) +
      xlab(expression("parameter settings (" * mu ~ ", " * theta * ", " * pi * ")")) +
      ylab("rejection rate") +
      
      # ggtitle(paste0("Effect in: ",gsub("Effect\\_","",delta.description$detail),
      #                "\nTest: ", test.full)) -> p
      ggtitle(
        if (is.null(title)) {
          paste0("Effect in: ",gsub("Effect\\_","",delta.description$detail))
        } else {
          title
        } ) -> p
    
  }
  p
}

pval.plot.null <- function(result.stat = result.stat.na, parameter_in_use = parameter,
                           i = 1, test = c("LB", "LN", "MAST", "KW",  "KW-II"), 
                           k.index=1:dim(parameter_in_use)[1], 
                           j.index = 1:dim(kappa)[1], ylim=c(0,0.3)
) {
  
  test.full = test.name[test == test.name[,1],2]  # full name
  delta.description = delta[i,]
  param.k = apply(parameter_in_use[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  
  for(ii in 1:length(test)){
    result.stat %>% 
      dplyr::filter_(.dots=paste0("i==",i)) %>%
      dplyr::filter(k %in% k.index) %>%
      dplyr::filter(j %in% j.index) %>%
      #mutate_(.dots=setNames(list(as.name(test)),"p.value")) %>%
      mutate_(.dots=setNames(list(as.name(test[ii])),"p.value")) -> tmp
    tmp$method = test[ii]
    if(ii==1){
      tmp.full = tmp
    }else{
      tmp.full = rbind(tmp.full, tmp)
    }
  }
  
  tmp.full$method_f = factor(tmp.full$method,
                             levels = c("LN", "LB", "MAST", "KW", "KW_II"),
                             labels = c("LN", "LB", "MAST", "KW", "KW_II"))
  tmp.full$batch_f = factor(tmp.full$batch,
                            levels = c("no batch effect", 
                                       "large(+,+,-) batch effect",
                                       "large(+,-,-) batch effect"),
                            labels = c("K1 (0, 0, 0)",
                                       "K3 (1, 1, -1)",
                                       "K5 (1, -1, -1)"))
  tmp.full$size_f = factor(tmp.full$size, levels = c("n = 80",  "n = 400"))
  
  tmp.full %>%
    ggplot(aes(factor(k), p.value,fill = batch_f)) +
    #ggplot(aes(factor(k), p.value, col=batch, group=batch, fill = batch)) +
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
    facet_grid(rows = vars(method_f), cols = vars(size_f)) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")-> p
  
  
  p 
}


pval.plot.power <- function(result.stat = result.stat.na, parameter_in_use = parameter,
                            i = c(2,3,4,6,8), test = c("LB", "LN", "MAST", "KW",  "KW-II"), 
                            k.index=1:dim(parameter_in_use)[1], 
                            j.index = 1:dim(kappa)[1], ylim=0:1,
                            sample_size = 80
) {
  
  test.full = test.name[test == test.name[,1],2]  # full name
  delta.description = delta[i,]
  param.k = apply(parameter_in_use[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  effects_index = i
  
  for(iii in 1:length(effects_index)){
    for(ii in 1:length(test)){
      result.stat %>% 
        dplyr:: filter(gsub("n = ", "", size) == sample_size) %>%
        dplyr::filter_(.dots=paste0("i==",effects_index[iii])) %>%
        dplyr::filter(k %in% k.index) %>%
        dplyr::filter(j %in% j.index) %>%
        mutate_(.dots=setNames(list(as.name(test[ii])),"p.value")) -> tmp
      tmp$method = test[ii]
      if(ii==1 & iii ==1){
        tmp.full = tmp
      }else{
        tmp.full = rbind(tmp.full, tmp)
      }
    }
  }
  
  tmp.full$method_f = factor(tmp.full$method,
                             levels = c("LN", "LB", "MAST", "KW", "KW-II"),
                             labels = c("LN", "LB", "MAST", "KW", "KW-II"))
  tmp.full$effect_f = factor(tmp.full$effect,
                             levels = c("Effect_mu(D>H)", "Effect_theta(D>H)", "Effect_pi(D<H)",
                                        "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"), 
                             labels = c(TeX("D2($\\mu_D$>$\\mu_{H}$)"), 
                                        TeX("D3($\\theta_D$>$\\theta_{H}$)"),
                                        TeX("D4($\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D6($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D8($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)")))
  tmp.full$batch_f = factor(tmp.full$batch,
                            levels = c("no batch effect", 
                                       "large(+,+,-) batch effect",
                                       "large(+,-,-) batch effect"),
                            labels = c("K1 (0, 0, 0)",
                                       "K3 (1, 1, -1)",
                                       "K5 (1, -1, -1)"))
  tmp.full %>%
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
  p 
}

pval.plot.pert <- function(result.stat = result.stat.na.pert, parameter_in_use = parameter,
                            i = c(1,2,3,4,6,8), test = c("LB"), 
                            k.index=1:dim(parameter_in_use)[1], 
                            j.index = 1:dim(kappa)[1], ylim=0:1,
                            pert = c(0, 0.25, 0.5), sample_size = 80
) {
  
  test.full = test.name[test == test.name[,1],2]  # full name
  delta.description = delta[i,]
  param.k = apply(parameter_in_use[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  effects_index = i
  
  for(iii in 1:length(effects_index)){
    for(ii in 1:length(pert)){
      result.stat %>% 
        dplyr::filter(gsub("n = ", "", size) == sample_size) %>%
        dplyr::filter_(.dots=paste0("i==",effects_index[iii])) %>%
        dplyr::filter(perturbation == pert[ii]) %>%
        dplyr::filter(k %in% k.index) %>%
        dplyr::filter(j %in% j.index) %>%
        mutate_(.dots=setNames(list(as.name(test)),"p.value")) -> tmp
      tmp$method = test
      if(ii==1 & iii ==1){
        tmp.full = tmp
      }else{
        tmp.full = rbind(tmp.full, tmp)
      }
    }
  }
  
  tmp.full$perturbation_f = factor(tmp.full$perturbation,
                                   levels = c(0, 0.25, 0.5),
                                   labels = c(TeX("perutrbation=0"),
                                              TeX("perutrbation=0.25"),
                                              TeX("perutrbation=0.5")))
  
  tmp.full$effect_f = factor(tmp.full$effect,
                             levels = c("Effect_null", "Effect_mu(D>H)", 
                                        "Effect_theta(D>H)", "Effect_pi(D<H)",
                                        "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"), 
                             labels = c(TeX("D1(null)"),
                                        TeX("D2($\\mu_D$>$\\mu_{H}$)"), 
                                        TeX("D3($\\theta_D$>$\\theta_{H}$)"),
                                        TeX("D4($\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D6($\\mu_D$>$\\mu_H$, $\\pi_D$<$\\pi_{H}$)"),
                                        TeX("D8($\\mu_D$>$\\mu_H$, $\\pi_D$>$\\pi_{H}$)")))
  tmp.full$batch_f = factor(tmp.full$batch,
                            levels = c("no batch effect", 
                                       "large(+,+,-) batch effect",
                                       "large(+,-,-) batch effect"),
                            labels = c("K1 (0, 0, 0)",
                                       "K3 (1, 1, -1)",
                                       "K5 (1, -1, -1)"))
  tmp.full %>%
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
    facet_grid(cols = vars(perturbation_f), rows = vars(effect_f), labeller = label_parsed) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")  -> p
  p 
}



pval.plot.power.full <- function(result.stat = result.stat.na, parameter_in_use = parameter,
                            i = 1:10, test = c("LB", "LN", "MAST", "KW",  "KW-II"), 
                            k.index=1:dim(parameter_in_use)[1], 
                            j.index = 1:dim(kappa)[1], ylim=0:1,
                            sample_size = 80
) {
  
  test.full = test.name[test == test.name[,1],2]  # full name
  delta.description = delta[i,]
  param.k = apply(parameter_in_use[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  effects_index = i
  
  for(iii in 1:length(effects_index)){
    for(ii in 1:length(test)){
      result.stat %>% 
        dplyr:: filter(gsub("n = ", "", size) == sample_size) %>%
        dplyr::filter_(.dots=paste0("i==",effects_index[iii])) %>%
        dplyr::filter(k %in% k.index) %>%
        dplyr::filter(j %in% j.index) %>%
        mutate_(.dots=setNames(list(as.name(test[ii])),"p.value")) -> tmp
      tmp$method = test[ii]
      if(ii==1 & iii ==1){
        tmp.full = tmp
      }else{
        tmp.full = rbind(tmp.full, tmp)
      }
    }
  }
  
  tmp.full$method_f = factor(tmp.full$method,
                             levels = c("LN", "LB", "MAST", "KW", "KW-II", "DESeq2"),
                             labels = c("LN", "LB", "MAST", "KW", "KW-II", "DESeq2"))
  tmp.full$effect_f = factor(tmp.full$effect,
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
  tmp.full$batch_f = factor(tmp.full$batch,
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
  tmp.full %>%
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


if (FALSE) {
  pval.plot(result.stat, i=1, test="Wg.glob", k.index= c(2,3,5,6,11,12,14,15))
  pval.plot(result.stat, i=1, test="Wg.glob", k.index= c(2,3,5,6,11,12,14,15), j.index=c(1,3), title="Wagner test")
  pval.plot(result.stat, i=11, test="LN")
  ggsave()
  
}

## To calculate the p-value
stat_na_p_val <- function(data){
  if(sum(is.na(data))/length(data) > 0.9){
    return(NA)
  }else{
    return(mean(ifelse(is.na(data), 1, data) <= sig, na.rm=TRUE))
  }
}