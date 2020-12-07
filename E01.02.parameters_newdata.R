library(dplyr); library(magrittr); library(ggplot2); library(gridExtra); library(cowplot)
library(latex2exp)
source("F00.00.generic.R")
source("F01.01.base.R")
zero.prob <- function (vec) {mean(vec == 0)}

type = "gene"
# type = "genebact"
# type = "bact"
for (type in c("genebact", "bact", "gene")) {
  print(type)
  type.full = switch(type, gene = "geneProp.marginal", 
                     genebact = "geneProp.joint", bact = "bactProp.marginal")
  nrm = "tpm5" #"rpk" "asin"
  if (nrm == "rpk") stop("RPK is not considered for this data.")
  
  ### 0.2 Data
  # Raw data of 118 subjects
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
  
  # Choose dataset and estimators according to the nrm!
  rm(data); gc()
  
  estr.ziln = function(yvec) {
    pp <- mean(yvec == 0)
    y.nz <- yvec[yvec != 0]
    
    c(mu = mean(y.nz),
      theta = var(y.nz) / mean(y.nz)^2,
      pi = pp)
  }
  estr.zinb = function(yvec) {
    ZINB.ML.time(yvec %>% round, notation="mtp", )
  }
  
  for (model in c("ziln", "zinb")) {
    cat("model = ", model, "\n")
    estr = switch(model, ziln = estr.ziln, zinb = estr.zinb)
    
    # disease groups
    grp <- unique(DataMeta$group)
    
    
    ### marginal sample (not considering batches and disease groups)
    zero.proportion <- apply(DataTPM, 1, zero.prob)
    genes.regular.index <- which(zero.proportion <= 0.9) # 29%
    
    ### Estimation 
    if(T){
      set.seed(1)
      samp <- sample(genes.regular.index, 
                     ifelse(type == "bact", length(genes.regular.index), 300))
      
      ### 2. conditional estimates
      
      cond.est <- 
        lapply(grp, function(g) {
          cat("group ", g, "\n")
          sapply(samp, function(x) {
            # if (x %% 10) cat(x, " ")
            estr(DataTPM[x, DataMeta$group == g])
          }) %>% t %>% 
            as.data.frame %>% 
            mutate(group = g, disease = substr(group, 1, 1), batch = substr(group, 2, 2))
        }) 
      cond.est <- do.call(rbind, cond.est)
      names(cond.est) <- gsub("\\.\\(Intercept\\)", "", names(cond.est))
      if (grepl("V", names(cond.est)[1])) names(cond.est)[1:3] = c("mu", "theta", "pi")
      saveRDS(cond.est, paste0("output/para_selection_est_NEWdata_", type, "_", model, "_", nrm,".rds"))
      
      cond.est %>% 
        dplyr::filter(theta < 150) %>% 
        ggplot(aes(mu, theta, col = pi)) +
        geom_point() +
        # coord_trans(y = "log", x = "log") +
        facet_grid(disease ~ batch)
      cond.est %>%
        dplyr::filter(theta < 150) %>% 
        group_by(group) %>% 
        summarize(mu.mean = mean(mu, na.rm = TRUE),
                  theta.mean = mean(theta, na.rm = TRUE),
                  pi.mean = mean(pi, na.rm = TRUE))
      
    }
    ### Final para plot for main
    
    ######## p1 (baseline)
    cond.est.pi3 <- cond.est %>% dplyr::filter (pi > 0.27 & pi < 0.33)
    cond.est.pi6 <- cond.est %>% dplyr::filter (pi > 0.57 & pi < 0.63)
    cond.est.pi9 <- cond.est %>% dplyr::filter (pi > 0.87 & pi < 0.93)
    cond.est.pi3$pi_id = 0.3
    cond.est.pi6$pi_id = 0.6
    cond.est.pi9$pi_id = 0.9
    
    cond.est.all = rbind(cond.est.pi3, cond.est.pi6, cond.est.pi9)
    cond.est.all$group = substr(cond.est.all$group,1,2)
    
    cond.est.all$pi_id_f = factor(cond.est.all$pi_id,
                                  levels = c(0.3, 0.6, 0.9),
                                  labels = c(TeX("$\\pi \\approx 0.3$"), 
                                             TeX("$\\pi \\approx 0.6$"),
                                             TeX("$\\pi \\approx 0.9$")))
    
  if (model == "zinb") {
    lim.x = ifelse (type == "bact", 100, 50)
    lim.y = 1e+5 # ifelse (type == "genebact", 1e+5, 1e+5)
  } else {
    lim.x = 50
    lim.y = 5
  }
    
    # 1st col of parameter selection
    cond.est.all %>%
      #{if (type != "bact") {dplyr::filter(., mu < 100, theta < 100)} else {.}} %>%
      ggplot(aes(mu, theta, shape = group, col = group)) +
      geom_point() +
      # scale_color_manual(name = "group",
      #                    values=c(D1 = "skyblue", 
      #                             D2 = "deepskyblue4",  
      #                             H1 = "tomato3", 
      #                             H2 = "tomato3")) +
      # scale_shape_manual(name = "group",values=c(17, 15, 2, 0)) +
      xlab(TeX('$\\mu$')) + ylab(TeX('$\\theta')) + xlim(c(0, lim.x)) + ylim(c(0, lim.y)) +
      ggtitle("baseline parameter estimates") +
      theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
      facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p1
    
    ####### p2 (delta)
    cond.est.delta.healthy <- 
      dplyr::filter(cond.est, group %in% c("H1","H2"))
    cond.est.delta.diseased <- 
      dplyr::filter(cond.est, group %in% c("D1","D2"))
    
    cond.est.delta <- 
      cbind(h = cond.est.delta.healthy, 
            d = cond.est.delta.diseased) %>%
      mutate(delta_pi = abs(qlogis(h.pi)- qlogis(d.pi)), 
             delta_mu = abs(log(h.mu) - log(d.mu)), 
             delta_theta = abs(log(h.theta) - log(d.theta)),
             batch = h.batch) %>%
      dplyr::select("batch", "delta_pi", "delta_mu", "delta_theta")
    
    cond.est.delta.pi02 <- cond.est.delta %>% dplyr::filter (delta_pi <= 0.2)
    cond.est.delta.pi05 <- cond.est.delta %>% dplyr::filter (delta_pi > 0.4 & delta_pi <= 0.6)
    cond.est.delta.pi09 <- cond.est.delta %>% dplyr::filter (delta_pi > 0.8 & delta_pi <= 1.0)
    cond.est.delta.pi02$pi_id = 0.2
    cond.est.delta.pi05$pi_id = 0.5
    cond.est.delta.pi09$pi_id = 0.9
    
    cond.est.delta.all = rbind(cond.est.delta.pi02, 
                               cond.est.delta.pi05,
                               cond.est.delta.pi09)
    rm(cond.est.delta.pi02, cond.est.delta.pi05, cond.est.delta.pi09)
    
    cond.est.delta.all$pi_id_f = factor(cond.est.delta.all$pi_id,
                                        levels = c(0.2, 0.5, 0.9),
                                        labels = c(TeX("$\\delta_{\\pi} \\in \\[0.0, 0.2\\]$"),
                                                   TeX("$\\delta_{\\pi} \\in (0.4, 0.6\\]$"),
                                                   TeX("$\\delta_{\\pi} \\in (0.8, 1.0\\]$")))
    
    
    cond.est.delta.all %>%
      # dplyr::filter(delta_theta < 50 & delta_mu < 5) %>%
      ggplot(aes(delta_mu, delta_theta, col = batch, shape = batch)) +
      geom_point() +
      scale_shape_manual(name = "batch",values=c(17, 15)) +
      xlab(TeX('$\\delta_\\mu$')) + 
      ylab(TeX('$\\delta_\\theta')) + 
      #ylim(if(model == "zinb") c(0, 5) else c(0, 5)) + # xlim(c(0, 2)) +
      ggtitle(TeX("(|$\\delta_{\\mu}$|, |$\\delta_{\\theta}$|, |$\\delta_{\\pi}$|) estimates")) +
      theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
      facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p2
    
    
    
    ######## p3 (kappa)
    cond.est.kappa.batch1 <- 
      dplyr::filter(cond.est, group %in% c("H1","D1"))
    cond.est.kappa.batch2 <- 
      dplyr::filter(cond.est, group %in% c("H2","D2"))
    
    cond.est.kappa <- 
      cbind(b1 = cond.est.kappa.batch1, 
            b2 = cond.est.kappa.batch2) %>%
      mutate(kappa_pi = abs(qlogis(b1.pi) - qlogis(b2.pi)), 
             kappa_mu = abs(log(b1.mu) - log(b2.mu)), 
             kappa_theta = abs(log(b1.theta) - log(b2.theta)),
             disease = ifelse(b1.disease == "D", "dieased", "healthy")) %>% 
      dplyr::select("disease", "kappa_pi", "kappa_mu", "kappa_theta")
    
    
    cond.est.kappa.pi03 <- cond.est.kappa %>% dplyr::filter (kappa_pi <= 0.4)
    cond.est.kappa.pi07 <- cond.est.kappa %>% dplyr::filter (kappa_pi > 0.6 & kappa_pi <= 0.8)
    cond.est.kappa.pi12 <- cond.est.kappa %>% dplyr::filter (kappa_pi > 1.1 & kappa_pi <= 1.3)
    cond.est.kappa.pi03$pi_id = 0.3
    cond.est.kappa.pi07$pi_id = 0.7
    cond.est.kappa.pi12$pi_id = 1.2
    
    cond.est.kappa.all = rbind(cond.est.kappa.pi03, cond.est.kappa.pi07, cond.est.kappa.pi12)
    
    cond.est.kappa.all$pi_id_f = factor(cond.est.kappa.all$pi_id,
                                        levels = c(0.3, 0.7, 1.2),
                                        labels = c(TeX("$\\kappa_{\\pi} \\in \\[0.0, 0.4\\]$"), 
                                                   TeX("$\\kappa_{\\pi} \\in (0.6, 0.8\\]$"),
                                                   TeX("$\\kappa_{\\pi} \\in (1.1, 1.3\\]$")))
    
    cond.est.kappa.all %>%
      # dplyr::filter(kappa_theta < 50 & kappa_mu < 5) %>%
      ggplot(aes(kappa_mu, kappa_theta, col = disease, shape = disease)) +
      geom_point() +
      scale_shape_manual(name = "disease",values=c(16,1)) +
      xlab(TeX('$\\kappa_\\mu$')) + 
      ylab(TeX('$\\kappa_\\theta')) + 
      # ylim(if(model == "zinb") c(0, 2) else c(0, 4)) + xlim(c(0, 3)) +
      ggtitle(TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates")) +
      theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
      facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p3
    
    ## combining altogether
    plot_grid(
      plotlist = list(p1, p2, p3),
      ncol = 3, nrow = 1, labels = list("A", "B", "C")
    ) -> p
    
    ggsave(paste0("figure/para_selection_NEWdata_", type, "_", model, "_", nrm,".png"),  p,  width = 10, height = 10)#, dpi = 200)
    #dev.off()
    
    save(cond.est, cond.est.delta, cond.est.kappa, p, p1, p2, p3, file = sprintf("output/parameter_distn_NEWdata_%s_%s_%s_%s.rda", model, type, "RNA", nrm))
  }
}

