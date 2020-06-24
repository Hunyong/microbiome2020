library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
library(latex2exp)
source("F00.00.generic.R")
source("F01.01.base.R")
zero.prob <- function (vec) {mean(vec == 0)}


model = "ziln" #"zinb"
type = "tpm5" #"rpk"

### 0.2 Data
  # Raw data of 118 subjects
  gene.marginal.RPK.DRNA <- readRDS("../Data-processed/data.geneRPK.marginal.DRNA.ZOE1.rds")
  excluded.subject <- gene.marginal.RPK.DRNA$meta$id %in% c(352, 420)
  DataMeta116 = gene.marginal.RPK.DRNA$meta[!excluded.subject,]
  DataMeta116 <-
    DataMeta116 %>% 
    mutate(group = paste0(ifelse(cariesfree == 1, "H", "D"), ifelse(batch.RNA == "170628", 1, 2)),
           group_batch =  ifelse(batch.RNA == "170628", 1, 2),
           group_disease = ifelse(cariesfree == 1, "H", "D"))
  
  RNA         = gene.marginal.RPK.DRNA$otu[,, 2]
  DataRPK116  = RNA[,colnames(RNA) %in% DataMeta116$id]
  ST          = apply(DataRPK116, 2, sum)
  mean(ST) # 5,551,718
  DataTPM116  = t(t(DataRPK116)/ST) * 5E+6
  DataComp    = t(t(DataRPK116)/ST)  # for LB tests
  
  # Choose dataset and estimators according to the type!
  Data = switch(type, rpk = DataRPK116, tpm5 = DataTPM116)
  rm(gene.marginal.RPK.DRNA, excluded.subject, RNA, DataRPK116, DataTPM116); gc()
  
  estr.ziln = function(yvec) {
    pp <- mean(yvec == 0)
    y.nz <- yvec[yvec != 0]
    
    c(mu = mean(y.nz),
      theta = var(y.nz) / mean(y.nz)^2,
      pi = pp)
  }
  estr.zinb = function(yvec) {
    ZINB.ML.time(yvec %>% round, notation="mtp")
  }
  
  for (model in c("ziln", "zinb")) {
    cat("model = ", model, "\n")
      estr = switch(model, ziln = estr.ziln, zinb = estr.zinb)
      
      # disease groups
      grp <- unique(DataMeta116$group)
      
      
    ### marginal sample (not considering batches and disease groups)
      zero.proportion <- apply(Data, 1, zero.prob)
      genes.regular.index <- which(zero.proportion <= 0.9) # 29%
      
    ### Estimation 
    if(T){
      set.seed(1)
      samp <- sample(genes.regular.index, 300)
      
      ### 1. marginal estimates
      # stat.ziln <- 
      #   sapply(samp, function(x) {
      #     yvec <- Data[x, ]  # TPM
      #     pp <- mean(yvec == 0)
      #     y.nz <- yvec[yvec != 0]
      #     
      #     c(pi = pp,
      #       mu = mean(y.nz),
      #       sig = sd(y.nz),
      #       theta = var(y.nz) / mean(y.nz)^2)
      #   }) %>% t %>% as.data.frame
      # 
      # ggplot(stat.ziln, aes(mu, theta, col = pi)) +
      #   geom_point() 
      # 
      # stat.ziln %>%
      #   # dplyr::filter(mu < 1000) %>%
      #   ggplot(aes(mu, theta, col = pi)) +
      #   coord_trans(y = "log", x = "log") +
      #   geom_point()
      # 
      
      ### 2. conditional estimates
      
      cond.est <- 
        lapply(grp, function(g) {
          cat("group ", g, "\n")
          sapply(samp, function(x) {
            # if (x %% 10) cat(x, " ")
            estr(Data[x, DataMeta116$group == g])
          }) %>% t %>% 
            as.data.frame %>% 
            mutate(group = g, disease = substr(group, 1, 1), batch = substr(group, 2, 2))
        }) 
      cond.est <- do.call(rbind, cond.est)
      names(cond.est) <- gsub("\\.\\(Intercept\\)", "", names(cond.est))
      saveRDS(cond.est, paste0("output/para_selection_est_", model, "_", type,".rds"))
      
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
      # 1st col of parameter selection
      cond.est.all %>%
        dplyr::filter(mu < 50, theta < 50) %>%
        ggplot(aes(mu, theta, shape = group, col = group)) +
        geom_point() +
        scale_color_manual(name = "group",
                           values=c(D1 = "skyblue", 
                                    D2 = "deepskyblue4",  
                                    H1 = "tomato3", 
                                    H2 = "tomato3")) +
        scale_shape_manual(name = "group",values=c(17, 15, 2, 0)) +
        xlab(TeX('$\\mu$')) + ylab(TeX('$\\theta')) + 
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
        mutate(delta_pi = pmax(h.pi/d.pi, d.pi/h.pi), 
               delta_mu = pmax(h.mu/d.mu, d.mu/h.mu), 
               delta_theta = pmax(h.theta/d.theta, d.theta/h.theta),
               batch = h.batch) %>%
        dplyr::select("batch", "delta_pi", "delta_mu", "delta_theta")
      
      cond.est.delta.pi12 <- cond.est.delta %>% dplyr::filter (delta_pi <= 1.2)
      cond.est.delta.pi15 <- cond.est.delta %>% dplyr::filter (delta_pi > 1.2 & delta_pi <= 1.5)
      cond.est.delta.pi20 <- cond.est.delta %>% dplyr::filter (delta_pi > 1.5 & delta_pi <= 2.0)
      cond.est.delta.pi12$pi_id = 1.2
      cond.est.delta.pi15$pi_id = 1.5
      cond.est.delta.pi20$pi_id = 2.0
    
      cond.est.delta.all = rbind(cond.est.delta.pi12, 
                                  cond.est.delta.pi15,
                                  cond.est.delta.pi20)
      rm(cond.est.delta.pi10, cond.est.delta.pi15, cond.est.delta.pi20)
    
      cond.est.delta.all$pi_id_f = factor(cond.est.delta.all$pi_id,
                                           levels = c(1.2, 1.5, 2.0),
                                           labels = c(TeX("$\\delta_{\\pi} \\in \\[1.0, 1.2\\]$"),
                                                      TeX("$\\delta_{\\pi} \\in (1.2, 1.5\\]$"),
                                                      TeX("$\\delta_{\\pi} \\in (1.5, 2.0\\]$")))
    
    
      cond.est.delta.all %>%
        dplyr::filter(delta_theta < 50 & delta_mu < 5) %>%
        ggplot(aes(delta_mu, delta_theta, col = batch, shape = batch)) +
        geom_point() +
        scale_shape_manual(name = "batch",values=c(17, 15)) +
        xlab(TeX('$\\delta_\\mu$')) + 
        ylab(TeX('$\\delta_\\theta')) + 
        ylim(if(model == "zinb") c(0, 50) else c(0, 30)) + xlim(c(1, 5)) +
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
        mutate(kappa_pi = pmax(b1.pi/b2.pi, b2.pi/b1.pi), 
               kappa_mu = pmax(b1.mu/b2.mu, b2.mu/b1.mu), 
               kappa_theta = pmax(b1.theta/b2.theta, b2.theta/b1.theta),
               disease = ifelse(b1.disease == "D", "dieased", "healthy")) %>% 
        dplyr::select("disease", "kappa_pi", "kappa_mu", "kappa_theta")
      
      
      cond.est.kappa.pi12 <- cond.est.kappa %>% dplyr::filter (kappa_pi <= 1.2)
      cond.est.kappa.pi15 <- cond.est.kappa %>% dplyr::filter (kappa_pi > 1.2 & kappa_pi <= 1.5)
      cond.est.kappa.pi20 <- cond.est.kappa %>% dplyr::filter (kappa_pi > 1.5 & kappa_pi <= 2.0)
      cond.est.kappa.pi12$pi_id = 1.2
      cond.est.kappa.pi15$pi_id = 1.5
      cond.est.kappa.pi20$pi_id = 2.0
      
      cond.est.kappa.all = rbind(cond.est.kappa.pi12, cond.est.kappa.pi15, cond.est.kappa.pi20)
      
      cond.est.kappa.all$pi_id_f = factor(cond.est.kappa.all$pi_id,
                                           levels = c(1.2, 1.5, 2.0),
                                           labels = c(TeX("$\\kappa_{\\pi} \\in \\[1.0, 1.2\\]$"), 
                                                      TeX("$\\kappa_{\\pi} \\in (1.2, 1.5\\]$"),
                                                      TeX("$\\kappa_{\\pi} \\in (1.5, 2.0\\]$")))
      
      cond.est.kappa.all %>%
        dplyr::filter(kappa_theta < 50 & kappa_mu < 5) %>%
        ggplot(aes(kappa_mu, kappa_theta, col = disease, shape = disease)) +
        geom_point() +
        scale_shape_manual(name = "disease",values=c(16,1)) +
        xlab(TeX('$\\kappa_\\mu$')) + 
        ylab(TeX('$\\kappa_\\theta')) + 
        ylim(if(model == "zinb") c(0, 50) else c(0, 30)) + xlim(c(1, 5)) +
        ggtitle(TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates")) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
        facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p3
    
    ## combining altogether
      grid.arrange(
        grobs = list(p1, p2, p3),
        ncol = 3, nrow = 1
      ) -> p
    
    ggsave(paste0("figure/para_selection_", model, "_", type,".png"),  p,  width = 10, height = 10)#, dpi = 200)
    #dev.off()
  }
  