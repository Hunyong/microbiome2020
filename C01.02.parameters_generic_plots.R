### plots, 3d plots, and ranges
library(xtable);
library(scatterplot3d)
library(cowplot)
library(dplyr)
library(latex2exp)
library(ggplot2)
nrm = "tpm5"
zoe = 2 # zoe = "IBD"
zoe.nm = if (zoe %in% 1:2) paste0("_zoe", zoe) else "_NEWDATA"

DRNA = "RNA"
DR.no = switch(DRNA, DNA = 1, RNA = 2)
DRNA.name = ifelse(DRNA == "DNA", "_DNA", "")

type = "gene"
model = "ziln"
for (type in c("gene", "genebact", "bact")) {
  for (model in c("ziln", "zinb")) {
    print(type); print(model)
    a <-
      data.frame(params = rep(c("theta", "delta", "kappa"), each = 3), 
                 params2 = c("mu", "th", "pi"), 
                 Q1 = NA, Q2 = NA, Q3 = NA)
    
    ######## p1 (baseline)
    cond.est = readRDS(paste0("output/para_selection_est", DRNA.name, zoe.nm, "_", type, "_", model, "_", nrm, ".rds"))
      
    # A. ranges
      a[1, c("Q1", "Q2", "Q3")] <- cond.est$mu %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      a[2, c("Q1", "Q2", "Q3")] <- cond.est$theta %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      a[3, c("Q1", "Q2", "Q3")] <- cond.est$pi %>% summary %>% round(2) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      
      n.est = sum(!is.na(cond.est$mu))
      cond.est$theta[cond.est$theta > 100] = NA
      cond.est$mu[cond.est$mu > 100] = NA
      n.irregular = n.est - sum(!is.na(cond.est$theta) & !is.na(cond.est$mu))
      if (n.irregular / n.est > 0.3) 
        warning(sprintf("There are more than %2.0f non regular thetas (%2.0f%%).", 
                        n.irregular, n.irregular/n.est * 100))

      lim1 = max(cond.est$mu, na.rm = TRUE) * 1.05
      lim2 = max(cond.est$theta, na.rm = TRUE) * 1.05

    # B. 3d plots
      p1_3d <- ~{
        lvls = c("D1", "D2", "H1", "H2")
        shp = c(17, 15, 2, 0)
        shp2 = cond.est$group %>% factor(levels = lvls) %>% as.numeric() %>% "["(shp, .)
        clr = c("skyblue", "deepskyblue4", "tomato3", "tomato3")
        clr2 = cond.est$group %>% factor(levels = lvls) %>% as.numeric() %>% "["(clr, .)
        
        Q3 = a[1:3, "Q3"]
        s3d = 
          scatterplot3d(cond.est[, c("pi", "theta", "mu")], angle = 55, 
                        main = "baseline parameter estimates",
                        xlab = expression(pi),
                        ylab = expression(theta),
                        zlab = expression(mu),
                        mar = c(2.5, 2.5, 2.5, 2.5),
                        xlim = c(0, 1), # pi
                        ylim = c(0, lim2), # theta
                        zlim = c(0, lim1), # mu
                        pch = shp2, color = clr2, box = TRUE)
      }
      
      
    
    # C. sliced plots
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
      int_breaks <- function(x, n = 5) {
        l <- pretty(x, n)
        l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
      }
      
      # 1st col of parameter selection
      cond.est.all %>%
        # {if (type != "bact") {dplyr::filter(., mu < 50, theta < 50)} else {.}} %>%
        ggplot(aes(mu, theta, shape = group, col = group)) +
        geom_point() +
        scale_color_manual(name = "group",
                           values=c(D1 = "skyblue", 
                                    D2 = "deepskyblue4",  
                                    H1 = "tomato3", 
                                    H2 = "tomato3")) +
        scale_shape_manual(name = "group",values=c(17, 15, 2, 0)) +
        xlim(c(0, lim1)) + ylim(c(0, lim2)) +
        scale_y_continuous(breaks = int_breaks) +
        xlab(TeX('$\\mu$')) + ylab(TeX('$\\theta')) + # xlim(c(0, 50)) +
        #ggtitle("baseline parameter estimates") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
        facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p1
      
      p1_all = plot_grid(p1_3d, p1, ncol = 1, nrow = 2, rel_heights = c(1, 2))
      
      
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
        dplyr::select("batch", "delta_pi", "delta_mu", "delta_theta") %>% 
        na.omit
      
    # A. ranges
      a[4, c("Q1", "Q2", "Q3")] <- cond.est.delta$delta_mu %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      a[5, c("Q1", "Q2", "Q3")] <- cond.est.delta$delta_theta %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      a[6, c("Q1", "Q2", "Q3")] <- cond.est.delta$delta_pi %>% "["(., is.finite(.)) %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      
      lim4 = min(max(cond.est.delta$delta_mu, na.rm = TRUE) * 1.05, a[4, "Q3"] * 5)
      lim5 = min(max(cond.est.delta$delta_theta, na.rm = TRUE) * 1.05, a[5, "Q3"] * 5)
      lim6 = min(max(cond.est.delta$delta_pi, na.rm = TRUE) * 1.05, a[6, "Q3"] * 5)
      
      
    # B. 3d plots
      p2_3d <- ~{
        lvls = 1:2
        shp = c(17, 15)
        shp2 = cond.est.delta$batch %>% factor(levels = lvls) %>% as.numeric() %>% "["(shp, .)
        clr = c("#F8766D", "#00BFC4")
        clr2 = cond.est.delta$batch %>% factor(levels = lvls) %>% as.numeric() %>% "["(clr, .)
        
        Q3 = a[4:6, "Q3"]
        s3d = 
          scatterplot3d(cond.est.delta[, c("delta_pi", "delta_theta", "delta_mu")], angle = 55, 
                        main = "disease effect estimates",
                        xlab = expression(delta[pi]),
                        ylab = expression(delta[theta]),
                        zlab = expression(delta[mu]),
                        mar = c(2.5, 2.5, 2.5, 2.5),
                        xlim = c(0, lim6), # pi
                        ylim = c(0, lim5), # theta
                        zlim = c(0, lim4), # mu
                        pch = shp2, color = clr2, box = TRUE)
      }
      
      
      
    # C. sliced plots
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
        # ylim(if(model == "zinb") c(0, 5) else c(0, 5)) + xlim(c(0, 2)) +
        ylim(c(0, lim5)) + xlim(c(0, lim4)) +
        scale_y_continuous(breaks = int_breaks) +
        #ggtitle(TeX("(|$\\delta_{\\mu}$|, |$\\delta_{\\theta}$|, |$\\delta_{\\pi}$|) estimates")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
        facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p2
      
      p2_all = plot_grid(p2_3d, p2, ncol = 1, nrow = 2, rel_heights = c(1, 2))
      
      
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
               disease = ifelse(b1.disease == "D", "diseased", "healthy")) %>% 
        dplyr::select("disease", "kappa_pi", "kappa_mu", "kappa_theta")
      
    # A. ranges
      a[7, c("Q1", "Q2", "Q3")] <- cond.est.kappa$kappa_mu %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      a[8, c("Q1", "Q2", "Q3")] <- cond.est.kappa$kappa_theta %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      a[9, c("Q1", "Q2", "Q3")] <- cond.est.kappa$kappa_pi %>% "["(., is.finite(.)) %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
      
      lim7 = min(max(cond.est.kappa$kappa_mu, na.rm = TRUE) * 1.05, a[7, "Q3"] * 5)
      lim8 = min(max(cond.est.kappa$kappa_theta, na.rm = TRUE) * 1.05, a[8, "Q3"] * 5)
      lim9 = min(max(cond.est.kappa$kappa_pi, na.rm = TRUE) * 1.05, a[9, "Q3"] * 5)
      
    # B. 3d plots
      p3_3d <- ~{
        lvls = c("diseased", "healthy")
        shp = c(16, 1)
        shp2 = cond.est.kappa$disease %>% factor(levels = lvls) %>% as.numeric() %>% "["(shp, .)
        clr = c("#F8766D", "#00BFC4")
        clr2 = cond.est.kappa$disease %>% factor(levels = lvls) %>% as.numeric() %>% "["(clr, .)
        
        Q3 = a[7:9, "Q3"]
        s3d = 
          scatterplot3d(cond.est.kappa[, c("kappa_pi", "kappa_theta", "kappa_mu")], angle = 55, 
                        # main = "batch effect estimates",
                        main = TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates"),
                        xlab = expression(kappa[pi]),
                        ylab = expression(kappa[theta]),
                        zlab = expression(kappa[mu]),
                        mar = c(2.5, 2.5, 2.5, 2.5),
                        xlim = c(0, lim9), # pi
                        ylim = c(0, lim8), # theta
                        zlim = c(0, lim7), # mu
                        pch = shp2, color = clr2, box = TRUE)
        # legend(s3d$xyz.convert(Q3[3]*3, Q3[2]*5, Q3[1]*5), legend = lvls,
        #        col =  clr, pch = shp, box.lwd = 0.1)
      }
      
      
    # C. sliced plots
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
      
      lim.x = a[7, "Q3"] * 10
      lim.y = a[8, "Q3"] * 10
      
      
      cond.est.kappa.all %>%
        # dplyr::filter(kappa_theta < 50 & kappa_mu < 5) %>%
        ggplot(aes(kappa_mu, kappa_theta, col = disease, shape = disease)) +
        geom_point() +
        scale_shape_manual(name = "disease",values=c(16,1)) +
        xlab(TeX('$\\kappa_\\mu$')) + 
        ylab(TeX('$\\kappa_\\theta')) + 
        # ylim(if(model == "zinb") c(0, 5) else c(0, 5)) + xlim(c(0, 2)) +
        ylim(c(0, lim8)) + xlim(c(0, lim7)) +
        scale_y_continuous(breaks = int_breaks) +
        # ggtitle(TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
        facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p3
      
      p3_all = plot_grid(p3_3d, p3, ncol = 1, nrow = 2, rel_heights = c(1, 2))
      
      ## combining altogether
      plot_grid(
        plotlist = list(p1_all, p2_all, p3_all),
        ncol = 3, nrow = 1, labels = list("A", "B", "C")
      ) -> p
      
      ggsave(paste0("figure/para_selection", DRNA.name, zoe.nm, "_", type, "_", model, "_", nrm,".png"),  p,  width = 10, height = 10)#, dpi = 200)
      #dev.off()
      
      save(cond.est, cond.est.delta, cond.est.kappa, p, p1, p2, p3, file = sprintf("output/parameter_distn%s_%s_%s_%s_%s.rda", zoe.nm, model, type, DRNA, nrm))
      
      
    ####### printing out the summary.  
      cat("    \n")
      cat(type, " ", model, "\n")
      # print(a)
      print(xtable(a), include.rownames = F)
  }
}
