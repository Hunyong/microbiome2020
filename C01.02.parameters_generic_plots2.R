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
n.samp = 300
n.samp.nm = paste0("samp", n.samp)

DRNA = "RNA"
DR.no = switch(DRNA, DNA = 1, RNA = 2)
DRNA.name = ifelse(DRNA == "DNA", "_DNA", "")

type = "gene"
model = "ziln"
for (type in c("gene", "genebact", "bact")) {
  # for (model in c("ziln")) {
  for (model in c("ziln", "zinb")) {
    print(type); print(model)
    a <-
      data.frame(params = rep(c("theta", "delta", "kappa"), each = 3), 
                 params2 = c("mu", "th", "pi"), 
                 Q1 = NA, Q2 = NA, Q3 = NA)
    
    ######## p1 (baseline)
    cond.est = readRDS(paste0("output/para_selection_est", DRNA.name, zoe.nm, "_", type, "_", model, "_", nrm, "_", n.samp.nm, ".rds"))
    
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
    
    for (D in c("all", "D", "H")) {
      # for (subfig in 0:1) {
        # 1. mu pi
        set.seed(1)
        cond.est %>%
          {if (D %in% c("D", "H")) {dplyr::filter(., disease == D)} else {.}} %>% 
          sample_frac(0.1) %>% 
          ggplot(aes(mu, pi)) +
          stat_density_2d(aes(fill = stat(density)), geom = 'raster', contour = FALSE) +
          scale_fill_viridis_c(direction = -1) +
          coord_cartesian(expand = FALSE) +
          # {if (subfig) {xlim(c(0, 5))} else NULL} +
          # {if (subfig) {ylim(c(0.6, 1))} else NULL} +
          xlab(TeX('$\\mu$')) + ylab(TeX('$\\pi')) +
          geom_jitter(height = 0, width = 0, shape = '.', col = 'black', alpha = 0.2) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
        
        # 2. mu theta
        set.seed(1)
        cond.est %>%
          sample_frac(0.1) %>% 
          ggplot(aes(mu, theta)) +
          stat_density_2d(aes(fill = stat(density)), geom = 'raster', contour = FALSE) +
          scale_fill_viridis_c(direction = -1) +
          coord_cartesian(expand = FALSE) +
          ylim (c(0,25)) +
          # {if (subfig) {xlim(c(0, 5))} else NULL} +
          # {if (subfig) {ylim(c(0.6, 1))} else NULL} +
          xlab(TeX('$\\mu$')) + ylab(TeX('$\\theta')) +
          geom_jitter(height = 0, width = 0, shape = '.', col = 'black', alpha = 0.2) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
      # }
    }
    
    
    
    ####### p2 (delta)
    cond.est.delta.healthy <- 
      dplyr::filter(cond.est, group %in% c("H1","H2"))
    cond.est.delta.diseased <- 
      dplyr::filter(cond.est, group %in% c("D1","D2"))
    cond.est.delta <-
      full_join(cond.est.delta.healthy, cond.est.delta.diseased, by = c("index", "batch")) %>% 
      mutate(delta_pi = qlogis(pi.x)- qlogis(pi.y), 
             delta_mu = log(mu.x) - log(mu.y), 
             delta_theta = log(theta.x) - log(theta.y),
             batch = batch) %>%
      dplyr::select("batch", "delta_pi", "delta_mu", "delta_theta") %>% 
      na.omit
    saveRDS(cond.est.delta, paste0("output/para_delta_est", DRNA.name, zoe.nm, "_", type, "_", model, "_", nrm, "_", n.samp.nm, ".rds"))
    
    # A. ranges
    a[4, c("Q1", "Q2", "Q3")] <- cond.est.delta$delta_mu %>% abs %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
    a[5, c("Q1", "Q2", "Q3")] <- cond.est.delta$delta_theta %>% abs %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
    a[6, c("Q1", "Q2", "Q3")] <- cond.est.delta$delta_pi %>% abs %>% "["(., is.finite(.)) %>% summary %>% round(1) %>% "["(c(2,3,5))#  %>% paste(collapse = ", ")
    
    lim4 = min(max(cond.est.delta$delta_mu %>% abs, na.rm = TRUE) * 1.05, a[4, "Q3"] * 5)
    lim5 = min(max(cond.est.delta$delta_theta %>% abs, na.rm = TRUE) * 1.05, a[5, "Q3"] * 5)
    lim6 = min(max(cond.est.delta$delta_pi %>% abs, na.rm = TRUE) * 1.05, a[6, "Q3"] * 5)
    
    p2 <- 
    cond.est.delta %>%
      # dplyr::filter(delta_theta < 50 & delta_mu < 5) %>%
      ggplot(aes(abs(delta_mu), abs(delta_pi), col = batch, shape = batch)) +
      geom_point() +
      scale_shape_manual(name = "batch",values=c(17, 15)) +
      stat_density_2d(aes(fill = stat(density)), geom = 'raster', contour = FALSE) +
      scale_fill_viridis_c(direction = -1) +
      xlab(TeX('$\\delta_\\mu$')) + 
      ylab(TeX('$\\delta_\\pi')) + 
      # ylim(if(model == "zinb") c(0, 5) else c(0, 5)) + xlim(c(0, 2)) +
      ylim(c(0, lim5)) + xlim(c(0, lim4)) +
      scale_y_continuous(breaks = int_breaks) +
      #ggtitle(TeX("(|$\\delta_{\\mu}$|, |$\\delta_{\\theta}$|, |$\\delta_{\\pi}$|) estimates")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
      facet_grid(rows = vars(pi_id_f), labeller = label_parsed)
    
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
    
    ggsave(paste0("figure/para_selection", DRNA.name, zoe.nm, "_", type, "_", model, "_", nrm, "_", n.samp.nm, ".png"),  p,  width = 10, height = 10)#, dpi = 200)
    #dev.off()
    
    save(cond.est, cond.est.delta, cond.est.kappa, p, p1, p2, p3, file = sprintf("output/parameter_distn%s_%s_%s_%s_%s_%s.rda", zoe.nm, model, type, DRNA, nrm, n.samp.nm))
    
    
    ####### printing out the summary.  
    cat("    \n")
    cat(type, " ", model, "\n")
    # print(a)
    print(xtable(a), include.rownames = F)
  }
}
