library(tidyverse)
library(latex2exp)
library(ggh4x) # for nested facetting
source("C01.02.simulation.setup.R")


#### power plot####

metricsplot <- function(model, size, width = 12, height = 8, metrics.name = "type1error", delta.base = TRUE, qval = FALSE,
                        fn = paste0(
                          "figure/", model, "_", metrics.name, "_size", size, if (!delta.base) "_effectSize(no_batch)",
                          if (include.null) "_with_null", if (qval) "_qval", ".pdf"
                        ),
                        res.tmp = TRUE, include.null = FALSE, stop.if.absent = TRUE) {
  require(tidyr)
  parameter <- switch(model,
    zinb = parameterNB,
    zig = parameterLN,
    ziln = parameterLN
  )
  k.index <- k.core
  param.k <- apply(parameter[k.index, -1], 1, function(x) paste0("(", paste(x, collapse = ",  "), ")"))

  if (delta.base) {
    i.rng <- c(2, 3, 4, 6, 8)
    disease.levels <- c(
      "Effect_mu(D>H)", "Effect_theta(D>H)", "Effect_pi(D<H)",
      "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"
    )
    disease.labels <- c(
      TeX(r"(D2 $\mu_D$>$\mu_{H}$)"),
      TeX(r"(D3 $\theta_D$>$\theta_{H}$)"),
      TeX(r"(D4 $\pi_D$<$\pi_{H}$)"),
      TeX(r"(D6 $\mu_D$>$\mu_H$, $\pi_D$<$\pi_{H}$)"),
      TeX(r"(D8 $\mu_D$>$\mu_H$, $\pi_D$>$\pi_{H}$)")
    )
    j.rng <- c(1, 5, 3)
    batch.levels <- c("no batch effect", "large(+,+,-) batch effect", "large(+,-,-) batch effect")
    batch.labels <- c("K1 (0, 0, 0)", "K3 (1, 1, -1)", "K5 (1, -1, -1)")
  } else {
    i.rng <- c(11, 13, 15, 2, 4, 6, 12, 14, 16)
    disease.levels <- c(
      "Effect_mu_small", "Effect_mu(D>H)", "Effect_mu_large",
      "Effect_pi_small", "Effect_pi(D<H)", "Effect_pi_large",
      "Effect_mu.pi_small", "Effect_mu(D>H).pi(D<H)", "Effect_mu.pi_large"
    )
    disease.labels <- c(
      TeX("D11. $\\delta = (0.5, 0, 0)$"),
      TeX("D2. $\\delta = (1, 0, 0)$"),
      TeX("D13. $\\delta = (2, 0, 0)$"),
      TeX("D13. $\\delta = (0, 0, -0.5)$"),
      TeX("D4. $\\delta = (0, 0, -1)$"),
      TeX("D14. $\\delta = (0, 0, -2)$"),
      TeX("D15. $\\delta = (0.5, 0, -0.5)$"),
      TeX("D6. $\\delta = (1, 0, -1)$"),
      TeX("D16. $\\delta = (2, 0, -2)$")
    )
    disease2.labels <- c(
      TeX("$\\mu$ effect (D2 and its variants)"),
      TeX("$\\mu$ effect (D2 and its variants)"),
      TeX("$\\mu$ effect (D2 and its variants)"),
      TeX("$\\pi$ effect (D4 and its variants)"),
      TeX("$\\pi$ effect (D4 and its variants)"),
      TeX("$\\pi$ effect (D4 and its variants)"),
      TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
      TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
      TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)")
    )
    j.rng <- 1
    batch.levels <- c("no batch effect")
    batch.labels <- c("K1 (0, 0, 0)")
  }
  if (include.null) {
    i.rng <- c(1, i.rng)
    disease.levels <- c("Effect_null", disease.levels)
    disease.labels <- c("D1 (null)", disease.labels)
  }
  dict <- list("sensitivity" = 1, "type1error" = 2, "FDR" = 3, "accuracy" = 4, "AUC" = 5)

  res <- NULL
  for (i in i.rng)
  {
    for (j in j.rng)
    {
      cat("\ni = ", i, "j = ", j, "k = ")
      for (k in k.index)
      {
        cat(k, " ")
        fn.tmp <- paste0("output/stat-n", size, "-pert0.5-signal0.1-", model, "-", i, ".", j, ".", k, ".rds")
        if (file.exists(fn.tmp)) {
          result <- readRDS(fn.tmp)
          if (qval) result.metrics <- data.frame(result$metrics[["qval"]] %>% t()); else result.metrics <- data.frame(result$metrics[["pval"]] %>% t())
          

          tmp <-
            result.metrics %>%
            mutate(
              "LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
              "i" = i, "j" = j, "k" = k,
              "batch" = as.character(result$setting$kappa[4]),
              "effect" = as.character(result$setting$delta[4])
            ) %>%
            dplyr::select(
              "LB", "LN", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS",
              "ANCOM", "i", "j", "k", "batch", "effect"
            )

          res <- rbind(res, tmp[dict[[metrics.name]], ])
        } else {
          if (stop.if.absent) stop("Not available")
          cat("(Not available) ")
        }
      }
    }
  }

  res <- res %>%
    gather(
      key = "method", value = "metrics",
      `LB`, `LN`, `MAST`, `KW`, `KW-II`, `DS2`, `DS2ZI`, `MGS`, `ANCOM`
    )
  res$method_f <- factor(res$method,
    levels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II"),
    labels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II")
  )
  res$batch_f <- factor(res$batch, levels = batch.levels, labels = batch.labels)
  res$effect_f <- factor(res$effect, levels = disease.levels, labels = disease.labels)
  if (!delta.base) res$effect2_f <- factor(res$effect, levels = disease.levels, labels = disease2.labels)
  # res[res$method %in% c("MGS", "ANCOM", "ANCOM.sz")
  #     & res$j != 1, "p.value"] <- NA #NA for MGS, ANCOM, and ANCOM.sz with batch effects
  res$k <- factor(res$k)
  if (res.tmp) res.tmp <<- res

  res %>%
    ggplot(aes(k, metrics, fill = batch_f)) +
    geom_bar(stat = "identity", position = position_dodge(width = .8)) +
    geom_hline(yintercept = 0.05, col = "black", linetype = 2) +
    ylim(c(0, 1)) +
    geom_point(aes(col = batch_f), position = position_dodge(width = .8), shape = 15, size = 0.1) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(labels = param.k) +
    {
      if (delta.base) {
        scale_fill_manual(
          name = TeX(r"(Batch effects $\kappa_\mu, \kappa_\theta, \kappa_{\pi}$)"),
          values = c(
            "K1 (0, 0, 0)" = "dodgerblue",
            "K3 (1, 1, -1)" = "chartreuse3",
            "K5 (1, -1, -1)" = "tomato1"
          )
        )
      } else {
        scale_fill_manual(values = c("K1 (0, 0, 0)" = "dodgerblue"))
      }
    } +
    {
      if (delta.base) {
        scale_color_manual(values = c(
          "K1 (0, 0, 0)" = "dodgerblue",
          "K3 (1, 1, -1)" = "chartreuse3",
          "K5 (1, -1, -1)" = "tomato1"
        ))
      } else {
        scale_color_manual(values = c("K1 (0, 0, 0)" = "dodgerblue"))
      }
    } +
    {
      if (!delta.base) guides(fill = FALSE, color = FALSE)
    } +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab(metrics.name) +
    guides(col = FALSE) +
    {
      if (!delta.base) facet_nested(method_f ~ effect2_f + effect_f, labeller = label_parsed)
    } +
    {
      if (delta.base) facet_grid(rows = vars(method_f), cols = vars(effect_f), labeller = label_parsed)
    } +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") -> p
  ggsave(file = fn, p, width = width, height = height)
  p
}


# modify to add the head points!!!!

metricsplot(model = "ziln", size = 400, stop.if.absent = FALSE, metrics.name = "sensitivity")
metricsplot(model = "ziln", size = 400, stop.if.absent = FALSE, metrics.name = "type1error")
metricsplot(model = "ziln", size = 400, stop.if.absent = FALSE, metrics.name = "FDR")
metricsplot(model = "ziln", size = 400, stop.if.absent = FALSE, metrics.name = "accuracy")
metricsplot(model = "ziln", size = 400, stop.if.absent = FALSE, metrics.name = "AUC")

metricsplot(model = "ziln", size = 80, stop.if.absent = FALSE, metrics.name = "sensitivity")
metricsplot(model = "ziln", size = 80, stop.if.absent = FALSE, metrics.name = "type1error")
metricsplot(model = "ziln", size = 80, stop.if.absent = FALSE, metrics.name = "FDR")
metricsplot(model = "ziln", size = 80, stop.if.absent = FALSE, metrics.name = "accuracy")
metricsplot(model = "ziln", size = 80, stop.if.absent = FALSE, metrics.name = "AUC")



metricsplot_single_effect <- function(model, size, width = 12, height = 8, metrics.c = c("sensitivity", "type1error", "FDR", "accuracy", "AUC"), input.effect = "Effect_mu(D>H)", delta.base = TRUE,
                                      fn = paste0(
                                        "figure/", model, "_", input.effect, "_size", size, if (!delta.base) "_effectSize(no_batch)",
                                        if (include.null) "_with_null", ".pdf"
                                      ),
                                      res.tmp = TRUE, include.null = FALSE, stop.if.absent = TRUE) {
  require(tidyr)
  parameter <- switch(model,
    zinb = parameterNB,
    zig = parameterLN,
    ziln = parameterLN
  )
  k.index <- k.core
  param.k <- apply(parameter[k.index, -1], 1, function(x) paste0("(", paste(x, collapse = ",  "), ")"))

  if (delta.base) {
    i.rng <- c(2, 3, 4, 6, 8)
    disease.levels <- c(
      "Effect_mu(D>H)", "Effect_theta(D>H)", "Effect_pi(D<H)",
      "Effect_mu(D>H).pi(D<H)", "Effect_mu(D>H),pi(D>H)"
    )
    disease.labels <- c(
      TeX(r"(D2 $\mu_D$>$\mu_{H}$)"),
      TeX(r"(D3 $\theta_D$>$\theta_{H}$)"),
      TeX(r"(D4 $\pi_D$<$\pi_{H}$)"),
      TeX(r"(D6 $\mu_D$>$\mu_H$, $\pi_D$<$\pi_{H}$)"),
      TeX(r"(D8 $\mu_D$>$\mu_H$, $\pi_D$>$\pi_{H}$)")
    )
    j.rng <- c(1, 5, 3)
    batch.levels <- c("no batch effect", "large(+,+,-) batch effect", "large(+,-,-) batch effect")
    batch.labels <- c("K1 (0, 0, 0)", "K3 (1, 1, -1)", "K5 (1, -1, -1)")
  } else {
    i.rng <- c(11, 13, 15, 2, 4, 6, 12, 14, 16)
    disease.levels <- c(
      "Effect_mu_small", "Effect_mu(D>H)", "Effect_mu_large",
      "Effect_pi_small", "Effect_pi(D<H)", "Effect_pi_large",
      "Effect_mu.pi_small", "Effect_mu(D>H).pi(D<H)", "Effect_mu.pi_large"
    )
    disease.labels <- c(
      TeX("D11. $\\delta = (0.5, 0, 0)$"),
      TeX("D2. $\\delta = (1, 0, 0)$"),
      TeX("D13. $\\delta = (2, 0, 0)$"),
      TeX("D13. $\\delta = (0, 0, -0.5)$"),
      TeX("D4. $\\delta = (0, 0, -1)$"),
      TeX("D14. $\\delta = (0, 0, -2)$"),
      TeX("D15. $\\delta = (0.5, 0, -0.5)$"),
      TeX("D6. $\\delta = (1, 0, -1)$"),
      TeX("D16. $\\delta = (2, 0, -2)$")
    )
    disease2.labels <- c(
      TeX("$\\mu$ effect (D2 and its variants)"),
      TeX("$\\mu$ effect (D2 and its variants)"),
      TeX("$\\mu$ effect (D2 and its variants)"),
      TeX("$\\pi$ effect (D4 and its variants)"),
      TeX("$\\pi$ effect (D4 and its variants)"),
      TeX("$\\pi$ effect (D4 and its variants)"),
      TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
      TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)"),
      TeX("$\\mu$ & $\\pi$ effect (D6 and its variants)")
    )
    j.rng <- 1
    batch.levels <- c("no batch effect")
    batch.labels <- c("K1 (0, 0, 0)")
  }
  if (include.null) {
    i.rng <- c(1, i.rng)
    disease.levels <- c("Effect_null", disease.levels)
    disease.labels <- c("D1 (null)", disease.labels)
  }
  dict <- list("sensitivity" = 1, "type1error" = 2, "FDR" = 3, "accuracy" = 4, "AUC" = 5)

  res <- NULL
  for (i in i.rng)
  {
    for (j in j.rng)
    {
      cat("\ni = ", i, "j = ", j, "k = ")
      for (k in k.index)
      {
        cat(k, " ")
        fn.tmp <- paste0("output/stat-n", size, "-pert0.5-signal0.1-", model, "-", i, ".", j, ".", k, ".rds")
        if (file.exists(fn.tmp)) {
          result <- readRDS(fn.tmp)
          if (qval) result.metrics <- data.frame(result$metrics[["qval"]] %>% t()); else result.metrics <- data.frame(result$metrics[["pval"]] %>% t())

          tmp <-
            result.metrics %>%
            mutate(
              "LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob,
              "i" = i, "j" = j, "k" = k,
              "batch" = as.character(result$setting$kappa[4]),
              "effect" = as.character(result$setting$delta[4]),
            ) %>%
            dplyr::select(
              "LB", "LN", "MAST", "KW", "KW-II", "DS2", "DS2ZI", "MGS",
              "ANCOM", "i", "j", "k", "batch", "effect"
            )

          for (metrics in metrics.c) {
            tmp.line <- data.frame(tmp[dict[[metrics]], ], "metrics" = metrics) %>% rename("KW-II" = KW.II)
            res <- rbind(res, tmp.line)
          }
        } else {
          if (stop.if.absent) stop("Not available")
          cat("(Not available) ")
        }
      }
    }
  }

  res <- res %>%
    gather(
      key = "method", value = "value",
      `LB`, `LN`, `MAST`, `KW`, `KW-II`, `DS2`, `DS2ZI`, `MGS`, `ANCOM`
    ) %>%
    filter(effect == input.effect)
  res$method_f <- factor(res$method,
    levels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II"),
    labels = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II")
  )
  res$batch_f <- factor(res$batch, levels = batch.levels, labels = batch.labels)
  res$effect_f <- factor(res$effect, levels = disease.levels, labels = disease.labels)
  res$metrics_f <- factor(res$metrics, levels = metrics.c, labels = metrics.c)
  if (!delta.base) res$effect2_f <- factor(res$effect, levels = disease.levels, labels = disease2.labels)
  # res[res$method %in% c("MGS", "ANCOM", "ANCOM.sz")
  #     & res$j != 1, "p.value"] <- NA #NA for MGS, ANCOM, and ANCOM.sz with batch effects
  res$k <- factor(res$k)
  if (res.tmp) res.tmp <<- res

  res %>%
    ggplot(aes(k, value, fill = batch_f)) +
    geom_bar(stat = "identity", position = position_dodge(width = .8)) +
    geom_hline(yintercept = 0.05, col = "black", linetype = 2) +
    ylim(c(0, 1)) +
    # geom_point(aes(col = batch_f), position = position_dodge(width = .8), shape = 15, size = 0.1) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(labels = param.k) +
    {
      if (delta.base) {
        scale_fill_manual(
          name = TeX(r"(Batch effects $\kappa_\mu, \kappa_\theta, \kappa_{\pi}$)"),
          values = c(
            "K1 (0, 0, 0)" = "dodgerblue",
            "K3 (1, 1, -1)" = "chartreuse3",
            "K5 (1, -1, -1)" = "tomato1"
          )
        )
      } else {
        scale_fill_manual(values = c("K1 (0, 0, 0)" = "dodgerblue"))
      }
    } +
    {
      if (delta.base) {
        scale_color_manual(values = c(
          "K1 (0, 0, 0)" = "dodgerblue",
          "K3 (1, 1, -1)" = "chartreuse3",
          "K5 (1, -1, -1)" = "tomato1"
        ))
      } else {
        scale_color_manual(values = c("K1 (0, 0, 0)" = "dodgerblue"))
      }
    } +
    {
      if (!delta.base) guides(fill = FALSE, color = FALSE)
    } +
    xlab(expression("baseline (" * mu ~ ", " * theta * ", " * pi * ")")) +
    ylab("value") +
    guides(col = FALSE) +
    {
      if (!delta.base) facet_nested(method_f ~ effect2_f + effect_f, labeller = label_parsed)
    } +
    {
      if (delta.base) facet_grid(rows = vars(method_f), cols = vars(metrics_f), scales = "free_y", labeller = label_parsed)
    } +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") -> p
  ggsave(file = fn, p, width = width, height = height)
  p
}


# modify to add the head points!!!!

metricsplot_single_effect(model = "ziln", size = 400, stop.if.absent = FALSE, input.effect = "Effect_mu(D>H)")
metricsplot_single_effect(model = "ziln", size = 400, stop.if.absent = FALSE, input.effect = "Effect_theta(D>H)")
metricsplot_single_effect(model = "ziln", size = 400, stop.if.absent = FALSE, input.effect = "Effect_pi(D<H)")
metricsplot_single_effect(model = "ziln", size = 400, stop.if.absent = FALSE, input.effect = "Effect_mu(D>H).pi(D<H)")
metricsplot_single_effect(model = "ziln", size = 400, stop.if.absent = FALSE, input.effect = "Effect_mu(D>H),pi(D>H)")

metricsplot_single_effect(model = "ziln", size = 80, stop.if.absent = FALSE, input.effect = "Effect_mu(D>H)")
metricsplot_single_effect(model = "ziln", size = 80, stop.if.absent = FALSE, input.effect = "Effect_theta(D>H)")
metricsplot_single_effect(model = "ziln", size = 80, stop.if.absent = FALSE, input.effect = "Effect_pi(D<H)")
metricsplot_single_effect(model = "ziln", size = 80, stop.if.absent = FALSE, input.effect = "Effect_mu(D>H).pi(D<H)")
metricsplot_single_effect(model = "ziln", size = 80, stop.if.absent = FALSE, input.effect = "Effect_mu(D>H),pi(D>H)")