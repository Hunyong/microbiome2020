
param.beta <- function(beta.obj) {
  mu     = beta.obj$mu.coefficients %>% plogis %>% as.numeric
  sig    = beta.obj$sigma.coefficients %>% exp %>% as.numeric
  aplusb = 1/sig^2 - 1  # see ?BE
  alpha  = mu * aplusb
  beta   = aplusb - alpha
  return(c(alpha = alpha, beta = beta, mu = mu, sig = sig))
}
param.gamma <- function(gamma.obj) {
  mu = gamma.obj$mu.coefficients %>% exp %>% as.numeric
  sig = gamma.obj$sigma.coefficients %>% exp %>% as.numeric
  # alpha = mu / beta
  return(c(mu = mu, sig = sig))
}
param.ln <- function(ln.obj) {
  mu = ln.obj$mu.coefficients %>% as.numeric
  sig = ln.obj$sigma.coefficients %>% exp %>% as.numeric
  return(c(mu = mu, sig = sig))
}
param.zinb <- function(zinb.obj) {
  mu  = zinb.obj$mu.coefficients %>% exp %>% as.numeric
  sig = zinb.obj$sigma.coefficients %>% exp %>% as.numeric
  pi  = zinb.obj$nu.coefficients %>% plogis %>% as.numeric
  return(c(mu = mu, sig = sig, pi = pi))
}

param.zinb2 <- function(zinb.obj) {
  if (class(zinb.obj)[1] != "zeroinfl") stop("Designed only for zeroinfl object.")
  mu  = zinb.obj$coefficients$count %>% exp %>% as.numeric
  sig = zinb.obj$theta %>% as.numeric
  pi  = zinb.obj$coefficients$zero %>% plogis %>% as.numeric
  return(c(mu = mu, sig = sig, pi = pi))
}

con = gamlss.control(n.cyc = 100, trace = FALSE)
est.beta <- function(data) { # data is positive relative abundance (tpm.i[tpm.i > 0]/const)
  full.beta  <- gamlss(data ~ 1, family = BE(sigma.link = "log"), control = con)
  param.beta(full.beta)
}
est.gamma <- function(data) { # data is positive abundance (otu.i[otu.i > 0])
  full.gamma  <- gamlss(data ~ 1, family = GA(sigma.link = "log"), control = con)
  param.gamma(full.gamma)
}
est.ln <- function(data) { # data is positive abundance (otu.i[otu.i > 0])
  full.ln  <- gamlss(data ~ 1, family = LOGNO(sigma.link = "log"), control = con)
  param.ln(full.ln)
}

ks.empirical <- function(data, model = "beta", return.est = FALSE, 
                         ks.pval = FALSE, lilliefors = TRUE, n.lilliefors = 100) {
  require(dplyr)
  est  = switch(model, beta = est.beta, gamma = est.gamma, ln = est.ln)
  pFUN = switch(model, beta = "pbeta", gamma = "pGA", ln = "pLOGNO")
  rFUN = switch(model, beta = rbeta, gamma = rGA, ln = rLOGNO)
  
  param = est(data)
  args = 
    if (model == "beta") {
      list(shape1 = param["alpha"], shape2 = param["beta"])
    } else {
      list(mu = param["mu"], sigma = param["sig"])
    }
  ks = do.call(ks.test, c(list(x = data, y = pFUN), args))
  result = c(ks.coef = as.numeric(ks$statistic), ks.pval =  if (ks.pval) as.numeric(ks$p.value))
  
  if (lilliefors) {
    cat(n.lilliefors, " Monte Carlo samples for the Lilliefors procedure\n")
    n = length(data)
    null.set =
      sapply(1:n.lilliefors, function(s) {
        if (s %% 10 == 0) cat(s, " ")
        x.sim = do.call(rFUN, c(list(n = n), args))
        res = try(ks.empirical(x.sim, model = model, lilliefors = FALSE, return.est = FALSE, ks.pval = FALSE)["ks.coef"])
        ifelse(class(res) == "try-error", NA, res)
      })
    cat("\n")
    # print(sort(null.set))
    lil.pval = 1 - mean(result["ks.coef"] > null.set, na.rm = TRUE)
    result = c(result, lil.pval = lil.pval)
  }
  
  if (return.est) {
    return(c(param, result))
  }
  return(result)
}

if (0) {
  # usual KS test (biased)
  ks.empirical(otu.i[otu.i>0]/const, model = "gamma", return.est = TRUE, lilliefors = FALSE)
  ks.empirical(otu.i[otu.i>0]/const, model = "gamma", return.est = TRUE, lilliefors = TRUE, n.lilliefors = 250)
}



qqplot1 <- function(values, distnFn, ks.pval, ..., title = "Beta", taxa = NULL) {
  values = sort(values)
  pval   = distnFn(values, ...)
  index  = (1:length(values) - 0.5)/length(values)
  kstest = paste0("Kolmogorov-Smirnov p = ", signif(ks.pval, 2))
  taxa   = if (!is.null(taxa)) paste0(" - ", taxa) else NULL
  
  data.frame(index = -log10(index), p = -log10(pval)) %>%
    ggplot(aes(index, p)) + 
    geom_line(col = "slateblue4", size = 1, alpha = 0.5) +
    geom_point() +
    annotate("text", x = 1, y = 0.25, label = kstest) +
    # xlim(c(0,max.q.x)) + ylim(c(0,max.q.y)) +
    xlim(c(0,2)) + ylim(c(0,2)) +
    xlab(expression(paste("Expected (",-log[10], " quantile)"))) + 
    ylab(expression(paste("Observed (",-log[10], " quantile)"))) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_abline(slope = 1, intercept = 0, col = "red", linetype = "dashed", size = 0.7) +
    ggtitle(paste0("QQ-plot under ", title, " distribution", taxa)) + theme_bw()
}

qqplot.zinb <- function(values, rng = seq(min(values), max(values), by = 1), 
                        estimates, taxa = NULL) {
  values = sort(values)
  ecdf   = sapply(rng, function(s) mean(values <= s))
  cdf    = pZINBI(rng, mu = estimates["mu"], sigma = estimates["sig"], nu = estimates["pi"])
  taxa   = if (!is.null(taxa)) paste0(" - ", taxa) else NULL
  
  data.frame(value = rng, ecdf = ecdf, cdf =  cdf) %>%
    ggplot(aes(value, ecdf)) +
    geom_col() +
    geom_line(aes(value, cdf), col = "red", size = 1, alpha = 0.5) +
    # geom_point(aes(value, cdf), col = "red", size = 1) +
    xlab("counts") + 
    ylab("cumulative distribution") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggtitle(paste0("Empirical (bars) v.s. model (red line) CDF under ZINB distribution", taxa)) + theme_bw()
}

