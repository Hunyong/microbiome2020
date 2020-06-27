### experiment - logistic.

# args = c(1, 1, 43, 2, 5, 80, 0, 1000)
args = commandArgs(trailingOnly=TRUE)  # passed from script
cat("The Command Arg is: ", str(args))
# i = as.numeric(args[1])  # 1..10    delta effect
j = as.numeric(args[2])  # 1..5     kappa effect
k = as.numeric(args[3])  # 1..34    baseline scenario
model = as.numeric(args[4])  # 1..3 generative model
perturb = as.numeric(args[5]) # 5, 3, 0
n = as.numeric(args[6])  # 80 800  sample size
save.stat.only = as.logical(args[7]) # 1, 0
n.gene = as.numeric(args[8]) # 1000 

if (is.na(save.stat.only)) save.stat.only = TRUE
if (is.na(n.gene)) n.gene = 1000

if (n == 80) {
  regular = 0.1 
} else if (n == 400) {
  regular = 0.02
}

if (is.null(model) | model == 1) {
  model = "zinb"
} else if (model == 2) {
  model = "zig"
} else if (model == 3) {
  model = "ziln"
}
n.sample = rep(round(n/4), 4) # sample size for H1, D1, H2, D2
n = sum(n.sample)             # correct n in case it is not a multiple of four.

if (is.null(perturb) | perturb == 0) {
  perturb = 0
} else {
  perturb = perturb/10
}

# for (i in 1:10) {
i = 1  
  
  cat("i: ", i,", j: ",j,", k: ",k,", model: ", model, ", perturb: ", perturb, "\n")
  cat("n = ", n,", stat.stat.only : ", save.stat.only,", n.gene: ",n.gene, "\n")
  
  
  ## To save the result
  # Check and create the folder
  save_path = paste0("output/")
  save_file.raw = paste0("output/raw-n", n, "-pert", perturb, "-", model, "-", i, ".", j, ".", k, ".rds")
  save_file.stat = paste0("output/stat-n", n, "-pert", perturb, "-", model, "-", i, ".", j, ".", k, ".rds")
  
  if (!dir.exists(save_path)) {message("No output folder detected. Creating one."); dir.create(save_path)}
  if (file.exists(save_file.stat)) {
    next  # next i
    # stop("done already.")
  }
  
  ### 0.1 library
  library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
  source("F00.00.generic.R")
  source("F01.01.base.R")
  source("F02.01.simulation.R")
  source("F01.02.models-base.R")
  source("F01.02.models.R")
  source("F01.02.summary.gamlss2.R")
  # devtools::install_github("RGLab/MAST");
  library(MAST)
  library(coin)
  
  # required parameters from...
  source("C01.02.simulation.setup.R")
  # sessionInfo()
  
  #parameter1; delta; kappa
  #(parameter = parameter3); 
  #(parameter = parameter4); 
  (parameter = switch(model, 
                      zinb = parameter4, 
                      zig = parameter5, 
                      ziln = parameter5))
  (delta = delta1)
  (kappa = kappa1)
  
  n.gene; n.sample; 
  print(test.dim <- method.stat %>% length)
  
  #
  
  
  tt(1)
  set.seed(i*10^3 + j*10^2 + k)
  # 1. parameter
  
  ## param.set = param (i, j, k) # list of H1, D1, H2, D2 (status-batch)
  param.set = param (i, j, k, baseParam = parameter, delta.table = delta, kappa.table = kappa)
  dat.args = list(n.sample = n.sample, n.gene=n.gene, 
                  scenario.delta = i, scenario.kappa = j, scenario.base = k,
                  baseParam = parameter,
                  delta.table = delta, 
                  kappa.table = kappa,
                  model = model, 
                  delta.perturb.prob = perturb)
  
  setting.summary <- 
    list(scenario = c(i = i, j = j, k = k),
         model = model,
         delta = delta[i, ], 
         kappa = kappa[j, ],
         baseParam = parameter[k, ],
         param.set = param.set,
         n.sample = n.sample, n.gene=n.gene,
         perturb = perturb,
         threshold = c(regular = cutoff, sig = sig))
  
  # 2. data
  data = do.call(r.sim, dat.args)
  data %<>% dplyr::filter(sampleSum > 0)
  cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
  
  #  logistic
  logistic <- 
    sapply(1:100, function(l) {
      dat.tmp = data.frame(y = data[, l], data[, 1001:1003])
      lmout = glm(y>0 ~ phenotype + batch, family = binomial, data = dat.tmp)
      summary(lmout)$coefficients["phenotypeH", c("z value", "Pr(>|z|)")]
    })
  logistic <- logistic %>% t %>% as.data.frame
  names(logistic) <- c("z", "p")
    ggplot(logistic) + geom_density(aes(z))
    ggplot(logistic) + geom_density(aes(p))
    mean(logistic$p <= 0.05, na.rm=T)
  
  #  LB
    LB2 <- function(data) {
      require(gamlss)
      n = dim(data)[1]
      data$y.prop <- data$y / data$sampleSum  #sampleSum[j] = sum_g y_g,j (g:gene, j:sample)
      bereg = try(gamlss(y.prop ~ phenotype + batch, nu.formula = ~ phenotype + batch,
                         family = BEZI(sigma.link = "log"), data = data,
                         control = gamlss.control(n.cyc = 1000, trace = FALSE)))
# bereg <<- bereg
      if (any(class(bereg) %in% "try-error")) return(out)
      
      bereg.summary <- suppressWarnings(summary.gamlss2(bereg))
      pheno.index <- grep("phenotype", rownames(bereg.summary))
      tab.tmp <- bereg.summary[pheno.index, c("t value", "Pr(>|t|)")] %>% as.matrix %>% as.vector
      names(tab.tmp) <- c("t.cnt", "t.disc", "p.cnt", "p.disc")
      tab.tmp
    }
    
    lb <- 
      sapply(1:100, function(l) {
        if (!l %% 10) cat(l, " ")
        dat.tmp = data.frame(y = data[, l], data[, 1001:1003])
        dat.tmp$y.prop <- dat.tmp$y / dat.tmp$sampleSum  #sampleSum[j] = sum_g y_g,j (g:gene, j:sample)
        lbout = try(LB2(dat.tmp))
        if (class(lbout)[1] == "try-error") rep(NA, 4) else lbout
      })
    
    lb <- lb %>% t %>% as.data.frame
    ggplot(lb) + geom_density(aes(t.disc))
    ggsave("figure/diagnostic_lb_disc.png")
    ggplot(lb) + geom_density(aes(p.disc))
    mean(lb$p.disc <= 0.05, na.rm=T)
    
    ggplot(lb) + geom_density(aes(t.cnt))
    ggsave("figure/diagnostic_lb_cnt.png")
    ggplot(lb) + geom_density(aes(t.cnt)) + xlim(c(-2,2))
    ggplot(lb) + geom_density(aes(p.cnt))
    mean(lb$p.cnt <= 0.05, na.rm=T)
    
    lmout <- 
      sapply(1:100, function(l) {
        dat.tmp = data.frame(y = data[, l], data[, 1001:1003])
        lmout = glm(y>0 ~ phenotype + batch, data = dat.tmp)
        summary(lmout)$coefficients["phenotypeH", c("t value", "Pr(>|t|)")]
      })
    lmout <- lmout %>% t %>% as.data.frame
    names(lmout) <- c("z", "p")
    ggplot(lmout) + geom_density(aes(z))
    ggsave("figure/diagnostic_lm.png")
    ggplot(lmout) + geom_density(aes(p))
    mean(lmout$p <= 0.05, na.rm=T)
    