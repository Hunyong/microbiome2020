cat("Microbiome - C02.01.simulation-Basic-ijk-batch.R\n")

args = commandArgs(trailingOnly=TRUE)  # passed from script
if (length(args) == 0) {
  warning("commandArgs() was not provided. Set as the default value.")
  args = c(i = 1, j = 1, k = 1, model = 1, perturb = 5, n = 80, save.stat.only = 1, n.gene = 1000)
}
cat("The Command Arg is: \n"); print(args)
# i = as.numeric(args[1])  # 1..10    delta effect
j = as.numeric(args[2])  # 1..5     kappa effect
k = as.numeric(args[3])  # 1..34    baseline scenario
model = as.numeric(args[4])  # 1..3 generative model
perturb = as.numeric(args[5]) # 5, 3, 0
n = as.numeric(args[6])  # 80 800  sample size
save.stat.only = as.logical(args[7]) # 1, 0
n.gene = as.numeric(args[8]) # 1000 
DS2.version = "zinb"

if (is.na(save.stat.only)) save.stat.only = TRUE
if (is.na(n.gene)) n.gene = 1000

if (n == 80) {
  prev.filter = 0.1 
} else if (n == 400) {
  prev.filter = 0.02
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

# for (model in c("zinb", "zig", "ziln")) {
  for (j in 1:5) {
    for (i in 1:10) {
      
      
      cat("i: ", i,", j: ",j,", k: ",k,", model: ", model, ", perturb: ", perturb, "\n")
      cat("n = ", n,", stat.stat.only : ", save.stat.only,", n.gene: ",n.gene, "\n")
      
      # bookkeeping
      # library(dplyr)
      # nm = nm1
      # nm = gsub("tmp_[0-9]*", sprintf("tmp_%s", i), nm)
      # nm = gsub("_s[0-9]*", "_s1", nm)
      # write.table(" ", nm)
      # ds.fatal = FALSE # by defalut FALSE
      # if (file.exists("tmp_bookkeeping.csv")) {
      #   book = read.csv("tmp_bookkeeping.csv") 
      #   book.s = 
      #     book %>% dplyr::filter(i == .GlobalEnv$i, j == .GlobalEnv$j, k == .GlobalEnv$k, 
      #                            model == .GlobalEnv$model, pert == .GlobalEnv$perturb, n == .GlobalEnv$n)
      # } else {
      #   book.s = data.frame(i = 0, j = 0, k = 0, model = "", pert = 0, n = 0, s = 0)
      # }
      
      ## To save the result
      # Check and create the folder
      save_path = paste0("output/")
      save_file.raw = paste0("output/raw-n", n, "-pert", perturb, "-", model, "-", i, ".", j, ".", k, ".rds")
      save_file.stat = paste0("output/stat-n", n, "-pert", perturb, "-", model, "-", i, ".", j, ".", k, ".rds")
      
      if (!dir.exists(save_path)) {message("No output folder detected. Creating one."); dir.create(save_path)}
      
      # update path and file name
      update_path = gsub("(output.*)/", "\\1_updated/", save_path)
      if (!dir.exists(update_path)) {message("No output folder detected. Creating one."); dir.create(update_path)}
      update_file.stat = gsub("(output.*)/", "\\1_updated/", save_file.stat)
      
      if (update_file.stat) {
        next
      }
      
      if (!file.exists(save_file.stat)) {
        cat("No file is present. Skipping!!!")
        saveRDS("", gsub("\\.rds", "_(NOT_DONE).rds", update_file.stat))
        next  # next i
      }
      
      ### 0.1 library
      library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
      source("F00.00.generic.R")
      source("F01.01.base.R")
      source("F02.01.simulation.R")
      source("F01.02.models-base.R")
      source("F01.02.models.R")
      source("F01.02.summary.gamlss2.R")
      library(MAST)
      library(coin)
      library(metagenomeSeq)
      
      # required parameters from...
      source("C01.02.simulation.setup.R")
      # sessionInfo()
      
      #parameter1; delta; kappa
      (parameter = switch(model, 
                          zinb = parameterNB, 
                          zig = parameterLN, 
                          ziln = parameterLN))
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
             threshold = c(regular = cutoff, sig = sig, prev.filter = prev.filter),
             DESeq.version = DS2.version)
      
      # 2. data
      data = do.call(r.sim, dat.args)
      data %<>% dplyr::filter(sampleSum > 0)
      
      # filtering
      nonzero.prop <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
      filtr = nonzero.prop >= prev.filter
      data[, which(!filtr)] = NA
      
      cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
      cat("Remaining genes after screening: ", sum(filtr), "out of ", length(filtr), ".\n")
      #if (any(class(try(readRDS(paste0("output/R0201sim181201/result.", i, ".", j, ".", k,".rds")))) %in% "try-error")) {
      
      result = readRDS(save_file.stat)
      result.tmp <-  tester.set.HD.batch(data, n.gene=n.gene, suppressWarnWagner = TRUE, # if not suppressed, the slurm-out file size explodes.
                                         LB.skip = T,LN.skip = T, MAST.skip = T,
                                         KW.skip = T, Wg.skip = T, De2.skip = F, WRS.skip = T,
                                         MGS.skip = T, 
                                         DS2.version = DS2.version)
      
      
      ### More MAST, DESeq2, MGS replicates (M = 10 in total)
      result.DS2 <- list(coef = list(), pval = list())
      result.DS2$coef[[1]] <- result.tmp$coef["DESeq2", ]
      result.DS2$pval[[1]] <- result.tmp$pval["DESeq2", ]
      result.DS2$nonzero.prop[[1]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
      # result.DS2$nonzero.prop[[1]] <- result$nonzero.prop
      
      # empty shells for the rest of the replicates
      for (s in 2:M) {
        result.DS2$coef[[s]] <- result.DS2$pval[[s]] <- rep(NA, n.gene)
      }
      message("More MAST, MGS, DESeq2 replicates (s):\n")
      for (s in 2:M) {
        cat("More replicate s = ", s, "\n")
        set.seed(s*10^5 + i*10^3 + j*10^2 + k)
        data = do.call(r.sim, dat.args)
        data %<>% dplyr::filter(sampleSum > 0)
        
        # filtering
        nonzero.prop.2 <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
        filtr.2 = nonzero.prop.2 >= prev.filter
        if (sum(filtr.2) < 10) next
        data[, which(!filtr.2)] = NA
        
        cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
        cat("Remaining genes after screening: ", sum(filtr.2), "out of ", length(filtr.2), ".\n")
        
        index.filtered = 
          apply(data[, 1:n.gene], 2, function(x) all(is.na(x))) %>% # filtered gene indices
          {which(!.)} %>% as.numeric
        index.meta = grepl("^y\\.", names(data)) %>% "!"(.) %>% which
        index.filtered.meta = c(index.filtered, index.meta)
        
        print("DESeq2")
        DS2 = switch(DS2.version, vanilla = DS2.vanilla, zinb = DS2.zinb)
        tmp.DS2 <- try({DS2(data[, index.filtered.meta])})
        if (class(tmp.DS2)[1] == "try-error") tmp.DS2 = matrix(NA, ncol = 2, dimnames = list(NULL, c("Estimate", "pval")))
        result.DS2$coef[[s]][index.filtered] <- tmp.DS2[, "Estimate"] #coef.
        result.DS2$pval[[s]][index.filtered] <- tmp.DS2[, "pval"] #pval.
        result.DS2$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
      }
      
      # plug in back the MAST result in the main object.
      result.DS2
      
      
      #### statistics
      # 2.1 NA replacement   # This is not actually needed, but for consistency of data (btw first replicate and the rest)
      
      
      # 3. Statistics
      # 3.1 statistics for all method (including first replicate of MAST, which is redundant)
      
      ## This part is not needed as it is done in the filtering step.
      # index.regular <- result$nonzero.prop >= prev.filter
      # result$pval[, !index.regular] <- NA    # removing irregular genes.
      
      # 3.3 statistics for MGS replicates
      ## This part is not needed as it is done in the filtering step.
      # for (m in 1:M) {
      #   index.regular.m <- result$MGS$nonzero.prop[[m]] >= prev.filter
      #   result$MGS$pval[[m]][!index.regular.m] <- NA
      # }
      
      result$stat["stat.power", "DESeq2"] <-
        result.DS2$pval %>% 
        sapply( function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))}) %>%
        mean(na.rm = TRUE)  # vector of three
      
      result$stat["stat.irregular", "DESeq2"] <- #stat based on NA p-val = 1.
        result.DS2$pval %>% 
        sapply(function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)}) %>%
        mean(na.rm = TRUE)
      
      result$stat["stat.na.prop", "DESeq2"] <-
        result.DS2$pval %>% 
        sapply( function(x)  {mean(is.na(x))}) %>%
        mean(na.rm = TRUE)
      gc()
      
      print(result$stat)
      
      saveRDS(result, update_file.stat)
      
      tt(2)
    }
  }
# }
