source("C01.02.simulation.setup.R"); i=8;j=3;k=7;model=3;perturb=5;n=80;save.stat.only=TRUE;n.gene = 1000
for (i in 7:8) {
for (k in c(7,9,25,27,43)) {
  cat("i, k = ",i, ", ", k, "\n")
  perturb=5
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
  
  {
    
    
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
    #devtools::install_github("RGLab/MAST");
    library(MAST)
    library(coin)
    library(metagenomeSeq)
    
    # required parameters from...
    source("C01.02.simulation.setup.R")
    # sessionInfo()
    
    #parameter1; delta; kappa
    # (parameter = switch(model, 
    #                     zinb = parameter4, 
    #                     zig = parameter5, 
    #                     ziln = parameter5))
    (parameter = switch(model, 
                        zinb = parameterNB2, 
                        zig = parameterLN2, 
                        ziln = parameterLN2))
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
    
    #if (any(class(try(readRDS(paste0("output/R0201sim181201/result.", i, ".", j, ".", k,".rds")))) %in% "try-error")) {
    
    # do the tests on the ramdon ZINB distribution we created
    result <- tester.set.HD.batch(data, n.gene=n.gene, suppressWarnWagner = TRUE, # if not suppressed, slurm-out file size explodes.
                                  LB.skip = F,LN.skip = F, MAST.skip = F,
                                  KW.skip = F, Wg.skip = F, De2.skip = F, WRS.skip = (j != 1),
                                  MGS.skip = (j != 1)) # if there are batch effects skip WRS and MGS.
    result$nonzero.prop <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
    
    ### More MAST replicates (M = 10 in total)
    result.MAST <- list(coef = list(), pval = list())
    result.MAST$coef[[1]] <- result$coef[c("MAST.nonz", "MAST.zero", "MAST.glob"), ]
    result.MAST$pval[[1]] <- result$pval[c("MAST.nonz", "MAST.zero", "MAST.glob"), ]
    result.MAST$nonzero.prop[[1]] <- result$nonzero.prop
    
    message("More MAST replicates (s):\n")
    for (s in 2:M) {
      cat("MAST replicate s = ", s, "\n")
      set.seed(s*10^5 + i*10^3 + j*10^2 + k)
      data = do.call(r.sim, dat.args)
      data %<>% dplyr::filter(sampleSum > 0)
      cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
      tmp.MAST <- MAST(data)
      result.MAST$coef[[s]] <- tmp.MAST[[1]][1:3, ] #coef. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
      result.MAST$pval[[s]] <- tmp.MAST[[2]][1:3, ] #pval. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
      result.MAST$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
    }
    
    # plug in back the MAST result in the main object.
    result$MAST <- result.MAST
    
    
    #### statistics
    # 2.1 NA replacement   # This is not actually needed, but for consistency of data (btw first replicate and the rest)
    index.MAST <- result$pval %>% rownames() %>% {grep("MAST", .)} #5, 6, 7
    # step 1. getting NA addresses
    na.index = result$pval[index.MAST[3],] %>% is.na %>% which
    # step 2. replacing with nonzero model values
    result$pval[index.MAST[3], na.index] = result$pval[index.MAST[2], na.index]
    # step 3. getting NA addresses again and replace with zero model values.
    na.index = result$pval[index.MAST[3], ] %>% is.na %>% which
    result$pval[index.MAST[3], na.index] = result$pval[index.MAST[1], na.index]
    # leftovers
    # result[[i]][[j]][[k]]$pval[7,na.index] %>% length %>%"/"(n.gene) %T>% print
    
    # 2.2 MAST NA replacement for each replicate
    result$MAST$pval <- 
      lapply(result$MAST$pval, function(s) {
        # step 1. getting NA addresses
        na.index = s[3,] %>% is.na %>% which
        # step 2. replacing with nonzero model values
        s[3,na.index] = s[2,na.index]
        # step 3. getting NA addresses again and replace with zero model values.
        na.index = s[3,] %>% is.na %>% which
        s[3,na.index] = s[1,na.index]
        s
      })
    
    
    ### More MGS replicates (M = 10 in total)
    result.MGS <- list(coef = list(), pval = list())
    result.MGS$coef[[1]] <- result$coef["MGS", ]
    result.MGS$pval[[1]] <- result$pval["MGS", ]
    result.MGS$nonzero.prop[[1]] <- result$nonzero.prop
    
    message("More MGS replicates (s):\n")
    for (s in 2:M) {
      cat("MGS replicate s = ", s, "\n")
      set.seed(s*10^5 + i*10^3 + j*10^2 + k)
      data = do.call(r.sim, dat.args)
      data %<>% dplyr::filter(sampleSum > 0)
      cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
      tmp.MGS <- mgs(data)
      result.MGS$coef[[s]] <- tmp.MGS[, "Estimate"] #coef.
      result.MGS$pval[[s]] <- tmp.MGS[, "pval"] #pval.
      result.MGS$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
    }
    
    result$MGS <- result.MGS
    #### statistics
    # 2.1 NA replacement   # This is not actually needed, but for consistency of data (btw first replicate and the rest)
    # index.MGS <- result$pval %>% rownames() %>% {grep("MGS", .)} #5, 6, 7
    
    # 2.2 MGS NA replacement for each replicate
    # result$MGS$pval <- 
    #   lapply(result$MGS$pval, function(s) {
    #     # step 1. getting NA addresses
    #     na.index = s %>% is.na %>% which
    #     # step 2. replacing with nonzero model values
    #     s[na.index] = s[2,na.index]
    #     # step 3. getting NA addresses again and replace with zero model values.
    #     na.index = s[3,] %>% is.na %>% which
    #     s[3,na.index] = s[1,na.index]
    #     s
    #   })
    
    # 3. Statistics
    # 3.1 statistics for all method (including first replicate of MAST, which is redundant)
    index.regular <- result$nonzero.prop >= regular
    #result$pval[index.twopart, !index.regular] <- NA    # removing irregular genes for two-part models only.
    result$pval[, !index.regular] <- NA    # removing irregular genes.
    
    
    stat.power = result$pval %>% 
      apply(1, function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))})
    stat.irregular = result$pval %>% 
      apply(1, function(x) {ifelse(sum(is.na(x))/length(x)>0.9, NA, mean(ifelse(is.na(x), 1, x) <= sig, na.rm=TRUE))})
    stat.na.prop = result$pval %>% apply(1, function(x) {mean(is.na(x))})
    
    # 3.2 statistics for MAST replicates
    for (m in 1:M) {
      index.regular.m <- result$MAST$nonzero.prop[[m]] >= regular
      result$MAST$pval[[m]][, !index.regular.m] <- NA
    }
    
    stat.power[c("MAST.nonz", "MAST.zero", "MAST.glob")] <-
      result$MAST$pval %>% 
      sapply(function(s) apply(s, 1, function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, 
                                                         mean(x<=sig, na.rm=TRUE))})) %>%
      apply(1, mean, na.rm = TRUE)  # vector of three
    
    stat.irregular[c("MAST.nonz", "MAST.zero", "MAST.glob")] <- #stat based on NA p-val = 1.
      result$MAST$pval %>% 
      sapply(function(s) apply(s, 1, function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)})) %>%
      apply(1, mean, na.rm = TRUE)
    
    stat.na.prop[c("MAST.nonz", "MAST.zero", "MAST.glob")] <-
      result$MAST$pval %>% 
      sapply(function(s) apply(s, 1, function(x)  {mean(is.na(x))})) %>%
      apply(1, mean, na.rm = TRUE)
    gc()
    
    
    # 3.3 statistics for MGS replicates
    for (m in 1:M) {
      index.regular.m <- result$MGS$nonzero.prop[[m]] >= regular
      result$MGS$pval[[m]][!index.regular.m] <- NA
    }
    
    stat.power["MGS"] <-
      result$MGS$pval %>% 
      sapply( function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))}) %>%
      mean(na.rm = TRUE)  # vector of three
    
    stat.irregular["MGS"] <- #stat based on NA p-val = 1.
      result$MGS$pval %>% 
      sapply(function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)}) %>%
      mean(na.rm = TRUE)
    
    stat.na.prop["MGS"] <-
      result$MGS$pval %>% 
      sapply( function(x)  {mean(is.na(x))}) %>%
      mean(na.rm = TRUE)
    gc()
    # Save the result
    stat.comb <- rbind(stat.power, stat.irregular, stat.na.prop)
    
    saveRDS(list(stat = stat.comb, setting = setting.summary), save_file.stat)
    if (!save.stat.only) saveRDS(result, save_file.raw)
    
    
    tt(2)
  }
}
}


a <- function(model,size){
  parameter = switch(model, 
                     zinb = parameterNB2, 
                     zig = parameterLN2, 
                     ziln = parameterLN2)
  j.index <- c(1,5,3)
  k.index = c(7,9,10,12,25,27,28,30,43,45,46,48)
  param.k = apply(parameter[k.index,-1], 1, function(x) paste0("(", paste(x, collapse=", "), ")"))
  ylim = c(0,1)
  for(i in c(2,3,4,6,8)) {
    
    for(j in j.index)
    {
      for(k in k.index)
      {
        if (!file.exists(paste0("output/stat-n",size,"-pert0.5-",model,"-",i,".",j,".",k,".rds")))
          cat("i= ",i,"j= ",j,"k= ",k, "\n")
        
      }
    }
  }
}
a(model = "ziln",size =400)

