cat("Microbiome - C02.01.simulation-Basic-ijk-batch.R\n")

### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
source("F00.00.generic.R")
source("F01.01.base.R")
source("F02.01.simulation.R")
source("F01.02.models-base.R")
source("F01.02.models.R")
source("F01.02.models-ANCOMBC.R")
source("F01.02.summary.gamlss2.R")
source("F02.01.metrics.R")

library(MAST)
library(coin)
library(metagenomeSeq)

# required parameters from...
source("C01.02.simulation.setup.R")
# sessionInfo()


args = commandArgs(trailingOnly=TRUE)  # passed from script
if (length(args) == 0) {
  warning("commandArgs() was not provided. Set as the default value.")
  args = c(i = 1, j = 1, k = 1, model = 3, perturb = 5, n = 80, save.stat.only = 0, n.gene = 1000, portion.signal = 0.1, replica = 1)
}
cat("The Command Arg is: \n"); print(args)
i = as.numeric(args[1])  # 1..10    delta effect
j = as.numeric(args[2])  # 1..5     kappa effect
k = as.numeric(args[3])  # 1..34    baseline scenario
model = as.numeric(args[4])  # 1..3 generative model
perturb = as.numeric(args[5]) # 5, 3, 0
n = as.numeric(args[6])  # 80 800  sample size
save.stat.only = as.logical(args[7]) # 1, 0
n.gene = as.numeric(args[8]) # 1000 
portion.signal = as.numeric(args[9]) # 0.1

replica = as.numeric(args[10])

# nm1 = sprintf("tmp_%s_%s_%s_%s_pert%1.1f_n%s_s%s.txt", 1, j, k, model, perturb, n, 1) # bookkeeping

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

# We do only up to j = 1, 2, ..., 10 for most of the settings. But we do it fully for some core settings (the ones in the figures of the main body of the paper)
j.terminal = ifelse(model == "ziln" & j == 1 & k %in% k.core & perturb == 0.5,
                    dim(delta1)[1], 10)  #n = both
cat("Do k upto ", j.terminal, "\n")

if (i == 0) rng = 1:j.terminal else rng = i:j.terminal
for (i in rng) {
  
    
    cat("i: ", i,", j: ",j,", k: ",k,", model: ", model, ", perturb: ", perturb, "\n")
    cat("n = ", n,", stat.stat.only : ", save.stat.only,", n.gene: ",n.gene, "portion.signal: ", portion.signal, "\n")
    cat("replica = ", replica, "\n")
    
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
    save_file.raw = paste0("output/raw-n", n, "-pert", perturb, "-signal", portion.signal, "-", model, "-", i, ".", j, ".", k,"-replica", replica, ".rds")
    save_file.stat = paste0("output/stat-n", n, "-pert", perturb, "-signal", portion.signal, "-", model, "-", i, ".", j, ".", k, "-replica", replica, ".rds")
    
    if (!dir.exists(save_path)) {message("No output folder detected. Creating one."); dir.create(save_path)}
    if (file.exists(save_file.stat)) {
      # # bookkeeping
      # if (file.exists(nm)) file.remove(nm)
      
      next  # next i
      # stop("done already.")
    }
    
    #parameter1; delta; kappa
    (parameter = switch(model, 
                        zinb = parameterNB, 
                        zig = parameterLN, 
                        ziln = parameterLN))
    (delta = delta1)
    (kappa = kappa1)
    
    n.gene; n.sample; 
    print(test.dim <- method.stat %>% length)
    cdf.cutoff = c((0:100)/500, 0.2 + 0.8 * (1:80)/80) # higher resolution for < 0.2.
    
    tt(1)
    set.seed(replica*10^4 + i*10^3 + j*10^2 + k)
    # 1. parameter
    
    ## param.set = param (i, j, k) # list of H1, D1, H2, D2 (status-batch)
    param.set = param (i, j, k, baseParam = parameter, delta.table = delta, kappa.table = kappa)
    dat.args = list(n.sample = n.sample, n.gene=n.gene, 
                    scenario.delta = i, scenario.kappa = j, scenario.base = k,
                    baseParam = parameter,
                    delta.table = delta, 
                    kappa.table = kappa,
                    model = model, 
                    delta.perturb.prob = perturb,
                    portion.signal = portion.signal)
    
    setting.summary <- 
      list(scenario = c(i = i, j = j, k = k),
           model = model,
           delta = delta[i, ], 
           kappa = kappa[j, ],
           baseParam = parameter[k, ],
           param.set = param.set,
           n.sample = n.sample, n.gene=n.gene,
           perturb = perturb,
           portion.signal = portion.signal,
           replica = replica,
           threshold = c(regular = cutoff, sig = sig, prev.filter = prev.filter))
    
    # 2. data
    data = do.call(r.sim, dat.args)
    data %<>% dplyr::filter(sampleSum > 0)
    n.signal = attr(data, "n.signal")
    
    # filtering
    prop <- apply(data[, 1:n.gene], 2, function(s) mean(round(s)))
    nonzero.prop <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
    filtr = (nonzero.prop >= prev.filter) * (prop>0)
    data[, which(!filtr)] = NA
    
    cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
    cat("Remaining genes after screening: ", sum(filtr), "out of ", length(filtr), ".\n")
    
    # do the tests on the ramdon ZINB distribution we created
    result <- tester.set.HD.batch(data, n.gene=n.gene)
                                  # suppressWarnWagner = TRUE, # if not suppressed, the slurm-out file size explodes.
                                  # LB.skip = F,LN.skip = F, MAST.skip = F,
                                  # KW.skip = F, Wg.skip = F, De2.skip = F, WRS.skip = F,
                                  # MGS.skip = F, ANCOM.skip = F, skip.cumulative = FALSE) # if there are batch effects skip WRS and MGS.
    result$nonzero.prop <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
    result$pval.cdf <- 
      result$pval %>% apply(1, function(x) {
        if (all(is.na(x))) rep(NA, length(cdf.cutoff)) else ecdf(x)(cdf.cutoff)
      }) %>% t
    attr(result$pval.cdf, "cutoff") = cdf.cutoff
    
    # Marking the cutoff rank declared as discovery by LEfSe.
    attr(result, "cutoff.LEfSe") = max(result$pval["LFE", ], na.rm = TRUE)
    
    ## Statistics needed for sensitivity and other metrics
    index.TP = 1:n.signal
    index.TN = (n.signal + 1):n.gene
    
    result$pval.cdf.TP <- 
      result$pval[, index.TP] %>% apply(1, function(x) {
        if (all(is.na(x))) rep(NA, length(cdf.cutoff)) else ecdf(x)(cdf.cutoff)
      }) %>% t
    attr(result$pval.cdf.TP, "cutoff") = cdf.cutoff
    
    result$pval.cdf.TN <- 
      result$pval[, index.TN] %>% apply(1, function(x) {
        if (all(is.na(x))) rep(NA, length(cdf.cutoff)) else ecdf(x)(cdf.cutoff)
      }) %>% t
    attr(result$pval.cdf.TN, "cutoff") = cdf.cutoff
    
    
    cdf.cutoff.q = c(0, 0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.2, 0.3)
    ## q-values for FDR
    result$qval = t(apply(result$pval, 1, function(x) p.adjust(x, method = "fdr")))
    result$qval.cdf.TP <- 
      result$qval[, index.TP] %>% apply(1, function(x) {
        if (all(is.na(x))) rep(NA, length(cdf.cutoff.q)) else ecdf(x)(cdf.cutoff.q)
      }) %>% t
    attr(result$qval.cdf.TP, "cutoff") = cdf.cutoff.q
    
    result$qval.cdf.TN <- 
      result$qval[, index.TN] %>% apply(1, function(x) {
        if (all(is.na(x))) rep(NA, length(cdf.cutoff.q)) else ecdf(x)(cdf.cutoff.q)
      }) %>% t
    attr(result$qval.cdf.TN, "cutoff") = cdf.cutoff.q
    
    
    ## Statistics needed for CATplot (concordance at top)
    result$ranks.TP <- 
      t(apply(result$pval[, index.TP], 1, function(x) ifelse(is.na(x), NA, order(x))))
    result$pval.metrics <- metrics(result$pval.cdf.TP, result$pval.cdf.TN, PN.rate = portion.signal, cutoff.LEF = cdf.cutoff[which.min(abs(cdf.cutoff - attr(result, "cutoff.LEfSe")))])
    result$qval.metrics <- metrics(result$qval.cdf.TP, result$qval.cdf.TN, PN.rate = portion.signal, cutoff.LEF = cdf.cutoff.q[which.min(abs(cdf.cutoff.q - attr(result, "cutoff.LEfSe")))])
    result.metrics <- list("pval" = result$pval.metrics, "qval" = result$qval.metrics)

   if (FALSE) # This section is ignored. Jan 19, 2022.
     { ### More MAST, DESeq2, MGS replicates (M = 10 in total)
          result.MGS <- result.DS2ZI <- result.DS2 <- result.MAST <- 
            result.ANCOM.sz <- result.ANCOM <- list(coef = list(), pval = list())
          
          result.MAST$coef[[1]] <- result$coef[c("MAST.nonz", "MAST.zero", "MAST.glob"), ]
          result.MAST$pval[[1]] <- result$pval[c("MAST.nonz", "MAST.zero", "MAST.glob"), ]
          result.MAST$nonzero.prop[[1]] <- result$nonzero.prop
          
          result.DS2$coef[[1]] <- result$coef["DS2", ]
          result.DS2$pval[[1]] <- result$pval["DS2", ]
          result.DS2$nonzero.prop[[1]] <- result$nonzero.prop
          
          result.DS2ZI$coef[[1]] <- result$coef["DS2ZI", ]
          result.DS2ZI$pval[[1]] <- result$pval["DS2ZI", ]
          result.DS2ZI$nonzero.prop[[1]] <- result$nonzero.prop
          
          
          result.MGS$coef[[1]] <- result$coef["MGS", ]
          result.MGS$pval[[1]] <- result$pval["MGS", ]
          result.MGS$nonzero.prop[[1]] <- result$nonzero.prop
          
          result.ANCOM.sz$coef[[1]] <- result$coef["ANCOM", ]
          result.ANCOM.sz$pval[[1]] <- result$pval["ANCOM", ]
          result.ANCOM.sz$nonzero.prop[[1]] <- result$nonzero.prop
          
          result.ANCOM$coef[[1]] <- result$coef["ANCOM", ]
          result.ANCOM$pval[[1]] <- result$pval["ANCOM", ]
          result.ANCOM$nonzero.prop[[1]] <- result$nonzero.prop
          
          # empty shells for the rest of the replicates
          for (s in 2:M) {
            result.MAST$coef[[s]] <- 
              result.MAST$pval[[s]] <- 
              matrix(NA, 3, n.gene, dimnames = list(c("MAST.nonz", "MAST.zero", "MAST.glob"), NULL))
            result.DS2$coef[[s]] <- result.DS2$pval[[s]] <- result.DS2ZI$coef[[s]] <- 
              result.DS2ZI$pval[[s]] <- result.MGS$coef[[s]] <- result.MGS$pval[[s]] <- 
              result.ANCOM.sz$coef[[s]] <- result.ANCOM.sz$pval[[s]] <- 
              result.ANCOM$coef[[s]] <- result.ANCOM$pval[[s]] <- rep(NA, n.gene)
          }
          message("More MAST, MGS, DESeq2, ANCOM replicates (s):\n")
          for (s in 2:M) {
            # # bookkeeping
            # if (file.exists(nm) & !ds.fatal) file.remove(nm)
            # nm = gsub("_s[0-9]*", sprintf("_s%s", s), nm)
            # write.table(" ", nm)
            
            cat("More replicate s = ", s, "\n")
            set.seed(replica*10^4 + s*10^5 + i*10^3 + j*10^2 + k)
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
            
            # MAST
            print("MAST")
            tmp.MAST <- try({MAST(data[, index.filtered.meta])})
            if (class(tmp.MAST)[1] == "try-error") tmp.MAST = matrix(NA, ncol = 2)
            result.MAST$coef[[s]][, index.filtered] <- tmp.MAST[[1]][1:3, ] #coef. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
            result.MAST$pval[[s]][, index.filtered] <- tmp.MAST[[2]][1:3, ] #pval. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
            result.MAST$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
  
            print("DESeq2 -- vanilla")
            # if (!s %in% book.s$s) {
            # DS2 = switch(DS2.version, vanilla = DS2.vanilla, zinb = DS2.zinb)
            DS2 = DS2.vanilla
            tmp.DS2 <- try({DS2(data[, index.filtered.meta])})
            if (class(tmp.DS2)[1] == "try-error") tmp.DS2 = matrix(NA, ncol = 2, dimnames = list(NULL, c("Estimate", "pval")))
            result.DS2$coef[[s]][index.filtered] <- tmp.DS2[, "Estimate"] #coef.
            result.DS2$pval[[s]][index.filtered] <- tmp.DS2[, "pval"] #pval.
            result.DS2$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
  
            print("DESeq2 -- zinbwave")
            DS2 = DS2.zinb
            tmp.DS2 <- try({DS2(data[, index.filtered.meta])})
            if (class(tmp.DS2)[1] == "try-error") tmp.DS2 = matrix(NA, ncol = 2, dimnames = list(NULL, c("Estimate", "pval")))
            result.DS2ZI$coef[[s]][index.filtered] <- tmp.DS2[, "Estimate"] #coef.
            result.DS2ZI$pval[[s]][index.filtered] <- tmp.DS2[, "pval"] #pval.
            result.DS2ZI$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
            #   ds.fatal = FALSE# otherwise it gives a fatal error, so this replicate is skipped.
            # } else {
            #   ds.fatal = TRUE# otherwise it gives a fatal error, so this replicate is skipped.
            # }
  
            print("MGS")
            tmp.MGS <- try({mgs(data[, index.filtered.meta])})
            if (class(tmp.MGS)[1] == "try-error") tmp.MGS = matrix(NA, ncol = 2)
            result.MGS$coef[[s]][index.filtered] <- tmp.MGS[, "Estimate"] #coef.
            result.MGS$pval[[s]][index.filtered] <- tmp.MGS[, "pval"] #pval.
            result.MGS$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
            
            print("ANCOM with the structural zero rule")
            tmp.ANC <- try({ANC(data[, index.filtered.meta], ignore.structural.zero = FALSE)})
            if (class(tmp.ANC)[1] == "try-error") tmp.ANC = matrix(NA, ncol = 2, dimnames = list(NULL, c("Estimate", "pval")))
            result.ANCOM.sz$coef[[s]][index.filtered] <- tmp.ANC[, "Estimate"] #coef.
            result.ANCOM.sz$pval[[s]][index.filtered] <- tmp.ANC[, "pval"] #pval.
            result.ANCOM.sz$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
            
            print("ANCOM")
            # tmp.ANC <- try({ANC(data[, index.filtered.meta])})
            # if (class(tmp.ANC)[1] == "try-error") tmp.ANC = matrix(NA, ncol = 2, dimnames = list(NULL, c("Estimate", "pval")))
            result.ANCOM$coef[[s]][index.filtered] <- tmp.ANC[, "Estimate"] #coef.
            result.ANCOM$pval[[s]][index.filtered] <- tmp.ANC[, "pval"] #pval.
            result.ANCOM$pval[[s]][is.infinite(result.ANCOM$coef[[s]])] <- NA
            result.ANCOM$nonzero.prop[[s]] <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0, na.rm = TRUE))
            
            gc()
          }
          
          # plug in back the MAST result in the main object.
          result$MAST  <- result.MAST
          result$DS2   <- result.DS2
          result$DS2ZI <- result.DS2ZI
          result$MGS   <- result.MGS
          result$ANCOM.sz <- result.ANCOM.sz
          result$ANCOM <- result.ANCOM
          
          
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
          
            result$pval.cdf[c("MAST.nonz", "MAST.zero", "MAST.glob"), ] <- 
              result$MAST$pval %>% Reduce(cbind, .) %>% apply(1, function(x) {
                if (all(is.na(x))) rep(NA, length(cdf.cutoff)) else ecdf(x)(cdf.cutoff)
              }) %>% t
            
            result$pval.cdf["DS2", ] <- result$DS2$pval %>% Reduce(c, .) %>% {if (all(is.na(.))) NA else ecdf(.)(cdf.cutoff)}
            result$pval.cdf["DS2ZI", ] <- result$DS2ZI$pval %>% Reduce(c, .) %>% {if (all(is.na(.))) NA else ecdf(.)(cdf.cutoff)}
            result$pval.cdf["MGS", ] <- result$MGS$pval %>% Reduce(c, .) %>% {if (all(is.na(.))) NA else ecdf(.)(cdf.cutoff)}
            result$pval.cdf["ANCOM.sz", ] <- result$ANCOM.sz$pval %>% Reduce(c, .) %>% {if (all(is.na(.))) NA else ecdf(.)(cdf.cutoff)}
            result$pval.cdf["ANCOM", ] <- result$ANCOM$pval %>% Reduce(c, .) %>% {if (all(is.na(.))) NA else ecdf(.)(cdf.cutoff)}
      }    
          
    # 3. Statistics
    # 3.1 statistics for all method (including first replicate of MAST, which is redundant)
    
    ## This part is not needed as it is done in the filtering step.
    # index.regular <- result$nonzero.prop >= prev.filter
    # result$pval[, !index.regular] <- NA    # removing irregular genes.
    
    
    stat.power = result$pval %>% 
      apply(1, function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))})
    stat.irregular = result$pval %>% 
      apply(1, function(x) {ifelse(sum(is.na(x))/length(x) > 0.9, NA, mean(ifelse(is.na(x), 1, x) <= sig, na.rm=TRUE))})
    stat.na.prop = result$pval %>% apply(1, function(x) {mean(is.na(x))})
    
    # 3.2 statistics for MAST replicates
      ## This part is not needed as it is done in the filtering step.
      # for (m in 1:M) {
      #   index.regular.m <- result$MAST$nonzero.prop[[m]] >= prev.filter
      #   result$MAST$pval[[m]][, !index.regular.m] <- NA
      # }
    if (FALSE) # This section is ignored. Jan 19, 2022.
    {  
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
      ## This part is not needed as it is done in the filtering step.
      # for (m in 1:M) {
      #   index.regular.m <- result$MGS$nonzero.prop[[m]] >= prev.filter
      #   result$MGS$pval[[m]][!index.regular.m] <- NA
      # }
      
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
      
      
      stat.power["DS2"] <-
        result$DS2$pval %>% 
        sapply( function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))}) %>%
        mean(na.rm = TRUE)  # vector of three
      
      stat.irregular["DS2"] <- #stat based on NA p-val = 1.
        result$DS2$pval %>% 
        sapply(function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)}) %>%
        mean(na.rm = TRUE)
      
      stat.na.prop["DS2"] <-
        result$DS2$pval %>% 
        sapply( function(x)  {mean(is.na(x))}) %>%
        mean(na.rm = TRUE)
      gc()
      
      
      stat.power["ANCOM.sz"] <-
        result$ANCOM.sz$pval %>% 
        sapply( function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))}) %>%
        mean(na.rm = TRUE)  # vector of three
      
      stat.irregular["ANCOM.sz"] <- #stat based on NA p-val = 1.
        result$ANCOM.sz$pval %>% 
        sapply(function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)}) %>%
        mean(na.rm = TRUE)
      
      stat.na.prop["ANCOM.sz"] <-
        result$ANCOM.sz$pval %>% 
        sapply( function(x)  {mean(is.na(x))}) %>%
        mean(na.rm = TRUE)
      gc()
      
      
      stat.power["ANCOM"] <-
        result$ANCOM$pval %>% 
        sapply( function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, mean(x<=sig, na.rm=TRUE))}) %>%
        mean(na.rm = TRUE)  # vector of three
      
      stat.irregular["ANCOM"] <- #stat based on NA p-val = 1.
        result$ANCOM$pval %>% 
        sapply(function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)}) %>%
        mean(na.rm = TRUE)
      
      stat.na.prop["ANCOM"] <-
        result$ANCOM$pval %>% 
        sapply( function(x)  {mean(is.na(x))}) %>%
        mean(na.rm = TRUE)
      gc()
    }
    
    
# Save the result
    stat.comb <- rbind(stat.power, stat.irregular, stat.na.prop)
    
    saveRDS(list(stat = stat.comb, cdf = result$pval.cdf, pval.cdf.TP = result$pval.cdf.TP, pval.cdf.TN = result$pval.cdf.TN, metrics = result.metrics, setting = setting.summary), save_file.stat)
    if (!save.stat.only) saveRDS(result, save_file.raw)
    
    # # bookkeeping
    # if (file.exists(nm) & !ds.fatal) file.remove(nm)
    # # nm = gsub("tmp_", "tmp_done_", nm)
    # # nm = gsub("_s[0-9]*", "", nm)
    # # write.table(" ", nm)
    
    tt(2)
}

