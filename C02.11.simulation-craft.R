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
source("F02.11.craftMaker.R")
library(MAST)
library(coin)
library(metagenomeSeq)

# required parameters from...
source("C01.02.simulation.setup.R")
# sessionInfo()


args = commandArgs(trailingOnly=TRUE)  # passed from script
if (length(args) == 0) {
  warning("commandArgs() was not provided. Set as the default value.")
  args = c(zoe = 1, j = 2, type = "gene", sim = 1, 
           n.signal = 100, n.gene = 1e+5, save.stat.only = 1)
}
cat("The Command Arg is: \n"); print(args)
# i = 1                       # delta effect
# k = 0                       # baseline scenario (Craft)
# model = 0                   # craft
# perturb = 0                 # Permuting is naturally considered.
# n = NA  # 80 800  sample size
zoe = as.numeric(args[1]) # 1..3  generative model (1 = ZOE 1, 2 = 2, 3 = IBD)
j = as.numeric(args[2])     # 0..1 kappa effect (No effect = permute, preserve the batch)
type = as.character(args[3])# gene genebact bact
sim = as.numeric(args[4]) # 1..100 simulation replicate id
n.signal = as.numeric(args[5]) # 100 300 1000
n.gene = as.numeric(args[6])  # Effective number of genes = 10K out of 100K (around 90% of them have prevalence rate < 10%. 50% of them are all zeros).
save.stat.only = as.logical(as.numeric(args[7])) # 1, 0

if (is.na(save.stat.only)) save.stat.only = TRUE
if (is.na(n.gene)) n.gene = 1e+5

prev.filter = 0.1
model = "craft"

# cat("ZOE data: ", zoe, ", j: ", j, ", sim: ", sim, ", model: ", model, "\n")
# cat("n.gene: ",n.gene, ", n.signal: " , n.signal, ", stat.stat.only : ", save.stat.only,"\n")

## To save the result
# Check and create the folder
save_path = paste0("output/")
save_file.raw = paste0("output/raw-craft", zoe, "-", type, "-nSig", n.signal, "-nGene", n.gene, "-j", j, "-rep", sim, ".rds")
save_file.stat = paste0("output/stat-craft", zoe, "-", type, "-nSig", n.signal, "-nGene", n.gene, "-j", j, "-rep", sim, ".rds")

if (!dir.exists(save_path)) {message("No output folder detected. Creating one."); dir.create(save_path)}
if (file.exists(save_file.stat)) stop("Already done")


type.full = switch(type, gene = "geneRPK.marginal", 
                   genebact = "geneRPK.joint", bact = "bactRPK.marginal")
zoe.nm = if (zoe %in% 1:2) paste0("_zoe", zoe) else "_NEWDATA"
zoe.nm2 = if (zoe %in% 1:2) paste0("ZOE", zoe) else "NEWDATA"

# Read data
if (zoe %in% 1:2) {
  fn       = sprintf("../Data-processed/data.%s.DRNA.%s.rds", type.full, zoe.nm2)
  fn.est   = sprintf("output/para_selection_est_zoe%s_%s_ziln_tpm5_samp%s.rds", zoe, type, "full")
  fn.delta = sprintf("output/para_delta_est%s_%s_%s_%s_samp%s.rds", zoe.nm, type, "ziln", "tpm5", "full")
  dat.raw     = readRDS(fn)
  # dat.est = readRDS(fn.est)
  cond.est.delta = readRDS(fn.delta)
  
  excluded.subject = dat.raw$meta$id %in% c(352, 420, 10083, 12623, 11259, 11790)
  dat.raw$otu = dat.raw$otu[, !excluded.subject, ]
  dat.raw$meta = dat.raw$meta[!excluded.subject, ]
  batchGrp = switch(zoe, "1" = "170628", "2" = "180530")
  
  DataMeta =
    dat.raw$meta %>% 
    mutate(group = paste0(ifelse(cariesfree == 1, "H", "D"), ifelse(batch.RNA == batchGrp, 1, 2)),
           group_batch =  ifelse(batch.RNA == batchGrp, 1, 2),
           group_disease = ifelse(cariesfree == 1, "H", "D"))
} else {
  fn <- sprintf("Nature2019data/data.%s.DRNA.IBD.rds", type.full)
  fn.est = sprintf("output/para_selection_est_NEWDATA_%s_ziln_tpm5_samp%s.rds", type, "full")
  fn.delta = sprintf("output/para_delta_est%s_%s_%s_%s_samp%s.rds", zoe.nm, type, "ziln", "tpm5", "full")
  dat.raw <- readRDS(fn)
  # dat.est <- readRDS(fn.est)
  cond.est.delta = readRDS(fn.delta)
  
  excluded.subject = dat.raw$meta$External.ID %in% "MSM9VZMA"
  dat.raw$otu = dat.raw$otu[, !excluded.subject]
  dat.raw$meta = dat.raw$meta[!excluded.subject, ]
  DataMeta = 
    dat.raw$meta %>% 
    mutate(id = External.ID,
           group_batch = ifelse(site_name %in% c("Cedars-Sinai", "MGH"), 1, 2),
           group_disease = ifelse(diagnosis %in% c("UC", "CD"), "D", "H"),
           group = paste0(group_disease, group_batch))
}

print(test.dim <- method.stat %>% length)
cdf.cutoff = c((0:100)/500, 0.2 + 0.8 * (1:80)/80) # higher resolution for < 0.2.

tt(1)
seed.no = sim * 10^2 + j*10^1
set.seed(seed.no)


### 0.2 Data
survivors = dat.raw$otu[,, 2] %>% apply(1, function(x) mean(x > 0, na.rm = TRUE) > prev.filter)
# for ZOE 1, sum(survivors) = 98879
data =
  craft (otu.matrix = dat.raw$otu[survivors,, 2], 
         meanEffect, batchVec = DataMeta$group_batch, 
         n.signal = n.signal, n.gene = n.gene, replace = F, 
         deltaMatrix = cond.est.delta, cut.delta.mu = 2, cut.delta.pi = NULL,
         prev.filter = prev.filter) 
# cond.est.delta$delta_mu %>% abs %>% quantile(0.95) ## 2.05
# cond.est.delta$delta_pi %>% abs %>% quantile(0.95) ## 1.85

setting.summary <- 
  data.frame(zoe = zoe,
             j = j,
             type = type,
             model = model,
             sim = sim,
             seed.no = seed.no,
             n.signal = n.signal,       # number of genes with signal
             n.gene = n.gene,           # number of genes included in craft data
             n.total = sum(survivors),  # number of genes being sampled from
             prev.filter = prev.filter)

attr(data, "setting") = setting.summary

# filtering
nonzero.prop <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))
filtr = nonzero.prop >= prev.filter
if (sum(!filtr)) data[, which(!filtr)] = NA
cat("Additional filtering after crafting: ", sum(!filtr), ".\n")

# cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
cat("Remaining genes after screening: ", sum(filtr), "out of ", length(filtr), ".\n")

# do the tests on the ramdon ZINB distribution we created
result <- tester.set.HD.batch(data, n.gene=n.gene, SONGBIRD.skip = TRUE) 
result$nonzero.prop <- apply(data[, 1:n.gene], 2, function(s) mean(s > 0))

# Marking the cutoff rank declared as discovery by LEfSe.
attr(result, "cutoff.LEfSe") = max(result$pval["LFE", ], na.rm = TRUE)

result$pval.cdf <- 
  result$pval %>% apply(1, function(x) {
    if (all(is.na(x))) rep(NA, length(cdf.cutoff)) else ecdf(x)(cdf.cutoff)
  }) %>% t
attr(result$pval.cdf, "cutoff") = cdf.cutoff


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
gc()

# Save the result
stat.comb <- rbind(stat.power, stat.irregular, stat.na.prop)

saveRDS(list(stat = stat.comb, cdf = result$pval.cdf, setting = setting.summary), save_file.stat)
if (!save.stat.only) saveRDS(result, save_file.raw)

# # bookkeeping
# if (file.exists(nm) & !ds.fatal) file.remove(nm)
# # nm = gsub("tmp_", "tmp_done_", nm)
# # nm = gsub("_s[0-9]*", "", nm)
# # write.table(" ", nm)

tt(2)

