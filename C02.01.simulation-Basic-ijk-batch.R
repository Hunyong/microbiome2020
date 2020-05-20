cat("Microbiome - C02.01.simulation-Basic-ijk-batch.R\n")

args = commandArgs(trailingOnly=TRUE)  # passed from script
cat("The Command Arg is: ", str(args))
i = as.numeric(args[1])  # 1..10    delta effect
j = as.numeric(args[2])  # 1..5     kappa effect
k = as.numeric(args[3])  # 1..34    baseline scenario
model = as.numeric(args[4])  # 1..3 generative model
perturb = as.numeric(args[5]) # 5, 3, 0
n = as.numeric(args[6])  # 80 800  sample size

if(F){
  i = 1  # 1..10    delta effect
  j = 1  # 1..5     kappa effect
  k = 1  # 1..34    baseline scenario
  model = 3  # 1..3 generative model
  perturb = 5 # 5, 3, 0
  n = 80  # 80 800  sample size
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

cat("i: ", i,", j: ",j,", k: ",k,", model: ", model, ", perturb: ", perturb, "\n")


## To save the result
# Check and create the folder
save_path = paste0("output/")
save_file = paste0("output/result-n", n, "-pert", perturb, "-", model, "-", i, ".", j, ".", k, ".rds")

if (!dir.exists(save_path)) {message("No output folder detected. Creating one."); dir.create(save_path)}
if (file.exists(save_file)) {stop("done already.")}

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

n.sim; n.sample; 
n.genes = 1e+4
print(test.dim <- method.stat %>% length)

#


tt(1)
set.seed(i*10^3 + j*10^2 + k)
# 1. parameter

## param.set = param (i, j, k) # list of H1, D1, H2, D2 (status-batch)
param.set = param (i, j, k, baseParam = parameter, delta.table = delta, kappa.table = kappa)
dat.args = list(n.sample = n.sample, n.genes=n.genes, scenario.delta = i, scenario.kappa = j, scenario.base = k,
                baseParam = parameter,
                delta.table = delta, 
                kappa.table = kappa,
                model = model, 
                delta.perturb.prob = perturb)
# 2. data
data = do.call(r.sim, dat.args)
data %<>% dplyr::filter(sampleSum > 0)
cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")

#if (any(class(try(readRDS(paste0("output/R0201sim181201/result.", i, ".", j, ".", k,".rds")))) %in% "try-error")) {

# do the tests on the ramdon ZINB distribution we created
result <- tester.set.HD.batch(data, n.sim=n.sim, suppressWarnWagner = TRUE, sig = 0.05,
                              LB.skip = F,LN.skip = F, MAST.skip = F,
                              KW.skip = F, Wg.skip = F, De2.skip = F) # if not suppressed, slurm-out file size explodes.
result$nonzero.prop <- apply(data[, 1:n.genes], 2, function(s) mean(s > 0))

### More MAST replicates (10 in total)
result.MAST <- list(coef = list(), pval = list())
result.MAST$coef[[1]] <- result$coef[c("MAST.nonz", "MAST.zero", "MAST.glob"), ]
result.MAST$pval[[1]] <- result$pval[c("MAST.nonz", "MAST.zero", "MAST.glob"), ]
result.MAST$nonzero.prop[[1]] <- result$nonzero.prop

message("More MAST replicates (s):\n")
for (s in 2:10) {
  cat("MAST replicate s = ", s, "\n")
  set.seed(s*10^5 + i*10^3 + j*10^2 + k)
  data = do.call(r.sim, dat.args)
  data %<>% dplyr::filter(sampleSum > 0)
  cat("sample size is ", dim(data)[1], "out of ", sum(n.sample), ".\n")
  tmp.MAST <- MAST(data, sig = 0.05)
  result.MAST$coef[[s]] <- tmp.MAST[[1]][1:3, ] #coef. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
  result.MAST$pval[[s]] <- tmp.MAST[[2]][1:3, ] #pval. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
  result.MAST$nonzero.prop[[s]] <- apply(data[, 1:n.genes], 2, function(s) mean(s > 0))
}

# plug in back the MAST result in the main object.
result$MAST <- result.MAST

# Save the result
saveRDS(result, save_file)
tt(2)
# R0201sim180718 #n.sim=10000,  R0201sim180711 #n.sim=1000
#}

