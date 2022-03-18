### functions for sensitivity, FDR, accuracy, etc.

#' @example
#' a = readRDS("output/stat-craft1-gene-nSig100-nGene10000-j2-rep1.rds")
#' metrics(a$cdf.TP, a$cdf.TN, PN.rate = a$setting$n.signal / a$setting$n.gene, cutoff = 0.05) 

process.zero <- function(metrics, value){
  for (i in c(1:length(1:length(metrics)))){
    metrics[i] <- max(metrics[i], value)
  }
  return(metrics)
}

metrics = function(cdf.TP, cdf.TN, PN.rate, cutoff.LEF, cutoff = 0.05) {
  sens = process.zero(cdf.TP[, attr(cdf.TP, "cutoff") == cutoff], 1e-10)    # Power
  fpr  = cdf.TN[, attr(cdf.TN, "cutoff") == cutoff]    # Type 1 error
  # cdf  = cdf.TP * PN.rate + cdf.TN * (1 - PN.rate)
  fdr  = fpr * (1 - PN.rate) / (sens * PN.rate + fpr * (1 - PN.rate))
  acc  = {sens * PN.rate + (1 - fpr) * (1 - PN.rate)}
  AUC  = auc(cdf.TP, cdf.TN)
  
  output = data.frame(sensitivity = sens, type1error = fpr, FDR = fdr, accuracy = acc, AUC = AUC)

  sens.LEF = max(cdf.TP[rownames(cdf.TP)=="LFE", attr(cdf.TP, "cutoff") == cutoff.LEF], 1e-10)
  fpr.LEF  = cdf.TN[rownames(cdf.TP)=="LFE", attr(cdf.TN, "cutoff") == cutoff.LEF]    # Type 1 error
  # cdf  = cdf.TP * PN.rate + cdf.TN * (1 - PN.rate)
  fdr.LEF  = fpr.LEF * (1 - PN.rate) / (sens.LEF * PN.rate + fpr.LEF * (1 - PN.rate))
  acc.LEF  = {sens.LEF * PN.rate + (1 - fpr.LEF) * (1 - PN.rate)}

  output["LFE", "sensitivity"] <- sens.LEF
  output["LFE", "type1error"] <- fpr.LEF
  output["LFE", "FDR"] <- fdr.LEF
  output["LFE", "accuracy"] <- acc.LEF
  attr(output, "legend") = "power = sensitivity, type 1 error = fpr"
  
  return(output)
}
auc = function(cdf.TP, cdf.TN) { 
  n.methods = dim(cdf.TP)[1]
  sapply(1:n.methods, function(l) auc.l(cdf.TP[l, ], cdf.TN[l, ]))
}
auc.l = function(cdf.TP.l, cdf.TN.l) {
  if (any(is.na(c(cdf.TP.l, cdf.TN.l)))) return(NA)
  n = length(cdf.TP.l)
  F1 = cdf.TP.l[-1] * 0.5 + cdf.TP.l[-n] * 0.5
  dF0 = cdf.TN.l[-1] - cdf.TN.l[-n]
  sum(F1 * dF0)
}