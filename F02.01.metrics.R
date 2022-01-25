### functions for sensitivity, FDR, accuracy, etc.

#' @example
#' metrics(a$cdf.TP, a$cdf.TN, PN.rate = a$setting$n.signal / a$setting$n.gene, cutoff = 0.05) 
metrics = function(cdf.TP, cdf.TN, PN.rate, cutoff = 0.05) {
  sens = cdf.TP[, attr(cdf.TP, "cutoff") == cutoff]
  fpr  = cdf.TN[, attr(cdf.TN, "cutoff") == cutoff]    # Type 1 error
  fdr  = fpr * (1 - PN.rate) / (sens * PN.rate + fpr * (1 - PN.rate))
  acc  = {sens * PN.rate + (1 - fdr) * (1 - PN.rate)}
  AUC  = auc(cdf.TP, cdf.TN)
  
  output = data.frame(sensitivity = sens, type1error = fpr, FDR = fdr, accuracy = acc, AUC = AUC)
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