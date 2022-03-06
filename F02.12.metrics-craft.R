### Functions for evaluating craft data
tab.metrics = function(zoe, type, n.signal, n.gene, reps = 1:10, suppressMsg = FALSE, cutoff = 0.05) {
  require(dplyr)
  require(abind)
  
  individual.tables = 
    lapply(reps, function(l) {
      fn = sprintf("output/stat-craft%d-%s-nSig%d-nGene%s-j2-rep%d.rds", zoe, type, n.signal, n.gene, l)
      if (!file.exists(fn)) {
        a = NULL
        if (!suppressMsg) cat("missing", fn, "\n")
        NULL
      } else {
        a = readRDS(fn)
        # otherwise the cutoff will stay the same (This new cutoff is only valid within the current l.)
        PN.rate = a$setting$n.signal / a$setting$n.gene
        res = metrics(a$cdf.TP, a$cdf.TN, PN.rate = PN.rate, cutoff = cutoff)
        
        # q-values for the FDR statistics.
        res2 = metrics(a$cdf.TP.q, a$cdf.TN.q, PN.rate = PN.rate, cutoff = cutoff)
        res[, "FDR"] = res2[, "FDR"]
        
        # For LEfSe
        res["LFE", "sensitivity"] = sens.LFE = attr(a$cdf.TP, "cutoff")[min(which(a$cdf.TP["LFE", ] == 1))]
        res["LFE", "type1error"]  = fpr.LFE  = attr(a$cdf.TN, "cutoff")[min(which(a$cdf.TN["LFE", ] == 1))]
        res["LFE", "accuracy"]    = sens.LFE * PN.rate + (1 - fpr.LFE) * (1 - PN.rate)
        res["LFE", "FDR"]         = res["LFE", "AUC"]  = NA
        
        return(res)
      }
    })
  if (all(sapply(individual.tables, is.null))) return(NULL)
  Reduce(function(...) abind(..., along = 3), individual.tables) %>%    # a list into a 3d array
    apply(1:2, mean, na.rm = TRUE) %>% round(3)                         # average over replicates
}
# tab.metrics(zoe = 3, type = "gene", n.signal = 100, n.gene = 10000) 
# tab.metrics(zoe = 3, type = "gene", n.signal = 1000, n.gene = 10000) 

