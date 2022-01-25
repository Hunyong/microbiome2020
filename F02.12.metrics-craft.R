### Functions for evaluating craft data

tab.metrics = function(zoe, type, n.signal, n.gene, reps = 1:10, suppressMsg = FALSE, cutoff = 0.05) {
  require(dplyr)
  require(abind)
  
  individual.tables = 
    lapply(reps, function(l) {
      fn = sprintf("output/stat-craft%d-%s-nSig%d-nGene%d-j2-rep%d.rds", zoe, type, n.signal, n.gene, l)
      if (!file.exists(fn)) {
        a = NULL
        if (!suppressMsg) cat("missing", fn, "\n")
        NULL
      } else {
        a = readRDS(fn)
        metrics(a$cdf.TP, a$cdf.TN, PN.rate = a$setting$n.signal / a$setting$n.gene, cutoff = cutoff)
      }
    })
  if (all(sapply(individual.tables, is.null))) return(NULL)
  Reduce(function(...) abind(..., along = 3), individual.tables) %>%    # a list into a 3d array
    apply(1:2, mean, na.rm = TRUE) %>% round(3)                         # average over replicates
}
# tab.metrics(zoe = 3, type = "gene", n.signal = 100, n.gene = 10000) 
# tab.metrics(zoe = 3, type = "gene", n.signal = 1000, n.gene = 10000) 
