craft = function(otu.matrix, meanEffect, batchVec, n.signal, n.gene,
                 replace = F, deltaMatrix, cut.delta.mu = NULL, cut.delta.pi = NULL,
                 prev.filter = 0.1) {
  
#### low prevalence genes are already filtered out in otu.matrix.
  p = dim(otu.matrix)[1]
  n = dim(otu.matrix)[2]
  effects = matrix(0, nrow = n.signal, ncol = 3, dimnames = list(1:n.signal, c("delta.mu", "delta.pi", "delta.theta")))
  
  ## 1. gene sampling (n.genes)
  if (1) {
    # Prevalance filtering is already done in otu.matrix!
    if (n.gene > p) {
      warning("The number of genes in the data being sampled from is less than n.gene. Sampling with replacement.")
      replace = TRUE
    }
    index.p = sample(1:p, n.gene, replace = replace)
    otu.matrix = otu.matrix[index.p, ]
    index.sig = 1:n.signal   # The first n.signal genes contain signal.
    
  } else {
    # When the otu.matrix without prevalance filtering was provided.
    prev.avg = 0
    while (prev.avg < 0.1) { # Should achieve average prev. rate >= 10%.
      index.p = sample(1:p, n.gene, replace = replace)
      otu.matrix.tmp = otu.matrix[index.p, ]
      prev.otu = apply(otu.matrix.tmp, 1, function(x) mean(x > 0, na.rm = TRUE))
      prev.avg = mean(prev.otu >= prev.filter)
    }
    otu.matrix = otu.matrix.tmp
    index.sig = 1:n.signal
  }
  otu.crafted = otu.matrix # The skeleton
  # significant genes are sampled from the prevalence rate >= filter.
  
  ## 2. permutation stratified by batch
  lvls = sort(unique(batchVec))
  batch.index = lapply(lvls, function(x) which(batchVec == x))
  group.index = lapply(batch.index, function(x) {d = rbinom(length(x), 1, 0.5); list(D = x[d == 1], H = x[d == 0])})
  group.index = list(D1 = group.index[[1]]$D, H1 = group.index[[1]]$H,
                     D2 = group.index[[2]]$D, H2 = group.index[[2]]$H)
  diseaseVec = rep("H", n)
  diseaseVec[unlist(group.index[c(1, 3)])] = "D"
  
  ## 3. filtering out the significant disease effects from the reference delta table.
  # 3.1 mu effects
  if (!is.null(cut.delta.mu)) {
    delta.mu.index = abs(deltaMatrix[, "delta_mu"]) >= cut.delta.mu
  } else {
    delta.mu.index = rep(TRUE, dim(deltaMatrix)[1])
  }
  # 3.2 pi effects
  if (!is.null(cut.delta.pi)) {
    delta.pi.index = abs(deltaMatrix[, "delta_pi"]) >= cut.delta.pi
  } else {
    delta.pi.index = rep(TRUE, dim(deltaMatrix)[1])
  }
  delta.index = delta.mu.index & delta.pi.index
  # 3.3 Leaving only the significant effects
  deltaMatrix = deltaMatrix[delta.index, ]
  deltaSamp = deltaMatrix[sample(1:dim(deltaMatrix)[1], size = n.signal, replace = TRUE),]
  
  ## 4.0 basic dimensions
  lens = sapply(group.index, length)
  
  ## 4. Adding disease effects
  for (i in index.sig) {
    
    # 4.1 delta_pi
    d.pi = deltaSamp[i, "delta_pi"]
    effects[i, "delta.pi"] = d.pi
    group.index.zero = lapply(group.index, function(x) otu.matrix[i, x] == 0)
    zero.prop = sapply(group.index.zero, mean) # D/H x batchs
    zero.prop.new = plogis(qlogis(zero.prop) + c(d.pi/2, -d.pi/2))
    zero.surplus = zero.prop * lens - sapply(1:4, function(group) rbinom(1, lens[group], zero.prop.new[group]))
    for (group in 1:4) {
      if (zero.surplus[group] > 0) { # replace zeros with some nonzeros from the vector
        subj.replaced = sample(group.index[[group]][group.index.zero[[group]]], size = zero.surplus[group], replace = FALSE)
        otu.crafted[i, subj.replaced] = otu.matrix[i, sample(group.index[[group]][!group.index.zero[[group]]], size = length(subj.replaced), replace = TRUE)]
        # group.index.zero[[group]][subj.replaced] = !group.index.zero[[group]][subj.replaced] # updating for the next steps (delta_mu)
      } else if (zero.surplus[group] < 0) {
        subj.replaced = sample(group.index[[group]][!group.index.zero[[group]]], size = -zero.surplus[group], replace = FALSE)
        otu.crafted[i, subj.replaced] = 0
        # group.index.zero[[group]][subj.replaced] = !group.index.zero[[group]][subj.replaced] # updating for the next steps (delta_mu)
      }
    }
    
    # 4.2 delta_mu
    d.mu = deltaSamp[i, "delta_mu"]
    effects[i, "delta.mu"] = d.mu
    otu.crafted[i, group.index[c(1,3)] %>% unlist] = otu.matrix[i,  group.index[c(1,3)] %>% unlist]  * exp(d.mu/2)
    otu.crafted[i, group.index[c(2,4)] %>% unlist] = otu.matrix[i,  group.index[c(2,4)] %>% unlist]  * exp(-d.mu/2)
  }
  
  ## 5. flip and add meta data
  dat = data.frame(t(otu.crafted), phenotype = factor(diseaseVec), batch = factor(batchVec))
  names(dat)[1:n.gene] = paste0("y.", 1:n.gene)
  dat$sampleSum = dplyr::select(dat, -phenotype, -batch) %>% apply(1, sum)
  attr(dat, "gene.names") = rownames(otu.crafted)
  attr(dat, "effects") = effects
  
  return(dat)
}