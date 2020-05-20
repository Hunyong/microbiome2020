#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(BiocManager)
n.test = 13 #"LB.nonz", "LB.zero", "LB.glob", "LN", "MAST.nonz", "MAST.zero", "MAST.glob", "KW", "Wg.nonz", "Wg.zero", "Wg.glob", "Reserved"
# if(!require(betareg)){
#   install.packages("betareg")
# }


### 0. testing a set of methods at once
# testing all methods using counts (or RPK) of a single gene (old: single gene. cannot do MAST).
# HD: binary phenotype (healthy-diseased)
tester.set.HD.batch <- function(data, n.sim = 10000, sig = 0.05, skeleton = FALSE, 
                                suppressWarnWagner = FALSE, LB.skip = FALSE, 
                                LN.skip = FALSE, MAST.skip = FALSE,
                                KW.skip = FALSE, Wg.skip = FALSE,
                                De2.skip = FALSE,
                                fin.correction = FALSE, skip.small.n = FALSE) {
  # description
  # data should have y and sampleSum    all n.sample x (n.sim(gene) + 3 (phenotype + batch + sampleSum))
  #          outcome (phenotype), nuisance (batch)
  # skeleton: returning only skeleton (for simulation structure)
  require(magrittr)
  
  # 0.0 skeleton #empty matrix
  test.names <- c("LB.nonz", "LB.zero", "LB.glob", "LN", "MAST.nonz", "MAST.zero", "MAST.glob", "KW", "Wg.nonz", "Wg.zero", "Wg.glob", "DESeq2","(Reserved)")
  mat.tmp <- matrix(NA, n.test, n.sim, dimnames = list(test.names, NULL)) # n.test=12, n.sim=1000
  result <- list(coef = mat.tmp, pval = mat.tmp, coef.fin = mat.tmp, pval.fin = mat.tmp)
  
  if (skeleton) {return(result)}
  
  # 0.1 data
  data2 <- data
  gsub("y\\.","",names(data)) %>% as.numeric %>% na.omit %>% as.numeric -> genes
  if (n.sim > genes[length(genes)]) stop(paste0("Only ", length(genes), " genes provided, while trying to do ", n.sim, " simulations."))
  genes = genes[1:min(n.sim,length(genes))]
  
  data = data.frame(data[,1:length(genes)], phenotype = data$phenotype, batch = data$batch)
  ## print(head(data))    (No longer need to check data)
  if (!"sampleSum" %in% names(data2)) {
    data$sampleSum <- data %>% dplyr::select(-phenotype, - batch) %>% apply(1, sum)
  } else {
    data$sampleSum <- data2$sampleSum
  }
  
  
  cat("1-3. Logistic Beta\n")
  #1-3. LB
  if (!LB.skip) {
    for (l in genes) {
      cat (" l = ",l," ")
      # if (l %% 200 == 0) {print("l = ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y)==0) {
        tmp <- data.frame(coef = rep(NA,3), pval = NA)
      } else {
        ## print(data[,l]) (For debug only)
        tmp <- LB(data.l)  #logistic beta
      }
      result[[1]][1:3, l] <- tmp[1:3,1] #coef. 1:3 corresponds to "LB.nonz", "LB.zero", "LB.glob"
      result[[2]][1:3, l] <- tmp[1:3,2] #pval. 1:3 corresponds to "LB.nonz", "LB.zero", "LB.glob"
    }
  } else {cat("LB is skipped")}
  
  
  #4. LN
  
  cat("\n4. Log normal\n l = ")
  if (!LN.skip) {
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y)==0) {
        tmp <- data.frame(coef = NA, pval = NA)
      } else {
        tmp <- LN(data.l, sig = sig)  #log normal
      }
      result[[1]][4, l] <- tmp[1,1] #coef.
      result[[2]][4, l] <- tmp[1,2] #pval.
    }
  } else {cat("LN is skipped")}
  
  #5-7. MAST
  cat("\n5-7. MAST\n")
  if(!MAST.skip){
    # tmp <- MAST(data, sig = sig)  #MAST
    # tmp <- data.frame(coef = rep(NA,3), pval = NA)    #MAST maybe not applicable
    tmp <- MAST(data, sig = sig)
    result[[1]][5:7, ] <- tmp[[1]][1:3,] #coef. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
    result[[2]][5:7, ] <- tmp[[2]][1:3,] #pval. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
  } else {cat("MAST is skipped")}
  
  #8. KW
  cat("\n8. Kruskal Wallis\n l = ")
  if(!KW.skip){
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y)==0) {
        tmp <- data.frame(coef = NA, pval = NA)
      } else {
        tmp <- KW(data.l, sig = sig)  #KW
      }
      result[[1]][8, l] <- tmp[1,1] #coef.
      result[[2]][8, l] <- tmp[1,2] #pval.
    }
  } else {cat("KW is skipped")}
  
  #9-11. Wagner
  cat("\n9-11. Wagner (2-part)\n l = ")
  if(!Wg.skip){
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y)==0) {
        tmp <- data.frame(coef = rep(NA,3), pval = NA)
      } else {
        tmp <- Wagner(data.l, sig = sig, zeroModel = "logistic", suppressWarning = suppressWarnWagner)
      }
      result[[1]][9:11, l] <- tmp[1:3,1] #coef.
      result[[2]][9:11, l] <- tmp[1:3,2] #pval.
    }
  } else {cat("Wg is skipped")}
  
  #12. De2
  cat("\n12. DESeq2\n")
  if(!De2.skip){
    tmp <- DS2(data, sig = sig)
    
    result[[1]][12, ] <- tmp[,1] #coef.
    result[[2]][12, ] <- tmp[,2] #pval.
    
  }
  
  
  #13. (reserved)
  cat("13. Reserved\n")
  # tmp <- data.frame(coef = NA, pval = NA)    #reserved for possible addition
  # result[[1]][13,1] <- tmp[1,1]
  # result[[2]][13,1] <- tmp[1,2]
  
  return(result)
}


if (FALSE) {# examples
  data = rZINB.sim(n.sample=rep(3,4),n.genes=30, 1,1,1)
  data %>% dplyr::select(-phenotype, - batch) %>% apply(1, sum) -> data$sampleSum
  data.1 = data.frame(y=data[,1], data[,c("phenotype", "batch", "sampleSum")])
  a = tester.set.HD.batch(data, n.sim=30)
}


### 1~3. LB
LB.old <- function (data, sig = 0.05, test = FALSE) {
  require(gamlss)
  data$y.prop = data$y / data$sampleSum  #sampleSum[j] = sum_g y_g,j (g:gene, j:sample)
  # print(head(data$y.prop))  
  # 1. two-part models
  bereg = try(gamlss(y.prop ~ phenotype + batch, nu.formula = ~ phenotype + batch,
                     family = BEZI(sigma.link = "log"), data = data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE)))
  
  if (any(class(bereg) %in% "try-error")) {
    return(matrix(NA, 3, 2, 
                  dimnames = list(c("LB.nonz", "LB.zero", "LB.glob"), c("Estimate", "pval"))))}
  
  out1 = summary(bereg)[c(2,6), c(1,4)]
  
  # Alternative code should be updated (b/c when vcov fails, qr should be used.)
  # bereg.mat = c(try(suppressWarnings(vcov(bereg, type = "all", 
  #                   robust = F, hessian.fun = "R")), silent = TRUE), 
  #               df.res = bereg$df.res)
  # coef <- bereg.mat$coef
  # pvalue <- 2 * pt(-abs(coef/bereg.mat$se), bereg.mat$df.res)
  # pheno.index = grep("phenotype", names(coef))  #location of phenotype
  # out1 = cbind(coef, pvalue)[pheno.index,]
  
  # 2. global test
  bereg.null = try(gamlss(y.prop ~ batch, nu.formula = ~ batch,
                          family = BEZI(sigma.link = "log"),  data = data,
                          control = gamlss.control(n.cyc = 100, trace = FALSE)))
  if (any(class(bereg.null) %in% "try-error")) {bereg.null <- list(G.deviance = NA, df.residual = NA)}
  
  chisq = bereg.null$G.deviance - bereg$G.deviance
  df = bereg.null$df.residual - bereg$df.residual
  out2 = matrix(c(chisq, 1 - pchisq(chisq, df)), 1, 2)
  
  # 3. stack up
  out = rbind(out1, out2)
  rownames(out) = c("LB.nonz", "LB.zero", "LB.glob")
  colnames(out) = c("Estimate", "pval")
  if (test == TRUE){
    print(out)
  }
  return(out)
}


### 1~3. LB
# instead of Boyangs LB.test(), use the previous LB()
# LB uses gamlss package (with log sig.link), while LB.test uses logistic + betareg.

LB <- function (data, fin.n = 0, fin.z = FALSE) {
  # fin.n: minimum nonzero counts for nonzero test
  # fin.z: logical. if true, correct beta hat of logistic regression
  require(gamlss)
  n = dim(data)[1]
  data$y.prop <- data$y / data$sampleSum  #sampleSum[j] = sum_g y_g,j (g:gene, j:sample)
  out <- matrix(NA, 3, 2, dimnames = list(c("LB.nonz", "LB.zero", "LB.glob"), c("Estimate", "pval")))
  drop.nonzero <- sum(data$y > 0, na.rm = TRUE) < fin.n
  
  # 1. two-part models
  if (!drop.nonzero) {
    bereg = try(gamlss(y.prop ~ phenotype + batch, nu.formula = ~ phenotype + batch,
                       family = BEZI(sigma.link = "log"), data = data,
                       control = gamlss.control(n.cyc = 100, trace = FALSE)))
    if (any(class(bereg) %in% "try-error")) return(out)
    # figure out which of mu and nu models are okay: index = either NULL, 1, 2, 1&2.
    index <- NULL
    if (!is.na(bereg$mu.coefficients["phenotypeH"])) {index <- c(1)}
    if (!is.na(bereg$nu.coefficients["phenotypeH"])) {index <- c(index, 2)}
    
    # extracting coef & pval matrix
    bereg.summary <- suppressWarnings(summary.gamlss2(bereg))
    pheno.index <- grep("phenotype", rownames(bereg.summary))
    tab.tmp <- bereg.summary[pheno.index, c("Estimate", "Pr(>|t|)", "Std. Error"), drop=FALSE]
    
    if (dim(tab.tmp)[1] != length(index)) stop("Something is wrong. number of coefficients do not match.")
    
    if (fin.z & (2 %in% index)) { #finite sample correction for zero models
      nu.index <- grep("(Intercept)", rownames(bereg.summary))[3]: dim(bereg.summary)[1]
      tab.nu.tmp <- bereg.summary[nu.index, c("Estimate", "Pr(>|t|)", "Std. Error"), drop=FALSE]
      
      pred <- plogis(bereg$nu.lp)
      X <- bereg$nu.x
      XDX.inv <- solve(t(X) %*% (X * pred))
      Q <- X %*% XDX.inv %*% t(X)
      xi <- diag(Q) * (pred - 0.5)
      bias <- XDX.inv %*% t(X) %*% (pred * (1-pred) * xi)
      
      tab.tmp[2, "Estimate"] <- (tab.nu.tmp[, "Estimate"] - bias)["phenotypeH",]
      tab.tmp[2, "Std. Error"] <- tab.nu.tmp["phenotypeH", "Std. Error"] * n/(n + bereg$nu.df)  # bereg$nu.df = 3 = intercept + Pheno + batch
      tab.tmp[2, "Pr(>|t|)"] <- 2 * pt(-abs(tab.tmp[2, "Estimate"] / tab.tmp[2, "Std. Error"]), df = bereg$df.residual)
    }
    out[index, ] <- tab.tmp[, 1:2]
    # out1 = summary.gamlss2(bereg)[c(2,6), c(1,4)]  # summary.gamlss2 in F01.02.summary.gamlss2.R
    
    
    # 2. global test
    ## LRT is not available when finite sample correction applies.
    # bereg.null = try(gamlss(y.prop ~ batch, nu.formula = ~ batch,
    #                         family = BEZI(sigma.link = "log"),  data = data,
    #                         control = gamlss.control(n.cyc = 100, trace = FALSE)))
    # if (any(class(bereg.null) %in% "try-error")) {bereg.null <- list(G.deviance = NA, df.residual = NA)}
    # 
    # chisq = bereg.null$G.deviance - bereg$G.deviance
    # df = bereg.null$df.residual - bereg$df.residual
    # out["LB.glob", ] = c(chisq, 1 - pchisq(chisq, df))
    chisq = sum((tab.tmp[, "Estimate"] / tab.tmp[, "Std. Error"])^2, na.rm = TRUE)
    df = length(index)
    out["LB.glob", ] = c(chisq, 1 - pchisq(chisq, df))
    
  } else {  # When the number of nonzero data is small.
    
    bereg = try(glm(ifelse(y.prop > 0, 1, 0) ~ phenotype + batch, family = binomial, data = data))
    if (any(class(bereg) %in% "try-error")) return(out)
    tab.tmp <- summary(bereg)$coef
    
    if (fin.z) { #finite sample correction for zero models
      pred <- bereg$fitted.values
      X <- model.matrix(bereg)
      XDX.inv <- solve(t(X) %*% (X * pred))
      Q <- X %*% XDX.inv %*% t(X)
      xi <- diag(Q) * (pred - 0.5)
      bias <- XDX.inv %*% t(X) %*% (pred * (1-pred) * xi)
      
      tab.tmp[, "Estimate"] <- tab.tmp[, "Estimate"] - bias
      tab.tmp[, "Std. Error"] <- tab.tmp[, "Std. Error"] * n/(2 * n - bereg$df.residual)  # df = 3 = intercept + Pheno + batch (n + k = n + (n - df.res))
      tab.tmp[, "Pr(>|z|)"] <- 2 * pt(-abs(tab.tmp[, "Estimate"] / tab.tmp[, "Std. Error"]), df = bereg$df.residual)
    }
    
    out["LB.zero", ] <- tab.tmp["phenotypeH", c("Estimate", "Pr(>|z|)")]
    # out["LB.nonz", ] <- NA
    out["LB.glob", ] <- out["LB.zero", ]
    
  }
  
  return(out)
}
if (FALSE) {
  data = rZINB.sim(n.sample=rep(3,4),n.genes=10, 1,1,1)
  
  a <- LB(data %>% mutate(y=y.1), fin.n = 0)
  b <- LB(data %>% mutate(y=y.1), fin.n = 10)
  
  a2 <- LB(data %>% mutate(y=y.1), fin.n = 0, fin.z = TRUE)
  b2 <- LB(data %>% mutate(y=y.1), fin.n = 10, fin.z = TRUE)
}



### 4. Log-normal
LN <- function (data, sig = 0.05, epsilon = 1) {
  # log-transformation
  data$log2y = log2(data$y + epsilon)
  
  # fitting a linear model
  out = lm(log2y ~ phenotype + batch, data = data)
  out = matrix(summary(out)$coef[2, c(1,4)], nrow = 1)
  colnames(out) = c("Estimate", "pval")
  
  return(out)
}


### 5-7. MAST TBD!!!
MAST <- function (data, sig = 0.05) {
  # devtools::install_github("RGLab/MAST")
  require (MAST)
  
  # whole-data-level test. not adequate for inidivdual-gene-level test.
  #print(303)  
  name = names(data)
  #print(305)
  gene = which(grepl("y\\.", name))
  #print(307)
  gene.name = gsub("y\\.", "", name[gene])
  #print(309)  
  # log-transformation
  data[,gene] = log2(data[,gene] + 1)
  # print(312)
  cData = data %>% transmute(wellKey = 1:n(), phenotype, batch)
  # print(314)  
  data = t(as.matrix(data[,gene]))
  # print(316)  
  dimnames(data) = list(gene.name, cData$wellKey)
  # print(318)  
  sca <- FromMatrix(data, cData = cData, 
                    fData = data.frame(primerid = gene.name))
  # print(321)  
  # Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca) > 0)
  # print(324)  
  colData(sca)$cngeneson <- scale(cdr2)
  # print(326)  
  zlm.out <- try(zlm( ~ phenotype + batch + cngeneson, sca))
  # saveRDS(zlm.out, "zlm.out.rds")
  # print(328)  
  if (any(class(zlm.out) %in% "try-error")) {
    result = matrix(NA, length(gene.name), 3, 
                    dimnames = list(gene.name, c("cont", "disc", "hurdle")))
    
    return(list(coef = result %>% t, p = result %>% t))
  }
  
  # print(337)  
  result <- data.frame(coefC = (zlm.out@coefC)[,2],  # second column = phenotypeH
                       coefD = (zlm.out@coefD)[,2],  # second column = phenotypeH
                       seC = zlm.out@vcovC %>% apply(3,diag) %>% "["(2,) %>% sqrt,
                       seD = zlm.out@vcovD %>% apply(3,diag) %>% "["(2,) %>% sqrt)
  # print(342)
  # saveRDS(result, "result342.rds")
  df = c(1,1,1) # phenotype = 1df
  result %<>% transmute(cont = coefC^2/seC^2,
                        disc = coefD^2/seD^2,
                        hurdle = cont + disc)
  # print(348)
  # saveRDS(result, "result348.rds")
  result.p = transmute(result,
                       pC = 1 - pchisq(cont, df[1]),
                       pD = 1 - pchisq(disc, df[2]),
                       pH = 1 - pchisq(hurdle, df[3]))
  # print(354)
  # saveRDS(result, "result354.rds")
  return(list(coef = result %>% t, p = result.p %>% t))
}


if (FALSE) {# example
  data %>% MAST() -> tmp.a
  waldTest(tmp.a[[1]], Hypothesis('`phenotypeH`'))
  tmp.a@coefC
  
  
  data <- data.frame(x=rnorm(500), z=rbinom(500, 1, .3))
  logit.y <- with(data, x*2 + z*2); mu.y <- with(data, 10+10*x+10*z + rnorm(500))
  y <- (runif(500)<exp(logit.y)/(1+exp(logit.y)))*1
  y[y>0] <- mu.y[y>0]
  data$y <- y
  fit <- zlm(y ~ x+z, data)
  summary.glm(fit$disc)
  summary.glm(fit$cont)
}


### 8. KW
KW <- function (data, sig = 0.05) {
  require(coin)
  # log-transformation not needed for KW
  
  # fitting a nonparametric model
  out = kruskal_test(y ~ phenotype | batch, data = data)
  out = matrix(c(statistic(out), pvalue(out)), nrow = 1)
  colnames(out) = c("Estimate", "pval")
  return(out)
}

### 9. Wagner
Wagner <- function (data, sig = 0.05, zeroModel = c("logistic", "t.test", "lm"),
                    suppressWarning = FALSE) {
  
  # Instead of proportion t-test, logistic regression is used to adjust for batch effect.
  # For nonzero model, modified WRS test (in coin package) is used. But it cannot handle small nonzero sample.
  
  # 1. zero model
  if (zeroModel[1] == "logistic") { # batch adjusted
    data.bin = data %>% mutate(y = ifelse(y>0, 1, 0))
    if (suppressWarning) {
      Z = suppressWarnings(glm(y~phenotype + factor(batch), data=data.bin, family="binomial") %>% summary)$coef[2,3:4] %>% as.numeric
    } else {
      Z = (glm(y~phenotype + factor(batch), data=data.bin, family="binomial") %>% summary)$coef[2,3:4] %>% as.numeric
    }
    
    #if phat = 0 or 1 (or all D and H are 0 or 1), then Z is automatically close to 0.
  } else if (zeroModel[1] == "t.test") {  # batch not adjusted!!!
    # table
    data %>% group_by(phenotype) %>% summarize(n = n(), n1 = sum(y!=0)) -> cont
    nD  = cont[1,2];  nH  = cont[2,2]
    nD1 = cont[1,3];  nH1 = cont[2,3]
    
    Z = abs(nD1/nD - nH1/nH) - .5/nD - .5/nH
    p.hat = (nD1 + nH1) / (nD + nH)
    Z = (Z / (p.hat*(1-p.hat)*(1/nD + 1/nH))^.5) %>% as.numeric
    
    if (p.hat*(1 - p.hat) == 0) {Z = 0}
    names(Z) <- NULL
    Z = c(Z, 2 - 2*pnorm(abs(Z %>% as.numeric)))
  } else if (zeroModel[1] == "lm") { # Just for comparison!!
    data.bin = data %>% mutate(y = ifelse(y>0, 1, 0))
    if (suppressWarning) {
      Z = suppressWarnings(lm(y~phenotype + factor(batch), data=data.bin) %>% summary)$coef[2,3:4, drop=FALSE]
    } else {
      Z = (lm(y~phenotype + factor(batch), data=data.bin) %>% summary)$coef[2,3:4, drop=FALSE]
    }
    
  } else stop("zeroModel was not correctly specified.")
  
  
  
  # return(c(z = Z, p = (1- pnorm(abs(as.numeric(Z))))*2))  
  
  # 2. nonzero model  
  # nonzero data
  data %>% dplyr::filter(y != 0) -> data.nonzero
  data.nonzero %>% dplyr::filter(phenotype=="D") %>% "$"("y") -> y.D
  data.nonzero %>% dplyr::filter(phenotype=="H") %>% "$"("y") -> y.H
  #print(y.D)
  
  # WRS test  
  # out = wilcox.test(x = y.D, y = y.H,
  #             alternative = c("two.sided"), correct = TRUE)
  # print(data.nonzero)
  # print("1. original WRS without batch")
  # wilcox.test(x = y.D, y = y.H,
  #             alternative = c("two.sided"), exact = FALSE, correct = TRUE) %>% print
  # print("2. modfified WRS without batch")
  # coin::wilcox_test(y ~ factor(phenotype), data=data.nonzero) %>% print
  
  # print("4. t-test WITH batch")
  # lm(y~phenotype+batch, data=data.nonzero) %>% print
  
  # print("3. modfified WRS WITH batch")
  W = try(coin::wilcox_test(y ~ factor(phenotype) | factor(batch), data=data.nonzero), silent = suppressWarning) 
  # silent=FALSE causes the slurm-output file size to explode.
  if (any(class(W) %in% "try-error")) {
    W = matrix(NA, 1, 2)
  } else {
    W = matrix(c(statistic(W), pvalue(W)), nrow = 1)
  }
  
  colnames(W) = c("Estimate", "pval")
  
  # put together
  # print(Z)
  # print(W)
  chi2 = c(Z^2 + ifelse(is.na(W), 0, W)^2)
  chi2 = matrix(c(chi2, 1-pchisq(chi2, df=2)),1,2)
  
  out = rbind(W, Z, chi2) # nonzero, zero, global
  rownames(out) = c("Wg.nonz", "Wg.zero", "Wg.glob")
  rownames(out) = c("LB.nonz", "LB.zero", "LB.glob")
  
  return(out)
}

if (FALSE) {# example
  data %>% mutate(y=y.1) %>% Wagner()
  Wagner(data.1)
  Wagner(data.2)
  
  data %>% mutate(y=y.1) %>% Wagner(zeroModel = "logistic")
  data %>% mutate(y=y.1) %>% Wagner(zeroModel = "t.test")
  data %>% mutate(y=y.1) %>% Wagner(zeroModel = "lm")
  glm(y.1~phenotype, data=data) %>% summary
}

### 12. DESeq2
DS2 <- function (data.l, sig = 0.05) {
  require(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = round(t(data.l[, grepl("^y", names(data.l))]),0),
                                colData = data.l[, !grepl("^y", names(data.l))],
                                design= ~ batch + phenotype)
  cts = counts(dds)
  # handle the case where all genes have at least one zero.
  geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds2 = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds3 <- DESeq(dds2)
  results(dds3)
  res <- results(dds3, name="phenotype_H_vs_D")
  out = matrix(c(res$log2FoldChange, res$pvalue), ncol = 2)
  colnames(out) = c("Estimate", "pval")
  return(out)
}
