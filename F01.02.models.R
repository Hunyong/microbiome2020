#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(BiocManager)
n.test = 22 #"LB.nonz", "LB.zero", "LB.glob", "LB.min", "LN", "MAST.nonz", "MAST.zero", "MAST.glob", "MAST.min", 
#"KW", "Wg.nonz", "Wg.zero", "Wg.glob", "Wg.min", "DS2", "DS2ZI, "WRS", "MGS", "ANCOM.sz", "ANCOM", "LFE", "ALDEX"
# if(!require(betareg)){
#   install.packages("betareg")
# }


### 0. testing a set of methods at once
# testing all methods using counts (or RPK) of a single gene (old: single gene. cannot do MAST).
# HD: binary phenotype (healthy-diseased)
tester.set.HD.batch <- function(data, n.gene = 10000, 
                                pvalue.only = FALSE, skeleton = FALSE, 
                                suppressWarnWagner = FALSE, LB.skip = FALSE, 
                                LN.skip = FALSE, MAST.skip = FALSE,
                                KW.skip = FALSE, Wg.skip = FALSE,
                                De2.skip = FALSE, WRS.skip = FALSE,
                                MGS.skip = FALSE, ANCOM.skip = FALSE,
                                LEfSe.skip = FALSE, ALDEX.skip = FALSE,
                                SONGBIRD.skip = FALSE,
                                skip.cumulative = FALSE,
                                skip.small.n = FALSE) {
  # DS2.version = c("vanilla", "zinb")
  # description
  # data should have y and sampleSum    all n.sample x (n.gene(gene) + 3 (phenotype + batch + sampleSum))
  #          outcome (phenotype), nuisance (batch)
  # skeleton: returning only skeleton (for simulation structure)
  # Sometime in 2019, DEseq2 test was added.
  # On May 20, 2020: WRS test has been added.
  require(magrittr)
  
  if (skip.cumulative) { # IF skip.cumulative = TRUE and specify and a method is specified as skipped, all its preceding methods are to be skipped.
    if (SONGBIRD.skip) ALDEX.skip = TRUE
    if (ALDEX.skip) LEfSe.skip = TRUE
    if (LEfSe.skip) ANCOM.skip = TRUE
    if (ANCOM.skip) MGS.skip = TRUE
    if (MGS.skip) WRS.skip = TRUE
    if (WRS.skip) De2.skip = TRUE
    if (De2.skip) Wg.skip = TRUE
    if (Wg.skip) KW.skip = TRUE
    if (KW.skip) MAST.skip = TRUE
    if (MAST.skip) LN.skip = TRUE
    if (LN.skip) LB.skip = TRUE
  }
  
  # 0.0 skeleton #empty matrix
  test.names <- c("LB.nonz", "LB.zero", "LB.glob", "LB.min", "LN", "MAST.nonz", "MAST.zero", "MAST.glob", "MAST.min", 
                  "KW", "Wg.nonz", "Wg.zero", "Wg.glob", "Wg.min", "DS2", "DS2ZI", "WRS", "MGS", "ANCOM.sz", "ANCOM", "LFE", "ALDEX")
  mat.tmp <- matrix(NA, n.test, n.gene, dimnames = list(test.names, NULL)) # n.test=15, n.gene=1000
  result <- list(coef = mat.tmp, pval = mat.tmp)
  if (skeleton) {return(result)}
  
  # 0.1 data
  data2 <- data
  genes <- gsub("^y\\.", "", names(data)) %>% as.numeric %>% na.omit %>% as.numeric
  if (n.gene > genes[length(genes)]) stop(paste0("Only ", length(genes), " genes provided, while trying to do ", n.gene, " simulations."))
  genes = genes[1:min(n.gene,length(genes))]
  
  index.filtered = 
    apply(data[, genes], 2, function(x) all(is.na(x))) %>% # filtered gene indices
    {which(!.)} %>% as.numeric
  index.meta = grepl("^y\\.", names(data)) %>% "!"(.) %>% which
  index.filtered.meta = c(index.filtered, index.meta)
  if (length(index.filtered) < 10) return(result)
  
  data = data.frame(data[,1:length(genes)], phenotype = data$phenotype, batch = data$batch)
  ## print(head(data))    (No longer need to check data)
  if (!"sampleSum" %in% names(data2)) {
    data$sampleSum <- data %>% dplyr::select(-phenotype, - batch) %>% apply(1, sum, na.rm = TRUE)
  } else {
    data$sampleSum <- data2$sampleSum
  }
  
  
  cat("1-3. Logistic Beta\n")
  #1-3. LB
  if (!LB.skip) {
    for (l in genes) {
      if (l %% 30 == 0) cat (" l = ",l," ")
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y, na.rm = TRUE)==0) {
        tmp <- data.frame(coef = rep(NA,3), pval = NA)
      } else {
        ## print(data[,l]) (For debug only)
        tmp <- LB(data.l)  #logistic beta
      }
      result[[1]][c("LB.nonz", "LB.zero", "LB.glob"), l] <- tmp[1:3, 1] #coef. 1:3 corresponds to "LB.nonz", "LB.zero", "LB.glob"
      result[[2]][c("LB.nonz", "LB.zero", "LB.glob"), l] <- tmp[1:3, 2] #pval. 1:3 corresponds to "LB.nonz", "LB.zero", "LB.glob"
      if (any(!is.na(tmp[1:2, 2]))) result[[2]]["LB.min", l] <- min(tmp[1:2, 2], na.rm = TRUE)
    }
  } else {cat("LB is skipped\n")}
  
  
  #4. LN
  
  cat("\n4. Log normal\n l = ")
  if (!LN.skip) {
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y, na.rm = TRUE)==0) {
        tmp <- data.frame(coef = NA, pval = NA)
      } else {
        tmp <- LN(data.l)  #log normal
      }
      result[[1]]["LN", l] <- tmp[1,1] #coef.
      result[[2]]["LN", l] <- tmp[1,2] #pval.
    }
  } else {cat("LN is skipped\n")}
  
  #5-7. MAST
  cat("\n5-7. MAST\n")
  if(!MAST.skip){
    # tmp <- MAST(data)  #MAST
    # tmp <- data.frame(coef = rep(NA,3), pval = NA)    #MAST maybe not applicable
    tmp <- try({MAST(data[, index.filtered.meta])})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
    result[[1]][c("MAST.nonz", "MAST.zero", "MAST.glob"), index.filtered] <- tmp[[1]][1:3,] #coef. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
    result[[2]][c("MAST.nonz", "MAST.zero", "MAST.glob"), index.filtered] <- tmp[[2]][1:3,] #pval. 1:3 corresponds to "MA.nonz", "MA.zero", "MA.glob"
    result[[2]]["MAST.min", index.filtered] <- pmin(tmp[[2]][1,], tmp[[2]][2,], na.rm = TRUE)
  } else {cat("MAST is skipped\n")}
  
  #8. KW
  cat("\n8. Kruskal Wallis\n l = ")
  if (!KW.skip){
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y, na.rm = TRUE)==0) {
        tmp <- data.frame(coef = NA, pval = NA)
      } else {
        tmp <- KW(data.l)  #KW
      }
      result[[1]]["KW", l] <- tmp[1,1] #coef.
      result[[2]]["KW", l] <- tmp[1,2] #pval.
    }
  } else {cat("KW is skipped\n")}
  
  #9-11. Two-part KW
  cat("\n9-11. Two-part KW (Wagner)\n l = ")
  if(!Wg.skip){
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y, na.rm = TRUE)==0) {
        tmp <- data.frame(coef = rep(NA,3), pval = NA)
      } else {
        tmp <- Wagner(data.l, zeroModel = "logistic", suppressWarning = suppressWarnWagner)
      }
      result[[1]][c("Wg.nonz", "Wg.zero", "Wg.glob"), l] <- tmp[1:3,1] #coef.
      result[[2]][c("Wg.nonz", "Wg.zero", "Wg.glob"), l] <- tmp[1:3,2] #pval.
      if (any(!is.na(tmp[1:2, 2]))) result[[2]]["Wg.min", l] <- min(tmp[1:2, 2], na.rm = TRUE)
    }
  } else {cat("Wg is skipped\n")}
  
  #12. De2
  cat("\n12. DESeq2 vanilla\n")
  if(!De2.skip){
    DS2 = DS2.vanilla
    tmp <- try({DS2(data[, index.filtered.meta])})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
    result[[1]]["DS2", index.filtered] <- tmp[,1] #coef.
    result[[2]]["DS2", index.filtered] <- tmp[,2] #pval.
  } else {cat("DS2-vanilla is skipped\n")}
  
  #13. De2 + ZINBwave
  cat("\n13. DS2 zinbwave\n")
  if(!De2.skip){
    DS2 = DS2.zinb
    tmp <- try({DS2(data[, index.filtered.meta])})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
    result[[1]]["DS2ZI", index.filtered] <- tmp[,1] #coef.
    result[[2]]["DS2ZI", index.filtered] <- tmp[,2] #pval.
  } else {cat("DS2-zinbwave is skipped\n")}
  
  cat("\n14. WRS\n l = ")
  if (!WRS.skip){
    for (l in genes) {
      # cat (l," ")
      if (l %% 200 == 0) {cat(l, " ")}
      data.l = data.frame(y=data[,l], data[,c("phenotype", "batch", "sampleSum")])
      if (sum(data.l$y, na.rm = TRUE)==0) {
        tmp <- data.frame(coef = NA, pval = NA)
      } else {
        tmp <- WRS(data.l)  #KW
      }
      result[[1]]["WRS", l] <- tmp[1,1] #coef.
      result[[2]]["WRS", l] <- tmp[1,2] #pval.
    }
  } else {cat("WRS is skipped\n")}
  
  #15. metagenomeSeq
  cat("\n15 metagenomeSeq\n")
  if(!MGS.skip){
    tmp <- try({mgs(data[, index.filtered.meta])})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
    
    result[[1]]["MGS", index.filtered] <- tmp[, "Estimate"] #coef.
    result[[2]]["MGS", index.filtered] <- tmp[, "pval"]     #pval.
  } else {cat("MGS is skipped\n")}
  
  #16. ANCOM-BC
  cat("16. ANCOM-BC\n")
  if(!ANCOM.skip){
    tmp <- try({ANC(data[, index.filtered.meta], ignore.structural.zero = FALSE)})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
    
    result[[1]]["ANCOM.sz", index.filtered] <- tmp[,1] #coef.
    result[[2]]["ANCOM.sz", index.filtered] <- tmp[,2] #pval.
    result[[1]]["ANCOM", index.filtered] <- tmp[,1] #coef.
    result[[2]]["ANCOM", index.filtered] <- tmp[,2] #pval.
    result[[2]]["ANCOM", is.infinite(result[[1]]["ANCOM", ])] = NA
  } else {cat("ANCOM-BC is skipped\n")}
  
  #17. LEfSe
  cat("17. LEfSe\n")
  if(!LEfSe.skip){
    tmp <- try({LEfSe(data[, index.filtered.meta])})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
    
    result[[1]]["LFE", index.filtered] <- tmp[,1] #coef.
    result[[2]]["LFE", index.filtered] <- tmp[,2] #rank/n.
  } else {cat("LEfSe is skipped\n")}
  
  #18. ALDEX2
  cat("18. ALDEX2\n")
  if(!ALDEX.skip){
    tmp <- try({aldex(data[, index.filtered.meta])})
    if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)

    result[[1]]["ALDEX", index.filtered] <- tmp[,1] #coef.
    result[[2]]["ALDEX", index.filtered] <- tmp[,2] #pval.
  } else {cat("ALDEX2 is skipped\n")}
  
  # #19. SONGBIRD
  # cat("19. SONGBIRD\n")
  # if(!SONGBIRD.skip){
  #   tmp <- try({songbird(data[, index.filtered.meta])})
  #   if (class(tmp)[1] == "try-error") tmp = matrix(NA, ncol = 2)
  #   
  #   result[[1]]["SBIRD", index.filtered] <- tmp[,1] #coef.
  #   result[[2]]["SBIRD", index.filtered] <- tmp[,2] #pval.
  # } else {cat("SONGBIRD is skipped\n")}
  
  return(result)
}


if (FALSE) {# examples
  data = rZINB.sim(n.sample=rep(3,4),n.genes=30, 1,1,1)
  data %>% dplyr::select(-phenotype, - batch) %>% apply(1, sum) -> data$sampleSum
  data.1 = data.frame(y=data[,1], data[,c("phenotype", "batch", "sampleSum")])
  a = tester.set.HD.batch(data, n.gene=30)
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
LN <- function (data, epsilon = 1) {
  # log-transformation
  data$log2y = log2(data$y + epsilon)
  
  # fitting a linear model
  out = lm(log2y ~ phenotype + batch, data = data)
  out = matrix(summary(out)$coef[2, c(1,4)], nrow = 1)
  colnames(out) = c("Estimate", "pval")
  
  return(out)
}


### 5-7. MAST TBD!!!
MAST <- function (data) {
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
KW <- function (data) {
  require(coin)
  # log-transformation not needed for KW
  
  # fitting a nonparametric model
  out = kruskal_test(y ~ phenotype | batch, data = data)
  out = matrix(c(statistic(out), pvalue(out)), nrow = 1)
  colnames(out) = c("Estimate", "pval")
  return(out)
}

### 9. Wagner
Wagner <- function (data, zeroModel = c("logistic", "t.test", "lm"),
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
  # print("1. original KW without batch")
  # kruskal.test(x = y.D, y = y.H,
  #             alternative = c("two.sided"), exact = FALSE, correct = TRUE) %>% print
  # print("2. modfified KW without batch")
  # coin::kruskal_test(y ~ factor(phenotype), data=data.nonzero) %>% print
  
  # print("4. t-test WITH batch")
  # lm(y~phenotype+batch, data=data.nonzero) %>% print
  
  # print("3. modfified KW WITH batch")
  W = try(coin::kruskal_test(y ~ factor(phenotype) | factor(batch), data=data.nonzero), silent = suppressWarning) 
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
  chi2 = Z[1]^2 + ifelse(is.na(W[1]), 0, W[1]) # W is a squared value already
  chi2 = matrix(c(Estimate = chi2, pval = 1-pchisq(chi2, df = 2 - is.na(Z[1]) - is.na(W[1]))), 
                1, 2)
  
  out = rbind(W, Z, chi2) # nonzero, zero, global
  rownames(out) = c("Wg.nonz", "Wg.zero", "Wg.glob")
  
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

### 12. DS2 ZINB-WAVE extension
DS2.zinb <- function (data.l) {
  ### Much part of this code is from 
  # https://github.com/mikelove/zinbwave-deseq2/blob/master/zinbwave-deseq2.knit.md
  # With lots of fatal errors from scran, this code does not have the scran components.
  
  require(DESeq2)
  require(zinbwave)
  require(BiocParallel)
  
  ### Extra filter for zero counts out too sparse genes
  col.otu = which(grepl("^y", names(data.l)))
  col.meta = which(!grepl("^y", names(data.l)))
  keepForTests <- colSums(round(data.l[, col.otu], 0) >= 1) >= 3  # This will be used again at the end.
  data.l <- data.l[, c(col.otu[keepForTests], col.meta)]
  
  ### Getting the ZINB-wave weights
  col.otu = which(grepl("^y", names(data.l))) # updated
  zinb <- DESeqDataSetFromMatrix(countData = round(t(data.l[, col.otu]),0),
                                 colData = data.l[, !grepl("^y", names(data.l))],
                                 design= ~ batch + phenotype)
  
  # # we need to reorganize the assays in the SumExp from splatter
  assay(zinb) <- as.matrix(assay(zinb))
  X = data.matrix(data.l[, c("phenotype", "batch")])
  
  # epsilon setting as recommended by the ZINB-WaVE integration paper
  zinb <- zinbwave(zinb, X = X, K=0, observationalWeights=TRUE,
                   BPPARAM=BiocParallel::SerialParam(), epsilon=1e12)
  
  dds <- DESeqDataSet(zinb, design= ~ batch + phenotype)
  
  dds <- DESeq(dds, sfType = "poscounts", test = "LRT", reduced = ~batch, minmu=1e-6, minRep=Inf)
  
  res <- results(dds, independentFiltering = FALSE, name = "phenotype_H_vs_D")
  
  # insert filtered out NA results
  out = matrix(NA, nrow = length(keepForTests), ncol = 2)
  out[keepForTests, 1] = res$log2FoldChange
  out[keepForTests, 2] = res$pvalue
  # out = matrix(c(res$log2FoldChange, res$pvalue), ncol = 2)
  colnames(out) = c("Estimate", "pval")
  return(out)
}

### 12. DESeq2 ZINB-WAVE extension (old version with scran)
DS2.zinb.old <- function (data.l) {
  tmp.dat <<- data.l
  saveRDS(tmp.dat, "tmpdat.rds")
  ### Much part of this code is from 
  # https://github.com/mikelove/zinbwave-deseq2/blob/master/zinbwave-deseq2.knit.md
  
  require(DESeq2)
  require(zinbwave)
  require(scran)
  require(BiocParallel)
  
  ### Extra filter for zero counts out too sparse genes
  col.otu = which(grepl("^y", names(data.l)))
  col.meta = which(!grepl("^y", names(data.l)))
  keepForTests <- colSums(round(data.l[, col.otu], 0) >= 1) >= 3  # This will be used again at the end.
  data.l <- data.l[, c(col.otu[keepForTests], col.meta)]
  
  ### Getting the ZINB-wave weights
  col.otu = which(grepl("^y", names(data.l))) # updated
  zinb <- DESeqDataSetFromMatrix(countData = round(t(data.l[, col.otu]),0),
                                 colData = data.l[, !grepl("^y", names(data.l))],
                                 design= ~ batch + phenotype)
  
  # # we need to reorganize the assays in the SumExp from splatter
  # nms <- c("counts", setdiff(assayNames(zinb), "counts"))
  # assays(zinb) <- assays(zinb)[nms]
  assay(zinb) <- as.matrix(assay(zinb))
  X = data.matrix(data.l[, c("phenotype", "batch")])
  # epsilon setting as recommended by the ZINB-WaVE integration paper
  zinb <- zinbwave(zinb, X = X, K=0, observationalWeights=TRUE,
                   BPPARAM=BiocParallel::SerialParam(), epsilon=1e12)
  
  dds <- DESeqDataSet(zinb, design= ~ batch + phenotype)
  
  # Use size factors from the scran package - Lots of error
  scr <- try({computeSumFactors(dds)})
  if (class(scr)[1] == "try-error") {
    # negative size factor is calculated. Filtering.
    keepForSizeComp <- rowSums(counts(dds) >= 3) >= 10
    if (sum(keepForSizeComp) < 30) keepForSizeComp <- rowSums(counts(dds) >= 3) >= 5
    scr <- try({computeSumFactors(dds[keepForSizeComp, ])})
    if (class(scr)[1] == "try-error")
      return(matrix(NA, ncol = 2, nrow = length(keepForTests), dimnames = list(NULL, c("Estimate", "pval"))))
  }
  sizeFactors(dds) <- sizeFactors(scr)
  
  dds <- try({DESeq(dds, test = "LRT", reduced = ~batch, minmu=1e-6, minRep=Inf)})
  if (class(dds)[1] == "try-error")
    return(matrix(NA, ncol = 2, nrow = length(keepForTests), dimnames = list(NULL, c("Estimate", "pval"))))
  # plotDispEsts(dds)
  
  keepForDispTrend <- rowSums(counts(dds) >= 3) >= 10
  # If the filtering results in too few samples, relax the threshold. If the alternative does not work, return NA.
  # if (sum(keepForDispTrend) < 30) keepForDispTrend <- rowSums(counts(dds) >= 3) >= 5 # However, this lower threshold causes a fatal error when DESe2 is implemented.
  if (sum(keepForDispTrend) < 30) {
    warning("Dispersion is not stably estimated. NA's are returned.")
    return(matrix(NA, ncol = 2, nrow = length(keepForTests), dimnames = list(NULL, c("Estimate", "pval"))))
  }
  dds2 <- estimateDispersionsFit(dds[keepForDispTrend,], fitType = "local")
  # parametric models still fail. Use the local estimation.
  # plotDispEsts(dds2, ylim=c(1e-3,1))
  
  dispersionFunction(dds) <- dispersionFunction(dds2)
  dds <- estimateDispersionsMAP(dds)
  dds <- nbinomLRT(dds, reduced = ~batch, minmu=1e-6)
  
  # filtering is already done.
  res <- results(dds, independentFiltering = FALSE, name = "phenotype_H_vs_D")
  
  # insert filtered out NA results
  out = matrix(NA, nrow = length(keepForTests), ncol = 2)
  out[keepForTests, 1] = res$log2FoldChange
  out[keepForTests, 2] = res$pvalue
  # out = matrix(c(res$log2FoldChange, res$pvalue), ncol = 2)
  colnames(out) = c("Estimate", "pval")
  return(out)
}
### 12. DESeq2 plain version.
DS2.vanilla <- function (data.l) {
  tmp.dat <<- data.l
  saveRDS(tmp.dat, "tmp.dat")  
  require(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = round(t(data.l[, grepl("^y", names(data.l))]),0),
                                colData = data.l[, !grepl("^y", names(data.l))],
                                design= ~ batch + phenotype)
  cts = counts(dds)
  # handle the case where all genes have at least one zero.
  geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds2 = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds3 <- try(DESeq(dds2))
  if (class(dds3)[1] == "try-error") {
    dds2 <- estimateDispersionsGeneEst(dds2)
    dispersions(dds2) <- mcols(dds2)$dispGeneEst
    dds3 <- try(DESeq(dds2))
    if (class(dds3)[1] == "try-error") {
      return(as.matrix(data.frame(Estimate = NA, pval = NA)))
    }
  }
  results(dds3)
  res <- results(dds3, name="phenotype_H_vs_D")
  out = matrix(c(res$log2FoldChange, res$pvalue), ncol = 2)
  colnames(out) = c("Estimate", "pval")
  return(out)
}


### 13. WRS
WRS <- function (data) {
  # fitting a nonparametric model
  out = wilcox.test(y ~ phenotype, data = data, exact = FALSE)
  out = matrix(c(out$statistic, out$p.value), nrow = 1)
  colnames(out) = c("Estimate", "pval")
  return(out)
}

### 14. metagenomeSeq
mgs.base <- function (data) {
  
  require (metagenomeSeq)
  
  # whole-data-level test. not adequate for inidivdual-gene-level test.
  #print(303)  
  name = names(data)
  #print(305)
  gene = which(grepl("y\\.", name))
  #print(307)
  gene.name = gsub("y\\.", "", name[gene])
  #print(309)  
  
  cData = data %>% transmute(samples = 1:n(), phenotype, batch)
  rownames(cData) = 1:dim(cData)[1]
  data = t(as.matrix(data[,gene]))
  dimnames(data) = list(gene.name, cData$samples)
  obj = newMRexperiment(counts = data, 
                        phenoData = AnnotatedDataFrame(cData), 
                        featureData = AnnotatedDataFrame(data.frame(primerid = gene.name, 
                                                                    row.names = gene.name)))
  datp = cumNormStat(data, pFlag = TRUE, main = "Trimmed lung data")
  obj = cumNorm(obj, p = datp)
  normFactors(obj)
  mod <- model.matrix(~ 1 + phenotype, data = cData)
  # mod <- model.matrix(~ 1 + phenotype + batch, data = cData) # "Can't analyze currently."
  mgsRes = fitFeatureModel(obj, mod)
  mgsRes = cbind(mgsRes@fitZeroLogNormal$logFC, mgsRes@pvalues)
  colnames(mgsRes) = c("Estimate", "pval")
  
  return(mgsRes)
}
mgs <- function(data) {
  tmp <- try(mgs.base(data))
  if (class(tmp)[1] == "try-error") {
    name = names(data)
    #print(305)
    gene = which(grepl("y\\.", name))
    gene.name = gsub("y\\.", "", name[gene])
    return(matrix(NA, length(gene), 2, dimnames= list(gene.name, c("Estimate", "pval"))))
  }
  tmp
}

ANC <- function(data, ignore.structural.zero = FALSE) {
  require(nloptr)
  col.otu = which(grepl("^y", names(data)))
  y.names = colnames(data)[col.otu]
  pre.process = 
    feature_table_pre_process (feature.table = t(data[, col.otu]), 
                               meta.data = cbind(ID = rownames(data), data[, -col.otu]), 
                               sample.var = "ID", group.var = "phenotype", 
                               zero.cut = 1, # This is ignored as we use the already-filtered data
                               #lib.cut = 100, # This may not play a role as well by design.
                               neg.lb = FALSE)
  out.ancom = 
    ANCOM_BC (feature.table = pre.process$feature.table, grp.name = pre.process$group.name, 
              grp.ind = pre.process$group.ind, struc.zero = pre.process$structure.zeros, 
              adj.method = "bonferroni", tol.EM = 1e-5, max.iterNum = 100, perNum = 1000, alpha = 0.05)
  out = matrix(NA, nrow = length(y.names), ncol = 2, dimnames = list(y.names, c("Estimate", "pval")))
  out[rownames(out.ancom$feature.table), 1] = out.ancom$res$W
  out[rownames(out.ancom$feature.table), 2] = out.ancom$res$p.val
  
  if (ignore.structural.zero) {out[is.infinite(out[, 1]), 2] = NA}
  return(out)
}


### Place holders for the new tests
LEfSe = function(data) {
  require(SummarizedExperiment)
  require(lefser)
  require(dplyr)
  name = names(data)
  gene = which(grepl("y\\.", name))
  gene.name = gsub("y\\.", "", name[gene])
  cData = data %>% transmute(samples = 1:n(), phenotype, batch)
  data = t(as.matrix(data[,gene]))
  row.names(data) = gene.name
  colData <- DataFrame(phenotype = as.character(cData$phenotype),
                       batch = as.character(cData$batch),
                       row.names = as.character(cData$samples))
  se <- SummarizedExperiment(assays=list(data = data),
                             colData=colData)
  res_lefse = lefser(se, groupCol = "phenotype", blockCol = "batch", 
                     kruskal.threshold = 0.5, wilcox.threshold = 0.5, lda.threshold = 0)
  res_lefse$rank = nrow(res_lefse) + 1 - rank(abs(res_lefse$scores))
  res_lefse$Names = res_lefse$Names %>% gsub("`", "", .) %>% as.character
  merged = left_join(data.frame(gene.name), res_lefse, by= c("gene.name" = "Names"))
  # merged$rank[which(is.na(merged$rank))] = max(merged$rank, na.rm = T) + 1
  # merged$rank[which(is.na(merged$rank))] = (max(merged$rank, na.rm = T) + length(merged$rank))/2
  merged$rank[which(is.na(merged$rank))] = NA # No rank provided for the rest.
  # ..... do tests here
  out = matrix(NA, nrow = length(gene.name), ncol = 2, dimnames = list(gene.name, c("Estimate", "pval")))
  out[, 1] = merged$scores # Insert the coefficients of n genes here!
  # For LEfSe, in the p-value slot, add the rank / the number of genes. (the top genes would get 1/n.genes)
  out[, 2] = merged$rank/length(gene.name) # Insert the p-values of n genes here!
  return(out)
}

aldex <- function(data) {
  library(dplyr)
  library(ALDEx2)
  name <- names(data)
  gene <- which(grepl("y\\.", name))
  gene.name <- gsub("y\\.", "", name[gene])
  cData <- data %>% transmute(samples = 1:n(), phenotype, batch)
  data <- t(as.matrix(data[, gene]))
  data <- round(data)
  row.names(data) <- gene.name

  covariates <- data.frame("phenotype" = cData$phenotype,
                          "batch" = cData$batch)
  mm <- model.matrix(~ phenotype + batch, covariates)

  x <- aldex.clr(data, mm, mc.samples=128, denom="all")
  glm.test <- aldex.glm(x)
  pval <- glm.test$`model.phenotypeH Pr(>|t|)`

  out <- matrix(NA, nrow = length(gene.name), ncol = 2, dimnames = list(gene.name, c("Estimate", "pval")))
  out[, 1] <- NA
  out[, 2] <- pval
  return(out)
}

songbird = function(data) {
  # ..... do tests here
  out = matrix(NA, nrow = length(y.names), ncol = 2, dimnames = list(y.names, c("Estimate", "pval")))
  out[, 1] = NA # Insert the coefficients of n genes here!
  out[, 2] = NA # Insert the p-values of n genes here!
  return(out)
}
