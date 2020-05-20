############################################################################################
### 0. library
  library(dplyr); library(magrittr); library(ggplot2)
  source("F00.00.generic.R")
  source("F01.01.base.R")
  source("F02.01.simulation.R")
  source("F01.02.models-base.R")
  source("F01.02.models.R")
  # devtools::install_github("RGLab/MAST");
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("limma")
  library(MAST); library(coin); library(gamlss)
  library(limma) #venn
  
  gene.marginal.RPK.RNA <- readRDS("../kimon/data/gene.marginal.RPK.RNA.118.rds") #server
  # gene.marginal.RPK.RNA <- readRDS("../code118/output118/gene.marginal.RPK.RNA.118.rds") #local
  outcome.118 <- readRDS("../kimon/data/outcome118.rds") #server
  # outcome.118 <- readRDS("../code118/output118/outcome118.rds") #local
  outcome.118$RNA.date %<>% as.factor
  outcome.118$ECC.bin = factor(ifelse(outcome.118$ECC==0,0,1))
  # screened genes and paths (420 and 352 indexed as 35 and 55)
  gene.marginal.RPK.RNA %>% names %>% "%in%"(c("RPK.420", "RPK.352"))  %>% which -> sample.exclude
  
  gene.index = screen.gene(gene.marginal.RPK.RNA[,-sample.exclude], gene.col=1)$feature$index
  gene.list = gene.marginal.RPK.RNA[gene.index,1]
  
  # min values
  (gene.marginal.RPK.RNA[,-1] %>% as.matrix %>% as.vector) -> a
  a[a>0] %>% sort %>% head  # 0.01927774 0.02180814 0.02209908
  
  # sample sum for BEZI
  sampleSum = apply(gene.marginal.RPK.RNA[,-1], 2, sum)
############################################################################################
  args = commandArgs(trailingOnly=TRUE)  # passed from Script_Boyang
  method = args[1] %>% as.numeric  # 1..5
  k = args[2] %>% as.numeric  # k = 1..11
print(c(method=method, k = k))

### 1. LM
  if (FALSE) {
    j <<- 1; cat("j=")
    gene.lm <- sapply(gene.index, function(i) {
      cat(j,".." ); j <<- j+1
      tmp <- data.frame(log.gene = log2(as.numeric(gene.marginal.RPK.RNA[i,-1]) + 0.019), id = outcome.118$id, 
                        phenotype = outcome.118$ECC.bin, batch = outcome.118$RNA.date)
      tmp <- tmp[!tmp$id %in% c(352, 420),]
      full <- lm(log.gene ~ phenotype  + batch, data=tmp)
      full %>% summary %>% coef %>% "[" (,c(1,4)) %>% as.vector
    }) # 6 x n.path matrix (coef: int + pheno + batch, pval: int + pheno + batch)
    gene.lm <- list ( tested = gene.marginal.RPK.RNA[gene.index,1],
                      coef = matrix(t(gene.lm[1:3,]), ncol=3, dimnames = list(gene.index, c("int", "phenotype", "batch"))),
                      pval = matrix(t(gene.lm[4:6,]), ncol=3, dimnames = list(gene.index, c("int", "phenotype", "batch"))))
    
    p.adjust(gene.lm$pval[,2], method = "BH") %>% 
      matrix(ncol=1, dimnames = list(gene.index, c("phenotype"))) %>% 
      data.frame() -> 
      gene.lm$pval.adj
    
    # 
    gene.lm$sigFeature = gene.list[gene.lm$pval.adj$phenotype < .05]
    
    
    saveRDS(gene.lm, "output/R03.01.LM.R.116.bin.rds")
  } else {
    # gene.lm <- readRDS("output/R03.01.LM.R.116.bin.rds")
  }


### 2. BEZI
  if (method==2) {
    j <<- 1; cat("j=")
    gene.index.k = gene.index[(k-1)*5000 + (1:5000)]
    gene.index.k = gene.index.k[!is.na(gene.index.k)]
print(length(gene.index.k))
    gene.BEZI <- sapply(gene.index.k, function(i) {
      cat(j,".." ); j <<- j+1
      tmp <- data.frame(y = gene.marginal.RPK.RNA[i,-1] %>% as.numeric, id = outcome.118$id, 
                        phenotype = outcome.118$ECC.bin, batch = outcome.118$RNA.date,
                        sampleSum = sampleSum)
      tmp <- tmp[!tmp$id %in% c(352, 420),]
print(dim(tmp))
      LB(tmp)  # 3 (zero, nonzero, global) x 2 (coef, pval) ... so c(3,6) to be extracted
    })
    gene.BEZI %<>% t
    rownames(gene.BEZI) = gene.index.k
    gene.BEZI = list(coef = gene.BEZI[,1:3],
                     pval = gene.BEZI[,4:6],
                     qval = gene.BEZI[,6] %>% p.adjust(method="BH"))
    gene.BEZI$sigFeature = gene.list[gene.BEZI$qval < 0.05]
    
    saveRDS(gene.BEZI, paste0("output/R03.01.BEZI.R.116.bin.", k,".rds"))
    
  } else {
   # gene.BEZI <- readRDS("output/R03.01.BEZI.R.116.bin.rds")  
    gene.BEZI <- lapply(1:11, function(k) readRDS(paste0("output/R03.01.BEZI.R.116.bin.",k,".rds"))) 
    gene.BEZI.com <- gene.BEZI[[1]]
    for (k in 2:11) {
      gene.BEZI.com$coef = rbind(gene.BEZI.com$coef, gene.BEZI[[k]]$ceof)
      gene.BEZI.com$pval = rbind(gene.BEZI.com$pval, gene.BEZI[[k]]$pval)
      gene.BEZI.com$qval = c(gene.BEZI.com$qval, gene.BEZI[[k]]$qval)
      gene.BEZI.com$sigFeature = c(gene.BEZI.com$sigFeature, gene.BEZI[[k]]$sigFeature)
    }
  }
gene.BEZI.com$qval = (gene.BEZI.com$pval[,3] %>% p.adjust(method="BH"))
gene.BEZI.com$sigFeature = gene.list[gene.BEZI.com$qval<.05]  
gene.BEZI <- gene.BEZI.com; rm(gene.BEZI.com)

### 3. MAST
  if (method == 3) {
print("MAST")
    tmp <- data.frame(y = gene.marginal.RPK.RNA[gene.index,-1] %>% t, id = outcome.118$id, 
                      phenotype = outcome.118$ECC.bin, batch = outcome.118$RNA.date,
                      sampleSum = sampleSum)
    tmp <- tmp[!tmp$id %in% c(352, 420),]
    gene.MAST = MAST(tmp)  # 3 (zero, nonzero, global) x 2 (coef, pval) ... so c(3,6) to be extracted
    saveRDS(gene.MAST, "output/R03.01.MAST.R.116.bin.rds") # to be safe
    gene.MAST = list(t(gene.MAST[[1]]), t(gene.MAST[[2]]))
    rownames(gene.MAST[[1]]) <- rownames(gene.MAST[[2]]) <- gene.list
    names(gene.MAST) = c("coef", "pval")
    gene.MAST$pval[which(is.na(gene.MAST$pval[,3])),3] = gene.MAST$pval[which(is.na(gene.MAST$pval[,3])),2]
    gene.MAST$qval = gene.MAST[[2]][,3] %>% p.adjust(method="BH")
    gene.MAST$sigFeature = gene.list[gene.MAST$qval < 0.05]
    # add gene info
    saveRDS(gene.MAST, "output/R03.01.MAST.R.116.bin.rds")
    
  } else {
   # gene.MAST <- readRDS("output/R03.01.MAST.R.116.bin.rds")  
   
  }
  
  
### 4. KW
  if (0) {
    j <<- 1; cat("j=")
    gene.KW <- sapply(gene.index, function(i) {
      cat(j,".." ); j <<- j+1
      tmp <- data.frame(y = gene.marginal.RPK.RNA[i,-1] %>% as.numeric, id = outcome.118$id, 
                        phenotype = outcome.118$ECC.bin, batch = outcome.118$RNA.date)
      tmp <- tmp[!tmp$id %in% c(352, 420),]
      KW(tmp)  # 3 (zero, nonzero, global) x 2 (coef, pval) ... so c(3,6) to be extracted
    })
    gene.KW %<>% t
    rownames(gene.KW) = gene.index
    gene.KW = list(coef = gene.KW[,1],
                     pval = gene.KW[,2],
                     qval = gene.KW[,2] %>% p.adjust(method="BH"))
    gene.KW$sigFeature = gene.list[gene.KW$qval < 0.05]
    saveRDS(gene.KW, "output/R03.01.KW.R.116.bin.rds")
    
  } else {
    #gene.KW <- readRDS("output/R03.01.KW.R.116.bin.rds")  
  }
  
### 5.Wagner
  if (method == 5) {
print("Wagner")
    j <<- 1; cat("j=")
    gene.Wg <- sapply(gene.index, function(i) {
      cat(j,".." ); j <<- j+1
      tmp <- data.frame(y = gene.marginal.RPK.RNA[i,-1] %>% as.numeric, id = outcome.118$id, 
                        phenotype = outcome.118$ECC.bin, batch = outcome.118$RNA.date)
      tmp <- tmp[!tmp$id %in% c(352, 420),]
      Wagner(tmp)  # 3 (zero, nonzero, global) x 2 (coef, pval) ... so c(3,6) to be extracted
    })
    gene.Wg %<>% t
    rownames(gene.Wg) = gene.index
    colnames(gene.Wg) = rep(c("LB.nonz","LB.zero", "LB.glob"),2)
    
    gene.Wg = list(coef = gene.Wg[,1:3],
                   pval = gene.Wg[,4:6],
                   qval = gene.Wg[,6] %>% p.adjust(method="BH"))
    gene.Wg$sigFeature = gene.list[gene.Wg$qval < 0.05]
    saveRDS(gene.Wg, "output/R03.01.Wg.R.116.bin.rds")
    
  } else {
    gene.Wg <- readRDS("output/R03.01.Wg.R.116.bin.rds")  
  }
  



### Summary of significance results of each model
### orders of sig gene/path
gene.sig.list = data.frame(id = gene.list, LN = 0, LB = 0, MAST = 0, KW = 0, Wg = 0)



# 1. binary coding
gene.list %in% gene.lm$sigFeature %>% ifelse(1,0) -> gene.sig.list$LN
gene.list %in% gene.BEZI$sigFeature %>% ifelse(1,0) -> gene.sig.list$LB
gene.list %in% gene.MAST$sigFeature %>% ifelse(1,0) -> gene.sig.list$MAST
gene.list %in% gene.KW$sigFeature %>% ifelse(1,0) -> gene.sig.list$KW
gene.list %in% gene.Wg$sigFeature %>% ifelse(1,0) -> gene.sig.list$Wg

## summary  
apply(gene.sig.list[,-1], 2, sum)  # 222 30 179 0 1


## Venn diagram for genes
png("../Document/plot118/P0509-venn-R.bin.png", width=300, height=300)
vennCounts(gene.sig.list[,c(2:3)]) -> tmp
vennDiagram(tmp, include = "both", 
            names = c("LL","LB", "MAST"), 
            cex = 1, counts.col = "red", main = "# significant genes for each model")
text(-1, -2, "(LM = NP = 0)", col="blue", cex=0.7)
dev.off()


### some results
  # three common genes
  which(apply(gene.sig.list[,-1], 1, sum)>=3) -> a
  gene.sig.list[a,1] -> a1
  # 30 LB sig genes
  gene.sig.list[gene.sig.list$LB>0,1] -> a2
  
  pval.a1 = rep(1, dim(gene.marginal.RPK.RNA)[1])
  sapply(a1, function(x) grep(x, gene.marginal.RPK.RNA[,1])) -> a1.index
  pval.a1[a1.index] = 
  sapply(a2, function(x) grep(x, gene.marginal.RPK.RNA[,1])) -> a2.index
  
  source("../code-function/F00.00.boxScatterPlot.R")
  box.scatter("Document/plot-LogisticBeta.pdf", p.val = gene.BEZI$qval, 
              outcome=outcome.118 %>% mutate(ECC= ifelse(ECC==0, 0, 1)),
              xtext = c("caries-free", "ECC"),
              data = gene.marginal.RPK.RNA[gene.BEZI$qval %>% names() %>% as.numeric,])
  
  
  