library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
library(latex2exp)
source("F00.00.generic.R")
source("F01.01.base.R")
zero.prob <- function (vec) {mean(vec == 0, na.rm = TRUE)}

zoe = 2  # zoe = "IBD"
type = "gene"
DRNA = "RNA"
nrm = "tpm5" #"rpk" "asin"
n.samp = "full" # or n.samp = 300

for (DRNA in c("RNA")) {
  for (zoe in (1:2))  {
    zoe.nm = if (zoe %in% 1:2) paste0("_zoe", zoe) else "_NEWDATA"
    zoe.nm2 = if (zoe %in% 1:2) paste0("ZOE", zoe) else "NEWDATA"
    
    for (type in c("genebact", "bact", "gene")) {
      print(type)
      type.full = switch(type, gene = "geneRPK.marginal", 
                         genebact = "geneRPK.joint", bact = "bactRPK.marginal")
      DR.no = switch(DRNA, DNA = 1, RNA = 2)
      DRNA.name = ifelse(DRNA == "DNA", "_DNA", "")
    
      ### 0.2 Data
      if (zoe %in% 1:2) {
        fn   = sprintf("../Data-processed/data.%s.DRNA.%s.rds", type.full, zoe.nm2)
        data = readRDS(fn)
        excluded.subject <- data$meta$id %in% c(352, 420, 10083, 12623, 11259, 11790)
        DataMeta = data$meta[!excluded.subject,]
        batchGrp = switch(zoe, "1" = "170628", "2" = "180530")
        DataMeta <-
          DataMeta %>% 
          mutate(group = paste0(ifelse(cariesfree == 1, "H", "D"), ifelse(batch.RNA == batchGrp, 1, 2)),
                 group_batch =  ifelse(batch.RNA == batchGrp, 1, 2),
                 group_disease = ifelse(cariesfree == 1, "H", "D"))
      } else {
        fn   = sprintf("../MicrobiomePaper2020/Nature2019data/data.%s.DRNA.IBD.rds", type.full)
        data = readRDS(fn)
        DataMeta = 
          data$meta %>% 
          mutate(id = External.ID,
                 group_batch = ifelse(site_name %in% c("Cedars-Sinai", "MGH"), 1, 2),
                 group_disease = ifelse(diagnosis %in% c("UC", "CD"), "D", "H"),
                 group = paste0(group_disease, group_batch))
      }
      
      RNA         = data$otu[,, DR.no]
      DataRPK116  = RNA[,colnames(RNA) %in% DataMeta$id]
      ST          = apply(DataRPK116, 2, sum, na.rm = TRUE)  # sample total
      mean(ST) # 21M, 0.99 for IBD
      n.genes     = dim(DataRPK116)[1]   # 12million for IBD
      scale1 = 10               # normalized relative abundance..  = average sampleSum / number of taxa.
      scale2 = scale1 * n.genes # actual scale (rel-abundance to compositional)
      
      DataTPM116  = t(t(DataRPK116)/ST) * scale2
      DataComp    = t(t(DataRPK116)/ST)  # for LB tests
      if (nrm == "asin") {
        DataTPM116 = asn(DataComp) * scale2
        DataComp = asn(DataComp)
      }
      
      # Choose dataset and estimators according to the nrm!
      Data = if (nrm == "rpk") DataRPK116 else DataTPM116
      rm(data, excluded.subject, RNA, DataRPK116, DataTPM116); gc()
      
      estr.ziln = function(yvec) {
        pp <- mean(yvec == 0)
        y.nz <- yvec[yvec != 0]
        
        c(mu = mean(y.nz),
          theta = var(y.nz) / mean(y.nz)^2,
          pi = pp)
        
      }
      estr.zinb = function(yvec) {
        ZINB.ML.time(yvec %>% round, parametrization="mtp")
      }
      
      for (model in c("ziln", "zinb")) {
        cat("model = ", model, "\n")
          estr = switch(model, ziln = estr.ziln, zinb = estr.zinb)
          
          # disease groups
          grp <- unique(DataMeta$group)
          
          
        ### marginal sample (not considering batches and disease groups)
          zero.proportion <- apply(Data, 1, zero.prob)
          genes.regular.index <- which(zero.proportion <= 0.9) # 29%
          
          
          
        ### Estimation 
        if(T){
          set.seed(1)
          if (is.null(n.samp) | n.samp == "full") {
            n.samp = "full"
            samp <- genes.regular.index
          } else {
            samp <- sample(genes.regular.index,
                           ifelse(type == "bact", length(genes.regular.index), n.samp))
          }
          n.samp.nm = paste0("samp", n.samp)
          
          ### 2. conditional estimates
          
          cond.est <- 
            lapply(grp, function(g) {
              cat("group ", g, "\n")
              sapply(samp, function(x) {
                # if (x %% 10) cat(x, " ")
                estr(Data[x, DataMeta$group == g] %>% na.omit)
              }) %>% t %>% 
                as.data.frame %>% 
                mutate(group = g, disease = substr(group, 1, 1), batch = substr(group, 2, 2),
                       index = samp)
            }) 
          if (any(sapply(cond.est, names)[1, ] == "V1")) {
            for (g in 1:length(cond.est)) {
              names(cond.est[[g]])[1:3] = c("mu", "theta", "pi")
            }
          }
          cond.est <- do.call(rbind, cond.est)
          names(cond.est) <- gsub("\\.\\(Intercept\\)", "", names(cond.est))
          saveRDS(cond.est, paste0("output/para_selection_est", DRNA.name, zoe.nm, "_", type, "_", model, "_", nrm, "_", n.samp.nm, ".rds"))
          
          ### >>> plots are outsourced to C01.02....plots.R
        }
      }
    }
  }
}
  