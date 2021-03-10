library(dplyr); library(tidyr)
source("FE01.01.data_reader.R")
### 0. Getting the subjects list
dat.marginal = readRDS(sprintf("Nature2019data/data.geneRPK.marginal.RNA.rds"))

### DNA
  fns = list.files(path = "Nature2019data/MGX_unzip", full.names = T)  # unnecessary to remove "./" ?
  fns.subset = sapply(dat.marginal$meta$External.ID, function(x) grep(paste0(x, "_[^TR]"), fns, value = TRUE))
  for (i in length(fns.subset):1) {
    if (length(fns.subset[[i]]) == 0) fns.subset[[i]] = NULL
  }
  fns.subset = unlist(fns.subset)
  length(fns) #1338
  length(fns.subset) # 104
  
  IBD.mgx <- jAnalyze(files = fns.subset, id = names(fns.subset), type="genefamilies")
  IBD.mgx$RPK %<>% as.character %>% as.numeric()
  
  IBD.mgx.wide = 
    IBD.mgx %>% 
    transmute(gene.bact = paste0(gene, " ", bacteria), gene, bacteria, RPK, id) %>%
    spread(key = id, value = RPK, fill = 0)
  
  # filling in missing subjects: Note the DNA data have more subjects
  # tmp <- names(IBD.mgx.wide)
  tmp <- dat.marginal$meta$External.ID[!dat.marginal$meta$External.ID %in% names(IBD.mgx.wide)]
  IBD.mgx.wide = cbind(IBD.mgx.wide, matrix(NA, nrow = dim(IBD.mgx.wide)[1], ncol = length(tmp), 
                                            dimnames = list(NULL, tmp)))
  
  # ordering the subjects - DNA
  tmp <- sapply(dat.marginal$meta$External.ID, function(x) which(names(IBD.mgx.wide)[-c(1:3)] == x))
  IBD.mgx.wide <- cbind(IBD.mgx.wide[, 1:3], IBD.mgx.wide[, tmp+3])
  
### RNA
  fns = list.files(path = "Nature2019data/MTX_unzip", full.names = T)  # unnecessary to remove "./" ?
  fns.subset = sapply(dat.marginal$meta$External.ID, function(x) grep(paste0(x, "_[^TR]"), fns, value = TRUE))
  for (i in length(fns.subset):1) {
    if (length(fns.subset[[i]]) == 0) fns.subset[[i]] = NULL
  }
  fns.subset = unlist(fns.subset)
  length(fns) #1338
  length(fns.subset) # 103
  
  IBD.mtx <- jAnalyze(files = fns.subset, id = names(fns.subset), type="genefamilies")
  IBD.mtx$RPK %<>% as.character %>% as.numeric()
  
  IBD.mtx.wide = 
    IBD.mtx %>% 
    transmute(gene.bact = paste0(gene, " ", bacteria), gene, bacteria, RPK, id) %>%
    spread(key = id, value = RPK, fill = 0)
  
  # filling in missing subjects: Note the DNA data have more subjects
  # tmp <- names(IBD.mgx.wide)
  tmp <- dat.marginal$meta$External.ID[!dat.marginal$meta$External.ID %in% names(IBD.mtx.wide)]
  IBD.mtx.wide = cbind(IBD.mtx.wide, matrix(NA, nrow = dim(IBD.mtx.wide)[1], ncol = length(tmp), 
                                            dimnames = list(NULL, tmp)))
  
  # ordering the subjects - DNA
  tmp <- sapply(dat.marginal$meta$External.ID, function(x) which(names(IBD.mtx.wide)[-c(1:3)] == x))
  IBD.mtx.wide <- cbind(IBD.mtx.wide[, 1:3], IBD.mtx.wide[, tmp+3])
  
  
### Matching the genes
  # how many overlapped?
  # geneRPK.full.DNA.wide <- readRDS(geneRPK.DNA.wide.fn)
  # geneRPK.full.RNA.wide <- readRDS(geneRPK.RNA.wide.fn)
  geneRPK.bact.RNAinDNA = which (IBD.mtx.wide[,1] %in% IBD.mgx.wide[,1])
  geneRPK.bact.names = IBD.mtx.wide[geneRPK.bact.RNAinDNA, 1]
  geneRPK.bact.DNAinRNA = which (IBD.mgx.wide[,1] %in% geneRPK.bact.names)
  geneRPK.bact.names.D = IBD.mgx.wide[geneRPK.bact.DNAinRNA, 1]
  
  geneRPK.bact.RNAonly = which (!IBD.mtx.wide[,1] %in% IBD.mgx.wide[,1])
  geneRPK.bact.RNAonly.names = IBD.mtx.wide[geneRPK.bact.RNAonly, 1]
  geneRPK.bact.DNAonly = which (!IBD.mgx.wide[,1] %in% geneRPK.bact.names)
  geneRPK.bact.DNAonly.names = IBD.mgx.wide[geneRPK.bact.DNAonly, 1]
  
  ###
  length(geneRPK.bact.names)           # 1138479    # 694,672 (wv)
  length(geneRPK.bact.names.D)         # 1138479    # 694,672 (wv)
  length(geneRPK.bact.RNAonly.names)   #   92616    #  72,868 (wv)
  length(geneRPK.bact.DNAonly.names)   # 1400655    # 170,699 (wv)
  
  # 1.0 data manipulation
  IBD.mtx.wide.sort <- IBD.mtx.wide[geneRPK.bact.RNAinDNA,]
  IBD.mtx.wide.sort %<>%  arrange(gene.bact)
  IBD.mgx.wide.sort <- IBD.mgx.wide[geneRPK.bact.DNAinRNA,]
  IBD.mgx.wide.sort %<>%  arrange(gene.bact)
  geneRPK.bact.names.sort <- IBD.mtx.wide.sort[,1:3]
  # make sure the subjects and geneRPKs are in the same order
  all(names(IBD.mtx.wide.sort)[-1] == names(IBD.mgx.wide.sort)[-1])
  all(IBD.mtx.wide.sort[,1] == IBD.mgx.wide.sort[,1])
  all(names(IBD.mtx.wide.sort)[-(1:3)] == dat.marginal$meta$External.ID)
  
  tmp <- IBD.mgx.wide[geneRPK.bact.DNAonly, ]
  tmp[,-(1:3)] <- 0
  IBD.mtx.wide.sort <- rbind(IBD.mtx.wide.sort,                    # shared
                             IBD.mtx.wide[geneRPK.bact.RNAonly, ], # RNA only
                             tmp)                                  # DNA only (zeros)
  
  tmp <- IBD.mtx.wide[geneRPK.bact.RNAonly, ]
  tmp[,-(1:3)] <- 0
  IBD.mgx.wide.sort <- rbind(IBD.mgx.wide.sort,                 # shared
                             tmp,                                     # RNA only (zeros)
                             IBD.mgx.wide[geneRPK.bact.DNAonly,])  # DNA only
  # sorting once more, after putting them altogether
  index <- order(IBD.mtx.wide.sort[,1])
  IBD.mtx.wide.sort <- IBD.mtx.wide.sort[index, ]
  IBD.mgx.wide.sort <- IBD.mgx.wide.sort[index, ]
  
  
  
  


# make an array of (overlapped geneRPK) x (sample) x (DNA+RNA)
# geneRPK.full.DNA.wide.sort <- readRDS("../Data-processed/tmp.data.geneRPK.full.DNA.wide.sort.ZOE2.rds")
IBD.full.DRNA <- array(NA, dim = c(dim(IBD.mtx.wide.sort)[1],
                                   dim(IBD.mtx.wide.sort)[2] - 3, 2), 
                       dimnames = list(IBD.mtx.wide.sort$gene.bact, 
                                       dat.marginal$meta$External.ID, 
                                       c("DNA", "RNA")))
IBD.full.DRNA[, , 1] <- IBD.mgx.wide.sort[, -(1:3)] %>% as.matrix
IBD.full.DRNA[, , 2] <- IBD.mtx.wide.sort[, -(1:3)] %>% as.matrix
IBD.full.DRNA <- list(otu = IBD.full.DRNA,
                      taxa = IBD.mtx.wide.sort[, 1:3],
                      meta = dat.marginal$meta)
# rm(geneRPK.full.RNA.wide.sort); gc()

# splitting marginal & joint data from full data
IBD.full.DRNA %$% otu[taxa$bacteria == "(TOTAL)",,] -> IBD.marginal.DRNA
IBD.full.DRNA %$% otu[taxa$bacteria != "(TOTAL)" & taxa$bacteria != "(GAP)",,] -> IBD.joint.DRNA
IBD.marginal.DRNA <- list(otu = IBD.marginal.DRNA,
                          taxa = IBD.full.DRNA$taxa[IBD.full.DRNA$taxa$bacteria == "(TOTAL)", ],
                          meta = dat.marginal$meta)
IBD.joint.DRNA <- list(otu = IBD.joint.DRNA,
                       taxa = IBD.full.DRNA$taxa[IBD.full.DRNA$taxa$bacteria != "(TOTAL)"& 
                                                   IBD.full.DRNA$taxa$bacteria != "(GAP)", ],
                       meta = dat.marginal$meta)
rownames(IBD.marginal.DRNA$taxa) <- rownames(IBD.joint.DRNA$taxa) <- NULL
saveRDS(IBD.full.DRNA, "Nature2019data/data.geneRPK.full.DRNA.IBD.rds")
saveRDS(IBD.marginal.DRNA, "Nature2019data/data.geneRPK.marginal.DRNA.IBD.rds")
saveRDS(IBD.joint.DRNA, "Nature2019data/data.geneRPK.joint.DRNA.IBD.rds")

# zero proportion
zp.D = IBD.marginal.DRNA$otu[,,1] %>% apply(1, function(x) mean(x == 0, na.rm = TRUE))
zp.R = IBD.marginal.DRNA$otu[,,2] %>% apply(1, function(x) mean(x == 0, na.rm = TRUE))
hist(zp.D[zp.D > 0])
hist(zp.R[zp.R > 0])
mean(zp.D) #87.8%; 
mean(zp.R) #96.3%
