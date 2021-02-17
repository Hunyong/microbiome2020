source("FE01.01.data_reader.R")
### 0. Getting the subjects list
dat.marginal = readRDS(sprintf("Nature2019data/data.geneRPK.marginal.RNA.rds"))
fns = list.files(path = "Nature2019data/MGX_unzip", full.names = T)  # unnecessary to remove "./" ?
fns.subset = sapply(dat.marginal$meta$External.ID, function(x) grep(paste0(x, "_[^TR]"), fns, value = TRUE))
for (i in length(fns.subset):1) {
  if (length(fns.subset[[i]]) == 0) fns.subset[[i]] = NULL
}
fns.subset = unlist(fns.subset)
length(fns) #1338
length(fns.subset) # 104

IBD.mgx <- jAnalyze(files = fns.subset, id = names(fns.subset), type="genefamilies")

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

# make an array of (overlapped geneRPK) x (sample) x (DNA+RNA)
# geneRPK.full.DNA.wide.sort <- readRDS("../Data-processed/tmp.data.geneRPK.full.DNA.wide.sort.ZOE2.rds")
IBD.full.DRNA <- array(NA, dim = c(dim(IBD.mgx.wide)[1],
                                   dim(IBD.mgx.wide)[2] - 3, 2), 
                       dimnames = list(IBD.mgx.wide$gene.bact, 
                                       dat.marginal$meta$External.ID, 
                                       c("DNA", "RNA")))
IBD.full.DRNA[, , 1] <- IBD.mgx.wide[, -(1:3)] %>% as.matrix
IBD.full.DRNA[, , 2] <- NA
IBD.full.DRNA <- list(otu = IBD.full.DRNA,
                      taxa = IBD.mgx.wide[, 1:3],
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

}

