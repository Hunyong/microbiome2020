### 0. raw data read-in
DRNA = "RNA"
  
  library(dplyr)
  nm = switch(DRNA, 
              DNA = "Nature2019data/ecs_relabDNA.tsv", 
              RNA = "Nature2019data/ecs_relabRNA.tsv") %>% print
  nm.genomics = switch(DRNA,
                       DNA = "metatranscriptomics",
                       RNA = "metagenomics") %>% print
  nm.meta = "Nature2019data/hmp2_metadata.csv"
  
  b <- read.delim(nm, sep = "\t")
  b2 <- read.delim(nm.meta, sep = ",")
  
  ## leaving relevant data only. (first visit, RNA data, ...)
  b.meta <- 
    b2 %>% 
    dplyr::filter(data_type == nm.genomics) %>% 
    dplyr::filter(External.ID %in% names(b)) %>%   # only those subjects in the RNA data remain.
    transmute(Participant.ID, External.ID, week_num, site_name,
              reads_human, Age.at.diagnosis, diagnosis, 
              race, Hispanic.or.Latino.Origin, sex) %>% 
    group_by(Participant.ID) %>% 
    mutate(visit = rank(week_num, ties.method = "first")) %>%  # there are ties for the first visits.
    dplyr::filter(visit == 1) %>% 
    dplyr::select(-visit) %>% 
    ungroup

    dim(b.meta) # 104 subjects x 8 variables (for RNA)
                # 109 subjects x 10 variables (for DNA)
    
    # rearrangement of the RNA data
    sample.list = b.meta$External.ID
    sample.index.in.b = sapply(sample.list, function(x) which(x == names(b)))

### 1. full data (gene total and gene-species)
  dat = 
    list(otu = b[, sample.index.in.b],
         meta = b.meta,
         taxa = data.frame(gene.bact = b$X..Gene.Family, gene = NA, gene.name = NA, bact = NA))
  
### 2. Separating out the joint and the marginal data
  dat$taxa$gene = gsub("\\:.*" , "", dat$taxa$gene.bact)
  dat$taxa$gene.name = gsub(".*: (.*)\\|.*" , "\\1", dat$taxa$gene.bact)
  dat$taxa$bact = ifelse(grepl("\\|", dat$taxa$gene.bact), 
                             gsub(".*\\|", "", dat$taxa$gene.bact) %>% gsub(".*.s__" , "", .), 
                             "(TOTAL)")
  
  dat.marginal <- dat.joint <- dat
  index.joint = dat$taxa$bact != "(TOTAL)"
  dat.joint$otu     <- dat.joint$otu [index.joint, ]
  dat.joint$taxa    <- dat.joint$taxa[index.joint, ]
  dat.marginal$otu  <- dat.marginal$otu [!index.joint, ]
  dat.marginal$taxa <- dat.marginal$taxa[!index.joint, ]

  saveRDS(dat, sprintf("Nature2019data/data.ecs_relab.geneProp.full.%s.rds", DRNA))
  saveRDS(dat.marginal, sprintf("Nature2019data/data.ecs_relab.geneProp.marginal.%s.rds", DRNA))
  saveRDS(dat.joint, sprintf("Nature2019data/data.ecs_relab.geneProp.joint.%s.rds", DRNA))


### 3. marginal bacteria data
  bact.list = dat.joint$taxa$bact %>% unique %>% sort
  bact.list = c(bact.list[bact.list != "unclassified"], "unclassified")
  bact.index = lapply(bact.list, function(x) which(x == dat.joint$taxa$bact))
  
  # skeleton
  dat.bact = list(otu = matrix(NA, length(bact.list), dim(b.meta)[1]),
                      meta = b.meta,
                      taxa = data.frame(bacteria = bact.list))
  
  # filling in the marginalized expressions per species
  for (i in seq_along(bact.index)) {
    dat.bact$otu[i, ] = apply(dat.joint$otu[bact.index[[i]], , drop = FALSE], 2, sum)
  }
  saveRDS(dat.bact, sprintf("Nature2019data/data.ecs_relab.bactProp.marginal.%s.rds", DRNA))
  