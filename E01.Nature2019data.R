### 0. raw data read-in
  b <- read.delim("Nature2019data/ecs_relabRNA.tsv", sep = "\t")
  b2 <- read.delim("Nature2019data/hmp2_metadata.csv", sep = ",")

  ## leaving relevant data only. (first visit, RNA data, ...)
  b.meta <- 
    b2 %>% 
    dplyr::filter(data_type == "metatranscriptomics") %>% 
    dplyr::filter(External.ID %in% names(b)) %>%   # only those subjects in the RNA data remain.
    transmute(Participant.ID, External.ID, week_num, site_name,
              reads_human, Age.at.diagnosis, diagnosis, 
              race, Hispanic.or.Latino.Origin, sex) %>% 
    group_by(Participant.ID) %>% 
    mutate(visit = rank(week_num, ties.method = "first")) %>%  # there are ties for the first visits.
    dplyr::filter(visit == 1) %>% 
    dplyr::select(-visit) %>% 
    ungroup

    dim(b.meta) # 104 subjects x 8 variables
    
    # rearrangement of the RNA data
    sample.list = b.meta$External.ID
    sample.index.in.b = sapply(sample.list, function(x) which(x == names(b)))

### 1. full data (gene total and gene-species)
  dat.RNA = 
    list(otu = b[, sample.index.in.b],
         meta = b.meta,
         taxa = data.frame(gene.bact = b$X..Gene.Family, gene = NA, gene.name = NA, bact = NA))
  
### 2. Separating out the joint and the marginal data
  dat.RNA$taxa$gene = gsub("\\:.*" , "", dat.RNA$taxa$gene.bact)
  dat.RNA$taxa$gene.name = gsub(".*: (.*)|.*" , "\\1", dat.RNA$taxa$gene.bact)
  dat.RNA$taxa$bact = ifelse(grepl("\\|", dat.RNA$taxa$gene.bact), 
                             gsub(".*\\|", "", dat.RNA$taxa$gene.bact) %>% gsub(".*.s__" , "", .), 
                             "(TOTAL)")
  
  dat.RNA.marginal <- dat.RNA.joint <- dat.RNA
  index.joint = dat.RNA$taxa$bact != "(TOTAL)"
  dat.RNA.joint$otu     <- dat.RNA.joint$otu [index.joint, ]
  dat.RNA.joint$taxa    <- dat.RNA.joint$taxa[index.joint, ]
  dat.RNA.marginal$otu  <- dat.RNA.marginal$otu [!index.joint, ]
  dat.RNA.marginal$taxa <- dat.RNA.marginal$taxa[!index.joint, ]

  saveRDS(dat.RNA, "Nature2019data/data.ecs_relab.geneProp.full.RNA.rds")
  saveRDS(dat.RNA.marginal, "Nature2019data/data.ecs_relab.geneProp.marginal.RNA.rds")
  saveRDS(dat.RNA.joint, "Nature2019data/data.ecs_relab.geneProp.joint.RNA.rds")


### 3. marginal bacteria data
  bact.list = dat.RNA.joint$taxa$bact %>% unique %>% sort
  bact.list = c(bact.list[bact.list != "unclassified"], "unclassified")
  bact.index = lapply(bact.list, function(x) which(x == dat.RNA.joint$taxa$bact))
  
  # skeleton
  dat.RNA.bact = list(otu = matrix(NA, length(bact.list), dim(b.meta)[1]),
                      meta = b.meta,
                      taxa = data.frame(bacteria = bact.list))
  
  # filling in the marginalized expressions per species
  for (i in seq_along(bact.index)) {
    dat.RNA.bact$otu[i, ] = apply(dat.RNA.joint$otu[bact.index[[i]], , drop = FALSE], 2, sum)
  }
  saveRDS(dat.RNA.bact, "Nature2019data/data.ecs_relab.bactProp.marginal.RNA.rds")
  