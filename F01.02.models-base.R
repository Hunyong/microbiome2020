###### F01.02.screen.gene.R  ##################################################################
###### Objective: filtering unexpressed genes/paths
###### Warning: ...
############################################################################################

screen.gene <- function(data, avg.detect = 2, gene.col = NULL) {
  # data: each column = sample, each row = genes
  # gene.col: the column where gene names are provided. If null, ignored.
  gene.name = data[,gene.col]
  if (!is.null(gene.col)) data = data[,-gene.col]
  n.genes = dim(data)[1]
  n.sample = dim(data)[2]
  
  useF = which(rowSums(data) > avg.detect * n.sample)
  print(paste(length(useF), "out of ", n.genes, " will be used."))
  list(sample.size = n.sample, 
       stat = c(use = length(useF), out.of = n.genes, percentage = round(length(useF)/n.genes,2)*100),
       feature = list(index = useF, names = if (is.null(gene.col)) NA else gene.name[useF]))
}
if (FALSE) {
  gene.marginal.RPK.RNA <- readRDS("output118/gene.marginal.RPK.RNA.118.rds")
  path.marginal.RNA <- readRDS("output118/path.marginal.RNA.118.rds")
  gene.list = screen.gene(gene.marginal.RPK.RNA, gene.col=1)$feature$names
  path.list = screen.gene(path.marginal.RNA, gene.col=1)$feature$names
}
