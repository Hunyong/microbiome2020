####################################################################################################
########## Functions - joint analysis
########## functions jointly looking into the data
####################################################################################################

### Read me: 
####################################################################################################
### 0. Library and working directory 
library(reshape2)
####################################################################################################
### 1. genefamilies
# 
jAnalyze <- function(type, direction = "long", ...) {
  # ... arguments:  study = 170421, id, type="genefamilies", info.file=humann2 
  # long: Yijk
  # hybrid: Yjk x i (gene-bacteria combination x subject)
  # array: j x k x i (3 dimension)
  data <- reader(..., type = type)
  if (class(data) == "data.frame") { data <- list(data) } else if (class(data) != "list") {
    stop ("The data is neither a data.frame nor a list.")
  }
  id = names(data)
  
  # for (i in 1:length(data)) {
  #   print(i)
  #   if (i == 15) print(data[[i]])
  #   .jAnalyze.ind(data[[i]], type = type, direction = direction)
  # }
  # joint-analysis by individual
  data <- lapply(data, .jAnalyze.ind, type = type, direction = direction)
  
  # Stack up and add id variable
  id <- rep(names(data), sapply(data, nrow))
  data <- as.data.frame(do.call("rbind", data))
  data$id <- id
  if (direction != "long") {
    stop("wide form: TBD")
    # long to wide
  }
  
  # cleansing
  rownames(data) <- NULL                   # remove rownames
  data[,1] <- gsub("\\:.*","", data[,1])  # remove description from gene names
  data[is.na(data)] <- 0                  # replace NA with 0
  return(data)
}

# marginalized analysis - 1. individual integration
# type: genefamilies, metaphlan, pathway, # margin: gene, pathway, bacteria
.jAnalyze.ind <- function(data, type = "genefamilies", direction = "long") {
  # data: first column = gene (or pathway) + bacteria, second = count (or percentage)
  if (dim(data)[1] == 0) return(data.frame())
  if (direction != "long") {
    stop("wide form: TBD")
  }
  
  if (grepl(type, "metaphlan2.tsv")) {
    type = "methaphlan"
    name = c("bacteria", "category")
    #unit = "percent"
  } else if (grepl(type, "pathabundance.tsv.pathcoverage.tsv")){
    type = "path"
    name = c("path", "bacteria")
    #unit = ifelse (grepl(type, "pathabundance"), "count" , "coverage")
  } else if (type == "bracken") {
    type = "bracken"
    name = c("bacteria", "category")
    #unit = ifelse (grepl("cpm", type), "CPM" , "RPK")
  } else {
    type = "gene"
    name = c("gene", "bacteria")
    #unit = ifelse (grepl("cpm", type), "CPM" , "RPK")
  }
  
  name[3] <- names(data)[2] # adding units
  data <- data.frame(a = data[,1], b = NA, d = data[,2])
  names(data) <- name
  data[,2] <- ifelse(grepl("\\|", data[,1]), gsub(".*\\|","",data[,1]), "(TOTAL)")
  data[,1] <- gsub("\\|.*", "", data[,1])
  
  # removed these lines (this causes inconsistency across subjects) on Oct 9, 2018
  #tmp <- table(data[,1])      # table of counts
  #tmp <- names(tmp)[tmp == 1] # if the counts is 1 (e.g. UNMAPPED), the (TOTAL) should remain, o/w to be removed.
  #data <- data[(data[,1] %in% tmp) | data[,2] != "(TOTAL)", ]
  
  return(data)
}

####################################################################################################
### 3. reader: reading files, given study, id, file type.
reader <- function (files, id, type = "genefamilies-cpm",
                    colClasses = c("character", "numeric"), 
                    col = if (type == "bracken") c(1,6) else "all", 
                    complete = !(type == "bracken")) {
  if (length(files) > 1) {
    data <- lapply (files, .read.delim.col, header=(type == "bracken"), 
                    skip = ifelse(type == "bracken", 0, 1), 
                    colClasses = colClasses, col = col, type = type, complete = complete, zipped = FALSE)
    names(data) <- id
  } else {
    data <- .read.delim.col(files, header=(type == "bracken"), 
                            skip = ifelse(type == "bracken", 0, 1),
                            colClasses = colClasses, col = col, type = type, complete = complete, zipped = FALSE)
  }
  return(data)
}
# wrapper function of read.delim2 (extracting specific columns only)
.read.delim.col <- function (fn, ..., col, type = NULL, complete = TRUE, zipped = FALSE) {
  if (zipped) {
    zip.file = gsub("\\@.*", "", fn)
    member.file = gsub(".*\\@", "", fn)
    fn <- unz(zip.file, member.file)
  }
  data <- read.csv(fn, ..., stringsAsFactors=FALSE, sep = "\t")
  # print(data)    
  if (complete) {  # this module was added on Oct 9, 2018 to fill the gap between total and sum(items).
    require(magrittr); require(dplyr)
    names(data)[1:2] <- c("V1", "V2")
    data$index = seq_len(nrow(data))
    tmp.marginal <- data %>% dplyr::filter(!grepl("\\|", V1))
    tmp.joint    <- data %>% 
      dplyr::filter(grepl("\\|", V1)) %>% 
      dplyr::mutate(V1 = gsub("\\|.*", "", V1)) %>% 
      dplyr::select(-index) %>%
      group_by(V1) %>% 
      summarize(V2 = sum(V2)) %>% 
      ungroup
    tmp <- dplyr::left_join(tmp.marginal, tmp.joint, by = "V1")
    tmp[is.na(tmp$V2.y), "V2.y"] <- 0  # forcing NA into 0
    tmp %<>% 
      transmute(V1 = paste0(V1, "|(GAP)"), V2 = V2.x - V2.y, index = index + .5) %>%
      filter(V2 >= 1)  #throw away small errors
    data <- rbind(data, tmp) %>% arrange(index) %>% dplyr::select(-index)
  }
  
  if (col[1] == "all") {col = 1:dim(data)[2]} else {col = which((1:dim(data)[2]) %in% col)}
  if (grepl(type, "bracken")) {
    names(data) <- gsub("new\\_est\\_reads", "reads", names(data))
  } else {
    if (!is.null(type)) {
      if (grepl(type, "genefamilies.")) {name = c("gene.bacteria", "RPK")
      } else if (grepl(type, "genefamilies-cpm.")) {name = c("gene.bacteria", "CPM")
      } else if (grepl(type, "pathabundance.")) {name = c("path.bacteria", "count")
      } else if (grepl(type, "pathabundance.")) {name = c("path.bacteria", "proportion")
      } else if (grepl(type, "metaphlan2.")) {name = c("bacteria", "percentage")
      }
    }
    names(data) <- name
  }
  data <- data[,col]
  data
}
