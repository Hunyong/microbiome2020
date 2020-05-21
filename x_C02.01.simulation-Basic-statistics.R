### library and settings  
source("C01.02.simulation.setup.R")
source("F02.01.simulation-analysis.R")
##install.packages("gamlss")
library(gamlss)
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
library(latex2exp)

model = "zinb"
model = "zig"
model = "ziln"
sig = 0.05
pert = 0.5
perturb = 0.5
n = 80; regular = 0.1 
n = 400; regular = 0.02  
cutoff = 0.9
test = c("LB", "LN", "MAST", "KW",  "KW-II", "DESeq2")
# e.g. for n=80, there should be at least 8 nonzeros for each gene.
#regular = 0.0 # e.g. for n=80, there should be at least 8 nonzeros for each gene.
if(TRUE){
  if(TRUE){ 
    save_path = paste0("output/")
    save_file <- paste0(save_path, "stat-n", n,"-pert",pert, "-reg", regular,"-", model, ".rds")
    fig_path = paste0("Document/plot/", Sys.Date(), "/")
    fig_file <- paste0(fig_path, "stat-n", n,"-pert",pert, "-reg", regular,"-", model)
    
    if (!dir.exists(save_path)) dir.create(file.path(save_path))
    if (!dir.exists(fig_path)) dir.create(fig_path)
    
    (parameter = switch(model, 
                        zinb = parameter4, 
                        zig = parameter5, 
                        ziln = parameter5))
    (delta = delta1); 
    (kappa = kappa1)
    
    length.delta = dim(delta)[1];
    length.kappa = dim(kappa)[1];
    length.parameter = dim(parameter)[1]
    length.whole = length.delta * length.kappa * length.parameter
  }
  ##
  
  #################### This is for loading Data ############
  # skeleton in case of no file
  tmp <- readRDS(paste0("output/result-n400-pert0.5-ziln-1.1.1.rds"))
  tmp <- readRDS(paste0("output/result-n400-pert0.5-zig-1.1.1.rds"))
  tmp[[1]][,] <- tmp[[2]][,] <- tmp[[3]][[1]] <- tmp[[3]][[2]] <- NA
  
  # generating empty slots
  a <- c(i=0, j=0, k=0, tmp$pval %>% apply(1, function(x) {mean(x<=0.05)}) )
  result.stat.na <- result.stat <- NA.proportion <-  
    matrix(a, nrow = length.whole, ncol = length(a), byrow = TRUE, 
           dimnames = list(NULL, names(a))) %>% as.data.frame
  rep.MAST <- length(tmp$MAST$pval)
  index.MAST <- tmp$pval %>% rownames() %>% {grep("MAST", .)} #5, 6, 7
  #index.twopart <- grep("(LB)|(MAST)|(Wg)", rownames(tmp$pval))
  
  ### STATISTICS ###
  row.index = 1
  for (i in 1:length.delta) {
    message("\ni = ", i, "\n")
    cat("(i, j, k) = ")
    for (j in 1:length.kappa) {
      for (k in 1:length.parameter) {
        cat("(", i, ", ", j , ", ", k, ") ")
        if (row.index > 0) {
          ## 1. Reading
          file.ijk <- paste0(save_path, "/result-n", n, "-pert", perturb, "-", model, "-", i, ".", j, ".", k, ".rds")
          if(file.exists(file.ijk)) {
            nofile <- FALSE
            rds <- readRDS(file.ijk)
          } else {
            nofile <- TRUE
            rds <- tmp
          }
          
          ## 2 MAST
          
          # 2.1 NA replacement   # This is not actually needed, but for consistency of data (btw first replicate and the rest)
          # step 1. getting NA addresses
          na.index = rds$pval[index.MAST[3],] %>% is.na %>% which
          # step 2. replacing with nonzero model values
          rds$pval[index.MAST[3], na.index] = rds$pval[index.MAST[2], na.index]
          # step 3. getting NA addresses again and replace with zero model values.
          na.index = rds$pval[index.MAST[3], ] %>% is.na %>% which
          rds$pval[index.MAST[3], na.index] = rds$pval[index.MAST[1], na.index]
          # leftovers
          # result[[i]][[j]][[k]]$pval[7,na.index] %>% length %>%"/"(n.gene) %T>% print
          
          # 2.2 MAST NA replacement for each replicate
          rds$MAST$pval <- 
            if (nofile) {NA} else {
              lapply(rds$MAST$pval, function(s) {
                # step 1. getting NA addresses
                na.index = s[3,] %>% is.na %>% which
                # step 2. replacing with nonzero model values
                s[3,na.index] = s[2,na.index]
                # step 3. getting NA addresses again and replace with zero model values.
                na.index = s[3,] %>% is.na %>% which
                s[3,na.index] = s[1,na.index]
                s
              })
            }
          
          
          # 3. Statistics
          # 3.1 statistics for all method (including first replicate of MAST, which is redundant)
          index.regular <- rds$nonzero.prop >= regular
          #rds$pval[index.twopart, !index.regular] <- NA    # removing irregular genes for two-part models only.
          rds$pval[, !index.regular] <- NA    # removing irregular genes.
          
          result.stat[row.index,] = 
            c(i=i, j=j, k=k, rds$pval %>% apply(1, function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, 
                                                                       mean(x<=sig, na.rm=TRUE))}))
          
          result.stat.na[row.index,] =
            c(i=i, j=j, k=k, rds$pval %>% apply(1, function(x) {ifelse(sum(is.na(x))/length(x)>0.9,
                                                                       NA,
                                                                       mean(ifelse(is.na(x), 1, x) <= sig, na.rm=TRUE))}) )
          
          stat_na_p_val <- function(data){
            if(sum(is.na(data))/length(data) > 0.9){
              
            }else{
              
            }
          }
          
          
          NA.proportion[row.index,] =
            c(i=i, j=j, k=k, rds$pval %>% apply(1, function(x) {mean(is.na(x))}))
          
          # 3.2 statistics for MAST replicates
          if (!nofile) {
            
            for (m in 1:10) {
              index.regular.m <- rds$MAST$nonzero.prop[[m]] >= regular
              rds$MAST$pval[[m]][, !index.regular.m] <- NA
            }
            
            result.stat[row.index, c("MAST.nonz", "MAST.zero", "MAST.glob")] <-
              rds$MAST$pval %>% 
              #sapply(function(s) apply(s, 1, function(x)  {mean(x<=sig, na.rm=TRUE)})) %>%  # 3 (tests) by 10 (replicates) matrix
              sapply(function(s) apply(s, 1, function(x) {ifelse(mean(is.na(x)) >= cutoff, NA, 
                                                                 mean(x<=sig, na.rm=TRUE))})) %>%
              apply(1, mean, na.rm = TRUE)  # vector of three
            
            
            
            result.stat.na[row.index, c("MAST.nonz", "MAST.zero", "MAST.glob")] <- #stat based on NA p-val = 1.
              rds$MAST$pval %>% 
              sapply(function(s) apply(s, 1, function(x)  {mean(ifelse(is.na(x), 1, x)<=sig, na.rm=TRUE)})) %>%
              apply(1, mean, na.rm = TRUE)
            
            NA.proportion[row.index, c("MAST.nonz", "MAST.zero", "MAST.glob")] <-
              rds$MAST$pval %>% 
              sapply(function(s) apply(s, 1, function(x)  {mean(is.na(x))})) %>%
              apply(1, mean, na.rm = TRUE)
          }
          gc()
          
        }  
        ## 4. updating row.index
        row.index = row.index + 1
        
      }
    }
  }
  
  ## 5. cosmetics 
  result.stat = result.stat %>% mutate(i, effect = delta[i, 4], j, batch = kappa[j,4])
  result.stat.na = result.stat.na %>% mutate(i, effect = delta[i, 4], j, batch = kappa[j,4])
  NA.proportion = NA.proportion %>% mutate(i, effect = delta[i, 4], j, batch = kappa[j,4])
  
  
  result <- list(base = data.frame(generative = model, significance = sig, n.gene = n.gene,
                                   num.replicates.MAST = rep.MAST), 
                 stat = result.stat,
                 stat.na = result.stat.na,
                 na.prop = NA.proportion)
  saveRDS(result, save_file)
}
if(FALSE){
  result.80 <- result
  result.400 <- result
  result_025 <- result
  result_000 <- result
  saveRDS(result.80, save_file)
  saveRDS(result.400, save_file)
  saveRDS(result_025, save_file)
  saveRDS(result_000, save_file)
}

###################### loading result part ###############
result <- readRDS(save_file)
result.stat.na <- result$stat.na
result.stat.na <- result$stat
# result.stat <- result$stat
# result.stat.na.80 <- result.80$stat.na
# result.stat.na.400 <- result.400$stat.na
# NA.proportion <- result$na.prop
result.stat.na$size = "n = 80"
result.stat.na$size = "n = 400"
##

if(FALSE){# This part is for diffferent pert
  result.stat.na.80$size = "n = 80"
  result.stat.na.400$size = "n = 400"
  result.stat.na = rbind(result.stat.na.80, result.stat.na.400)
  
  result.stat.na.025 <- result_025$stat.na
  result.stat.na.000 <- result_000$stat.na
  result.stat.na.025$size = "n = 80"
  result.stat.na.000$size = "n = 80"
  result.stat.na.80$perturbation = 0.5
  result.stat.na.025$perturbation = 0.25
  result.stat.na.000$perturbation = 0
  result.stat.na.pert = rbind(result.stat.na.80, result.stat.na.025, result.stat.na.000)
}


result.stat.na %>% 
  mutate ("LB" = LB.glob, "MAST" = MAST.glob, "KW-II" = Wg.glob) %>%
  dplyr::select(-c(4,5,6,8,9,10,12,13,14,16)) -> result.stat.na

result.stat.na.pert %>% 
  mutate ("LB" = LB.glob) %>%
  dplyr::select(-c(4,5,6,7,8,9,10,11,12,13,14,15)) -> result.stat.na.pert


##################### !!! plot for main article !!! ##############################

test = c("LB", "LN", "MAST", "KW",  "KW-II", "DESeq2")
#3.1. updated null plot
a <- pval.plot.null(result.stat.na, 
                    parameter_in_use = parameter,
                    i=1, test=test, 
                    k.index= c(1,3,4,6,7,9,10,12,19,21,22,24),
                    j.index= c(1,5,3), ylim = c(0,0.3) )
ggsave(paste0(fig_path, "main_", "_null_reduced_all_methods" ,".png"), a, width = 10, height=12)

#3.2. updated power plot - size 80
b_80 <- pval.plot.power(result.stat.na, 
                        parameter_in_use = parameter,
                        i = c(2,3,4,6,8),
                        k.index= c(1,3,4,6,7,9,10,12,19,21,22,24),
                        j.index= c(1,5,3), ylim = 0:1,
                        sample_size = 80)
ggsave(paste0(fig_path, "main_", "_power_reduced_all_methods_80" ,".png"), b_80, width = 20, height=12)

#3.2. updated power plot - size 400
b_400 <- pval.plot.power(result.stat.na, 
                         parameter_in_use = parameter,
                         i = c(2,3,4,6,8),
                         k.index= c(1,3,4,6,7,9,10,12,19,21,22,24),
                         j.index= c(1,5,3), ylim = 0:1,
                         sample_size = 400)
ggsave(paste0(fig_path, "main_", "_power_reduced_all_methods_400" ,".png"), b_400, width = 20, height=12)

#4. updated power plot - LB with pert = 0, 0.25, 0.5
c <- pval.plot.pert(result.stat.na.pert, 
                    parameter_in_use = parameter,
                    i = c(1,2,3,4,6,8),
                    k.index= c(1,3,4,6,7,9,10,12,19,21,22,24),
                    j.index= c(1,5,3), ylim = 0:1,
                    pert = c(0, 0.25, 0.5), sample_size = 80)
ggsave(paste0(fig_path, "main_", "_power_LB_perutrbation" ,".png"), c, width = 20, height=12)


#5.1.1 full power plot - ZILN size 80
d1_80 <- pval.plot.power.full(result.stat.na, 
                              parameter_in_use = parameter,
                              i = c(1,2,3,4,5,6,7,8,9,10),
                              k.index= 1:dim(parameter)[1],
                              j.index= 1:5, ylim = 0:1,
                              sample_size = 80)
ggsave(paste0(fig_path, "ZILN_full_80",".png"), d1_80, width = 20, height=21.6)

#5.1.2 full power plot - ZILN size 400
d1_400 <- pval.plot.power.full(result.stat.na, 
                               parameter_in_use = parameter,
                               i = c(1,2,3,4,5,6,7,8,9,10),
                               k.index= 1:dim(parameter)[1],
                               j.index= 1:5, ylim = 0:1,
                               sample_size = 400)
ggsave(paste0(fig_path, "ZILN_full_400",".png"), d1_400, width = 20, height=21.6)

#5.2.1 full power plot - ZINB size 80
d2_80 <- pval.plot.power.full(result.stat.na, 
                              parameter_in_use = parameter, test = test,
                              i = c(1,2,3,4,5,6,7,8,9,10),
                              k.index= 1:dim(parameter)[1],
                              j.index= 1:5, ylim = 0:1,
                              sample_size = 80)
ggsave(paste0(fig_path, "ZINB_full_80",".png"), d2_80, width = 20, height=21.6)

#5.2.2 full power plot - ZINB size 400
d2_400 <- pval.plot.power.full(result.stat.na, 
                               parameter_in_use = parameter, test = test,
                               i = c(1,2,3,4,5,6,7,8,9,10),
                               k.index= 1:dim(parameter)[1],
                               j.index= 1:5, ylim = 0:1,
                               sample_size = 400)
ggsave(paste0(fig_path, "ZINB_full_400",".png"), d2_400, width = 20, height=21.6)

#5.3.1 full power plot - ZIG size 80
d3_80 <- pval.plot.power.full(result.stat.na, 
                              parameter_in_use = parameter, test = test,
                              i = c(1,2,3,4,5,6,7,8,9,10),
                              k.index= 1:dim(parameter)[1],
                              j.index= 1:5, ylim = 0:1,
                              sample_size = 80)
ggsave(paste0(fig_path, "ZIG_full_80",".png"), d3_80, width = 20, height=21.6)

#5.3.2 full power plot - ZIG size 400
d3_400 <- pval.plot.power.full(result.stat.na, 
                               parameter_in_use = parameter, test = test,
                               i = c(1,2,3,4,5,6,7,8,9,10),
                               k.index= 1:dim(parameter)[1],
                               j.index= 1:5, ylim = 0:1,
                               sample_size = 400)
ggsave(paste0(fig_path, "ZIG_full_400",".png"), d3_400, width = 20, height=21.6)

##################### !!! original plot for full display!!! ##############################

if(FALSE){
  for (test in method.stat) {
    # a <- lapply(c(length.delta + 1, 1:length.delta), 
    #             function(i) pval.plot(result.stat, parameter_in_use = parameter, i=i, test=test))
    a <- lapply(1:length.delta, 
                function(i) pval.plot(result.stat, parameter_in_use = parameter, i=i, test=test))
    a <- marrangeGrob(a, nrow = 5, ncol = 2, top = NULL)
    #ggsave(paste0("Document/plot/P1101.zinb.", test, ".png"), a, width = 10, height=12)
    ggsave(paste0(fig_file, "-" ,test,".png"), a, width = 10, height=12)
  }
  for (test in method.stat) {
    # a <- lapply(c(length.delta+1, 1:length.delta), 
    a <- lapply(1:length.delta, 
                function(i) pval.plot(NA.proportion, parameter_in_use = parameter,
                                      i=i, test=test, title=paste0("NA proportion of ", test)))
    a <- marrangeGrob(a, nrow = 5, ncol = 2, top = NULL)
    #ggsave(paste0("Document/plot/P1101", model, ".", test,".png"), a, width = 10, height=12)
    ggsave(paste0(fig_file, "-_na_prop-", test,".png"), a, width = 10, height=12)
  }
}

##################### !!! reduced plots need to be modified !!! ##############################
## reduced plots
if(FALSE){
  #1. null effects
  a <- lapply(c("LB.glob", "LN", "MAST.glob", "KW",  "Wg.glob"), 
              function(test) pval.plot(result.stat.na, 
                                       parameter_in_use = parameter,
                                       i=1, test=test, 
                                       k.index= c(2,3,5,6,11,12,14,15),
                                       j.index= c(1,3), ylim = c(0,0.3),
                                       title = test.name[test == test.name[,1],2] %>% gsub(" \\- .*", "", .) ))
  a <- marrangeGrob(a, nrow=5, ncol=1, top = NULL)
  ggsave(paste0(fig_path, "result.", model, "_reduced" , "0null.png"), a, width = 5, height=12)
  
  #1.5. more null effects
  a <- lapply(c("LB.glob", "LN", "MAST.glob", "KW",  "Wg.glob"), 
              function(test) pval.plot(result.stat.na, 
                                       parameter_in_use = parameter,
                                       i=1, test=test, 
                                       k.index= c(1:32),
                                       j.index= c(1:5), ylim = c(0,0.3),
                                       title = test.name[test == test.name[,1],2] %>% gsub(" \\- .*", "", .) ))
  a <- marrangeGrob(a, nrow=5, ncol=1, top = NULL)
  ggsave(paste0(fig_path, "main_", model, "_null_reduced_" ,n,".png"), a, width = 5, height=12)
  
  #2. differential effects
  for (test in c("LB.glob", "LN", "MAST.glob", "KW",  "Wg.glob")) {
    a <- lapply(2:7, function(i) pval.plot(result.stat.na, 
                                           parameter_in_use = parameter,
                                           i=i, test=test, 
                                           k.index= c(2,3,5,6,11,12,14,15),
                                           j.index= c(1,3)))
    a <- marrangeGrob(a, nrow=3, ncol=2, top = NULL)
    ggsave(paste0(fig_path, "result.",model, "_reduced.", test,".png"), a, width = 10, height=12)
  }
  
  #2.1. reduced LN effects
  for (test in c("LN")) {
    a <- lapply(c(6, 2, 10, 3, 8), function(i) pval.plot(result.stat.na, 
                                                         parameter_in_use = parameter,
                                                         i=i, test=test, 
                                                         k.index= c(1:32),
                                                         j.index= c(1:5)))
    a <- marrangeGrob(a, nrow=5, ncol=1, top = NULL)
    ggsave(paste0(fig_path, "main_", model, "_", test, "_reduced_",n,".png"), a, width = 5, height=12)
  }
  
  #2.2. reduced LB, MAST effects
  for (test in c("LB.glob","MAST.glob")) {
    a <- lapply(c(6, 7, 3), function(i) pval.plot(result.stat.na, 
                                                  parameter_in_use = parameter,
                                                  i=i, test=test, 
                                                  k.index= c(1:32),
                                                  j.index= c(1:5)))
    a <- marrangeGrob(a, nrow=3, ncol=1, top = NULL)
    ggsave(paste0(fig_path, "main_", model, "_", test, "_reduced_",n,".png"), a, width = 5, height=7.2)
  }
  
  
  #2.3. reduced KW effects
  for (test in c("KW")) {
    a <- lapply(c(6, 2, 8, 3), function(i) pval.plot(result.stat.na, 
                                                     parameter_in_use = parameter,
                                                     i=i, test=test, 
                                                     k.index= c(1:32),
                                                     j.index= c(1:5)))
    a <- marrangeGrob(a, nrow=4, ncol=1, top = NULL)
    ggsave(paste0(fig_path, "main_", model, "_", test, "_reduced_",n,".png"), a, width = 5, height=9.6)
  }
  
  #2.4. reduced Wg.glob effects
  for (test in c("Wg.glob")) {
    a <- lapply(c(6, 8), function(i) pval.plot(result.stat.na, 
                                               parameter_in_use = parameter,
                                               i=i, test=test, 
                                               k.index= c(1:32),
                                               j.index= c(1:5)))
    a <- marrangeGrob(a, nrow=2, ncol=1, top = NULL)
    ggsave(paste0(fig_path, "main_", model, "_", test, "_reduced_  ",n,".png"), a, width = 5, height=4.8)
  }
}



# NA proprtions are almost the same across methods
NA.proportion %>% dplyr::filter(j%in%c(1,3), k%in% c(2,3,5,6,11,12,14,15)) %>% dplyr::select(i, j,  k, LB.glob, MAST.glob, LN, KW, Wg.glob, effect, batch)

if (FALSE) {
  result.stat.1 %>% dplyr::filter(i==1) # null effect            
  # LB suffers from zero-inflation(k=10~27)
  # LN and KW is robust to batch effects, but LB is not.
  
  result.stat.1 %>% dplyr::filter(i==2) # mean shift (nonzero)
  # LN has generally high power than KW
  # For high zero-inflation (p=0.95), LN and KW suffers (power = alpha),
  # For extreme zero-inflation (p=0.99), LN and KW has no power
  # KW suffers more than LN as zero-inflation gets higher.
  
  # When batch effect gets larger, overall  power decreases, but KW suffers more.
  
  
  result.stat.1 %>% dplyr::filter(i==3) # scale effect (keeping mean the same)
  result.stat.1 %>% dplyr::filter(i==4) # zero inflation
  result.stat.1 %>% dplyr::filter(i==5) # mean + scale
  result.stat.1 %>% dplyr::filter(i==6) # mean + zero inflation
  result.stat.1 %>% dplyr::filter(i==7) # scale + zero inflation
}

