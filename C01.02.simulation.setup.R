### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra); require(tidyr)
source("F00.00.generic.R")
source("F01.01.base.R")
source("F01.02.models.R")

### 1.0 simulation parameters  
# dataset-wise parameters
# n = 120; nD <- nT <- nH <- 40; 

# simulation-wise parameters
M = 10 # Mast replicates
cutoff = 0.9; sig = 0.05 # % NA < 0.9, p <= 0.05
method.stat = tester.set.HD.batch(n.gene=5, skeleton=TRUE)[[1]] %>% rownames
# LB.nonz LB.zero LB.glob LN MAST.nonz MAST.zero MAST.glob KW Wg.nonz Wg.zero Wg.glob DESeq2 (Reserved)
method = gsub("\\..*$","",method.stat) %>% unique
# LB LN MAST KW Wagner DESeq2 (spare)

### 2.0 distribution parameters
# parameter1 = basic scenarios

expand.grid(m=c(2, 3, 5, 10), t=c(0.5, 1), p=c(.3, .5, .6, .9, .95)) %>% # normal scenarios 1-40
  rbind(expand.grid(m=c(2, 3, 5), t=c(5), p=c(.9, .95))) %>%  # extreme scenario 41-46
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameter1

expand.grid(m=c(2, 3, 5, 10), t=c(0.5, 1), p=c(.65, .7, .75, .8, .85 )) %>% # normal scenarios 
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameter2

expand.grid(m=c(2, 3, 5, 10), t=c(0.5, 1), p=c(.3, .5, .6 ,.7 ,.8, .9, .95)) %>% # normal scenarios 1-40
  rbind(expand.grid(m=c(2, 3, 5), t=c(5), p=c(.9, .95))) %>%  # extreme scenario 41-46
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameter3

expand.grid(m=c(2, 5, 10), t=c(0.5, 1), p=c(.3, .6, .65, .7, .75, .8, .85, .9, .95)) %>% # normal scenarios 1-30
  #rbind(expand.grid(m=c(2, 5), t=c(5), p=c(.9, .95))) %>%  # extreme scenario 31-34
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameter4


expand.grid(m=c(1, 2, 3), t=c(0.5, 1), p=c(.3, .6 , .8, .9, .95)) %>% # normal scenarios 1-30 (t = sigma)
  rbind(expand.grid(m=c(2), t=c(0.1), p=c(.9, .95))) %>%  # extreme scenario 31-32
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameterLN1

expand.grid(m=c(2, 4, 6), t=c(1, 3), p=c(.3, .6, .65, .7, .75, .8, .85, .9, .95)) %>% # normal scenarios 1-30 (t = sigma)
  #rbind(expand.grid(m=c(6), t=c(5), p=c(.9, .95))) %>%  # extreme scenario 31-32
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameter5


expand.grid(m=c(1, 10, 20), t=c(0.2, 1), p=c(.3, .6, .65, .7, .75, .8, .85, .9, .95)) %>% # normal scenarios 1-54 (t = sigma)
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameterLN2


expand.grid(m=c(1, 10, 20), t=c(1, 5), p=c(.3, .6, .65, .7, .75, .8, .85, .9, .95)) %>% # normal scenarios 1-54
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameterNB2



# multipliers

# delta = differential expression
matrix(c(0, 0, 0,
         1, 0, 0,  0, 1,  0,    0, 0, -1,
         1, 1, 0,  1, 0, -1,    0, 1, -1,
         1, 0, 1,  1, -1, 0,    0, -1, -1), byrow = TRUE,
       nrow = 10, dimnames=list(1:10, c("m", "t", "p"))) %>% 
  data.frame ->
  delta1
delta1$detail = paste0("Effect_", c("null", "mu(D>H)", "theta(D>H)", "pi(D<H)", 
                                   "mu(D>H).theta(D>H)", "mu(D>H).pi(D<H)", "theta(D>H).pi(D<H)",
                                   "mu(D>H),pi(D>H)","mu(D>H).theta(D<H)","theta(D<H).pi(D<H)"))

matrix(c(0, 0, 0,   .5, -.5, -.5,   1, -1, -1,   .5, .5, -.5,   1, 1, -1), byrow = TRUE,
       nrow = 5, dimnames=list(1:5, c("m", "t", "p"))) %>% 
  data.frame -> 
  kappa1
kappa1$detail = paste(c("no","small(+,-,-)", "large(+,-,-)", "small(+,+,-)", "large(+,+,-)"), "batch effect")

