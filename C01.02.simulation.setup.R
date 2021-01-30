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
# LB.nonz LB.zero LB.glob LN MAST.nonz MAST.zero MAST.glob KW Wg.nonz Wg.zero Wg.glob DESeq2 MGS (Reserved)
method = gsub("\\..*$","",method.stat) %>% unique
# LB LN MAST KW Wagner DESeq2 MGS (spare)

### 2.0 distribution parameters
# parameter1 = basic scenarios

expand.grid(m=c(1, 10, 50), t=c(0.5, 2), p=c(.3, .6, .65, .7, .75, .8, .85, .9, .95)) %>% # normal scenarios 1-54 (t = sigma)
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameterLN


expand.grid(m=c(1, 10, 50), t=c(1, 5), p=c(.3, .6, .65, .7, .75, .8, .85, .9, .95)) %>% # normal scenarios 1-54
  data.frame() %>%
  mutate(no = 1:n()) %>%       # add scenario numbers (no)
  dplyr::select(no, everything())  -> # reorder columns
  parameterNB



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

