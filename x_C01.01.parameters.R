### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra); library(cowplot)
library(MASS); library(boot); library(pscl); library(R.utils); library(latex2exp)
source("F00.00.generic.R")
source("F01.01.base.R")

### 0.2 Data
if (zoe == 1) {
  gene.marginal.RPK.DRNA <- readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE1.rds")
} else if (zoe == 2) {
  gene.marginal.RPK.DRNA <- readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE2.rds")
}
excluded.subject <- gene.marginal.RPK.DRNA$meta$id %in% c(352, 420, 10083, 11210, 11259, 11790, 12623)
DataMeta = gene.marginal.RPK.DRNA$meta[!excluded.subject,]
n.test   = 1000

# Raw data of 118 subjects
RNA     = gene.marginal.RPK.DRNA$otu[,, 2]
DNA     = gene.marginal.RPK.DRNA$otu[,, 1]
DataRPKRNA  = RNA[,colnames(RNA) %in% DataMeta$id]
DataRPKDNA  = DNA[,colnames(DNA) %in% DataMeta$id]
ST_RNA    = apply(DataRPKRNA, 2, sum)
ST_DNA= apply(DataRPKDNA, 2, sum)
mean(ST_RNA) # 5,551,718
DataTPM <- t(t(DataRPKRNA)/ST_RNA) * 5E+6
DataTPMDNA = t(t(DataRPKDNA)/ST_DNA) * 5E+6
rm(gene.marginal.RPK.DRNA, excluded.subject, RNA, DNA)

############################################################################################
# basic statistics and distribution
############################################################################################
dim(DataTPM) # 342,006 genes
DataRPKRNA %>% apply(2, sum) %>% mean # 5,551,718 average 'total RPKs per sample'

zero.proportion =
  DataRPKRNA %>% 
  apply(1, function(x) mean(x == 0)) 
zero.proportion %>% mean # avg zero proportion = 0.879

zero.proportion.DNA =
  DataRPKDNA %>% 
  apply(1, function(x) mean(x == 0)) 
zero.proportion.DNA %>% mean # avg zero proportion = 0.684
mean(zero.proportion.DNA <= 0.2) # 16%
mean(zero.proportion <= 0.2) # 3%

mean.var = 
  tibble(nz.mean.DNA = DataRPKDNA %>% apply(1, function(x) mean(x[x>0])),
         nz.var.DNA =  DataRPKDNA %>% apply(1, function(x)  var(x[x>0])),
         overdispersion.DNA = nz.var.DNA/nz.mean.DNA^2,
         nz.mean.RNA = DataRPKRNA %>% apply(1, function(x) mean(x[x>0])),
         nz.var.RNA =  DataRPKRNA %>% apply(1, function(x)  var(x[x>0])),
         overdispersion.RNA = nz.var.RNA/nz.mean.RNA^2)

mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.mean.DNA, nz.var.DNA, col = overdispersion.DNA)) + geom_point() + theme_classic() + 
  xlab("nonzero mean (DNA)") + ylab("nonzero variance (DNA)") +
  xlim(c(0, 100)) + ylim(c(0, 100))
# ggsave("figure/C0101overdispersion.png")

mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.mean.RNA, nz.var.RNA, col = overdispersion.RNA)) + geom_point() + theme_classic() + 
  xlab("nonzero mean (RNA)") + ylab("nonzero variance (RNA)") +
  xlim(c(0, 100)) + ylim(c(0, 100))

mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(overdispersion.DNA, overdispersion.RNA)) + geom_point() + theme_classic()
  # xlab("nonzero mean (RNA)") + ylab("nonzero variance (RNA)") +

mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.mean.DNA, nz.mean.RNA, col = overdispersion.RNA)) + geom_point() + theme_classic() + 
  # xlab("nonzero mean (RNA)") + ylab("nonzero variance (RNA)") +
  xlim(c(0, 100)) + ylim(c(0, 100))

mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.var.RNA/nz.mean.RNA, nz.var.DNA/nz.mean.DNA)) + geom_point() + theme_classic() + 
  xlim(c(0, 100)) + ylim(c(0, 100))

tmp.prop = with(mean.var, mean(nz.var.RNA/nz.mean.RNA^2 > nz.var.DNA/nz.mean.DNA^2, na.rm = TRUE))
set.seed(1)
mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.var.RNA/nz.mean.RNA^2, nz.var.DNA/nz.mean.DNA^2)) + geom_point() + theme_classic() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  xlab(TeX("$var (RNA| +) / \\mu (RNA| +)^2$")) +
  ylab(TeX("$var (DNA| +) / \\mu (DNA| +)^2$")) +
  annotate("text", x = 10, y = 5, label = paste0(round(tmp.prop*100, 1), "%"), col = "#F8766D") +
  annotate("text", y = 10, x = 5, label = paste0(round(100 - tmp.prop*100, 1), "%"), col = "#619CFF")
ggsave("figure/C0101overdispersion1.png")

tmp.prop = with(mean.var, mean(nz.var.RNA/nz.mean.RNA > nz.var.DNA/nz.mean.DNA, na.rm = TRUE))
set.seed(1)
mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.var.RNA/nz.mean.RNA, nz.var.DNA/nz.mean.DNA)) + geom_point() + theme_classic() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  xlim(c(0, 1000)) + ylim(c(0,1000)) +
  xlab(TeX("$var (RNA| +) / \\mu (RNA| +)$")) +
  ylab(TeX("$var (DNA| +) / \\mu (DNA| +)$")) +
  annotate("text", x = 900, y = 500, label = paste0(round(tmp.prop*100, 1), "%"), col = "#F8766D") +
  annotate("text", y = 900, x = 500, label = paste0(round(100 - tmp.prop*100, 1), "%"), col = "#619CFF")
ggsave("figure/C0101overdispersion2.png")

mean.var %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.var.RNA/nz.mean.RNA, nz.var.DNA/nz.mean.DNA)) + geom_point() + theme_classic() +
  geom_abline(slope = 1, intercept = 0, col = "red")
mean(mean.var$nz.var.RNA/mean.var$nz.mean.RNA > mean.var$nz.var.DNA/mean.var$nz.mean.DNA * 0.867 , na.rm = TRUE) # only 25% var(RNA) > var(DNA)
# 0.867 is the scale adjustment (the ratio below) for fair comparison btw DNA and RNA.

mean.var$nz.mean.RNA %>% mean(na.rm = T) # 47.17
mean.var$nz.mean.DNA %>% mean(na.rm = T) # 54.38
mean(mean.var$nz.mean.RNA, na.rm = T)/mean(mean.var$nz.mean.DNA, na.rm = T) # 0.867

mean.var$nz.mean.DNA %>% mean(na.rm=T)
mean.var$nz.mean.RNA %>% mean(na.rm=T)

zp.DNA <-
  zero.proportion.DNA %>% 
  data.frame(x = .) %>% 
  ggplot(aes(x)) + geom_histogram() + theme_classic() + 
  xlab("proportion of zeros") + ylab("frequency") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_y_continuous(limits = c(0, 1.6e+5)) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5,vjust=0.5)) +
  geom_vline(xintercept = 0.9, color = "red", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.95, y = .5e+5 * 1.2, col = "red",
           label = paste0(round(mean(zero.proportion.DNA >= 0.9) * 100, 0), "%")) +
  geom_vline(xintercept = 0.8, color = "orange", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.85, y = .5e+5 * .75, col = "orange",
           label = paste0(round(mean(zero.proportion.DNA >= 0.8) * 100, 0), "%")) +
  geom_vline(xintercept = 0.7, color = "cadetblue4", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.75, y = .5e+5 * .5, col = "cadetblue4",
           label = paste0(round(mean(zero.proportion.DNA >= 0.7) * 100, 0), "%")) +
  geom_vline(xintercept = 0.6, color = "darkslategray", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.65, y = .5e+5 * .35, col = "darkslategray",
           label = paste0(round(mean(zero.proportion.DNA >= 0.6) * 100, 0), "%"))
ggsave(zp.DNA, "figure/C0101zeroproportionDNA.png")

zp.RNA <-
  zero.proportion %>% 
  data.frame(x = .) %>% 
  ggplot(aes(x)) + geom_histogram() + theme_classic() + 
  xlab("proportion of zeros") + ylab(" ") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_y_continuous(limits = c(0, 1.6e+5)) +
  theme(axis.text.y=element_blank()) +
  geom_vline(xintercept = 0.9, color = "red", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.95, y = 1e+5 * 1.2, col = "red",
           label = paste0(round(mean(zero.proportion >= 0.9) * 100, 0), "%")) +
  geom_vline(xintercept = 0.8, color = "orange", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.85, y = 1e+5 * .75, col = "orange",
           label = paste0(round(mean(zero.proportion >= 0.8) * 100, 0), "%")) +
  geom_vline(xintercept = 0.7, color = "cadetblue4", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.75, y = 1e+5 * .5, col = "cadetblue4",
           label = paste0(round(mean(zero.proportion >= 0.7) * 100, 0), "%")) +
  geom_vline(xintercept = 0.6, color = "darkslategray", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.65, y = 1e+5 * .35, col = "darkslategray",
           label = paste0(round(mean(zero.proportion >= 0.6) * 100, 0), "%"))
ggsave(zp.RNA, "figure/C0101zeroproportion.png")

zp <-
  plot_grid(zp.DNA, zp.RNA + ylab(NULL), labels = "AUTO", rel_widths = c(1.12, 1), label_x = c(0.2, 0.1))
save_plot("figure/C0101zeroproportion.png", plot = zp, base_width = 9, base_height = 5)

zero.proportion %>% mean
zero.proportion %>% 
  cut(breaks = c(-0.1, 0.5, 0.8, 0.9, 0.95, 0.99, Inf)) %>% 
  table %>% {./length(zero.proportion)}


############################################################################################
# variances of DNA and RNA
############################################################################################
variances = data.frame(var.DNA = apply(DataRPKDNA, 1, var), var.RNA = apply(DataRPKRNA, 1, var),
                       mean.DNA = apply(DataRPKDNA, 1, mean), mean.RNA = apply(DataRPKRNA, 1, mean))
variances =
  variances %>% 
  mutate(disp.DNA = ifelse(mean.DNA < 1e-3, NA, var.DNA/mean.DNA),
         disp.RNA = ifelse(mean.RNA < 1e-3, NA, var.RNA/mean.RNA),
         RNAtoDNA = disp.RNA/ disp.DNA)
variances %>% 
  filter(var.DNA < 10000, var.RNA < 10000) %>% 
  ggplot(aes(var.DNA, var.RNA)) + geom_point()
variances %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(disp.DNA, disp.RNA)) + geom_point() + xlim(c(0, 2000)) + ylim(c(0, 2000))
variances$RNAtoDNA %>% {function(s) mean(s > 1, na.rm=TRUE)}(.)

var.nz <- function(x) var(x[x>0])
variances.nz = data.frame(var.DNA = apply(DataRPKDNA, 1, var.nz), var.RNA = apply(DataRPKRNA, 1, var.nz))
variances.nz %>% 
  filter(var.DNA < 10000, var.RNA < 10000) %>% 
  ggplot(aes(var.DNA, var.RNA)) + geom_point()



############################################################################################
#Neew to throw out DNA info                        
#filter batches
############################################################################################


# # filtering 2 lowly-expressed subjects (352, 420)
# gene.marginal.RPK.DRNA %>% select(-RPK.352, -RPK.420) -> DataRPKRNA
# outcome %>% filter(!id %in% c(352,420)) -> DataMeta
# rm(gene.marginal.RPK.RNA, outcome)

# disease groups (0:Healthy, 1:Treated. 2:Diseased)
HD = which(DataMeta$cariesfree == 0) #disease
H0 = which(DataMeta$cariesfree == 1) #healty

B1 = which(DataMeta$batch.RNA == "170628")
B2 = which(DataMeta$batch.RNA == "170718")

H0B1 = which(DataMeta$ECC == 0 & DataMeta$batch.RNA == "170628")
D2B1 = which(DataMeta$ECC == 2 & DataMeta$batch.RNA == "170628")
H0B2 = which(DataMeta$ECC == 0 & DataMeta$batch.RNA == "170718")
D2B2 = which(DataMeta$ECC == 2 & DataMeta$batch.RNA == "170718")



if(F){
  # subject_HD = DataMeta[ DataMeta$cariesfree == 0,]
  # subject_H0 = DataMeta[ DataMeta$cariesfree == 1,]
  # sample_HD = DataTPM[,colnames(DataTPM) %in% subject_HD$id]
  # sample_H0 = DataTPM[,colnames(DataTPM) %in% subject_H0$id]
  # 
  # subject_B1 = DataMeta[ DataMeta$batch.RNA == 170628,]
  # subject_B2 = DataMeta[ DataMeta$batch.RNA == 170718,]
  # sample_B1 = DataTPM[,colnames(DataTPM) %in% subject_B1$id]
  # sample_B2 = DataTPM[,colnames(DataTPM) %in% subject_B2$id]
  
  subject_H1 = DataMeta[ DataMeta$cariesfree == 1 & DataMeta$batch.RNA == "170628",]
  subject_H2 = DataMeta[ DataMeta$cariesfree == 1 & DataMeta$batch.RNA == "170718",]
  subject_D1 = DataMeta[ DataMeta$cariesfree == 0 & DataMeta$batch.RNA == "170628",]
  subject_D2 = DataMeta[ DataMeta$cariesfree == 0 & DataMeta$batch.RNA == "170718",]
  set.seed(1)
  samp = sample(1:nrow(DataTPM),1000)
  sample_H1 = DataTPM[samp,colnames(DataTPM) %in% subject_H1$id]
  sample_H2 = DataTPM[samp,colnames(DataTPM) %in% subject_H2$id]
  sample_D1 = DataTPM[samp,colnames(DataTPM) %in% subject_D1$id]
  sample_D2 = DataTPM[samp,colnames(DataTPM) %in% subject_D2$id]
}
#ZINB.ML takes too long


dataExamine <- function(data){
  mu = character(0)
  theta = character(0)
  pi = character(0)
  check = character(0)
  i_mu = 1; i_theta = 1; i_pi = 1;i_check = 1;
  
  for(i in 1:dim(data)[1]){
    if(i%%10==0){print(i)}
    y.vector = as.vector(data[i,])
    if(max(y.vector)>0){
      result = ZINB.ML.time(y.vector)
    }
    if(!is.nan(result[1]) & !is.nan(result[2]) & !is.nan(result[3]) &
       !is.na(result[1]) & !is.na(result[2]) & !is.na(result[3]) ){
      mu[i_mu] = result[1]
      i_mu = i_mu+1;
      
      theta[i_theta] = result[2]
      i_theta = i_theta+1;
      
      pi[i_pi] = result[3]
      i_pi = i_pi+1;
    }
  }
  return (list(as.numeric(mu),as.numeric(theta),as.numeric(pi)))
}



if(F){
  result_H1 = dataExamine(sample_H1)
  
  result_H2 = dataExamine(sample_H2)
  
  result_D1 = dataExamine(sample_D1)
  
  result_D2 = dataExamine(sample_D2)
  
  saveRDS(result_H1, "parameters/result_H1.rds")
  saveRDS(result_H2, "parameters/result_H2.rds")
  saveRDS(result_D1, "parameters/result_D1.rds")
  saveRDS(result_D2, "parameters/result_D2.rds")
  
  ####Analyze H1
  ## mu
  summary(unlist(result_H1[1]))
  boxplot(unlist(result_H1[1]), main = "mu of H1")
  boxplot(unlist(result_H1[1])[unlist(result_H1[1])<1000], main = "mu of H1(>1000 omitted)")
  boxplot(unlist(result_H1[1])[unlist(result_H1[1])<100], main = "mu of H1(>100 omitted)")
  boxplot(unlist(result_H1[1])[unlist(result_H1[1])<40], main = "mu of H1(>40 omitted)")
  plot(hist(unlist(result_H1[1])[unlist(result_H1[1])<40]), main = "mu of H1(>40 omitted)", xlab = "mu")
  ## theta
  
  ## pi
  
  ####Analyze H2
  
  ####Analyze D1
  
  ####Analyze D2
  
  # Comparing between HD and H0, mu_D > mu_H, theta_D > theta_H, pi_D < pi_H
  # Comparing between B1 and B2, mu_1 < mu_2, theta_1 > theta_2, pi_1 > pi_2
}

### 1. parameter estimates from real data (ZINB): baseline parameters
if (FALSE) {
  param = data.frame(gene.id = NA, mu=NA, theta=NA, pi=NA)
  n.gene = dim(DataTPM)[1]
  
  tt(1); k=1
  for (i in 1:n.gene) {
    if (i %% 100) next  #only every 100 other genes
    if (!i %% 3000) cat(i," out of ",n.gene, "genes (", round(i/n.gene*100,2),"%), expected time = ",
                        (Sys.time()-time.tmp)/i*n.gene, "minutes.")
    param[k,] = c(i, DataTPM[i,-1] %>% round %>% ZINB.ML.time(notation="mtp"))
    k = k+1
  }
  tt(2)  #14 mins for 2,386 genes
  saveRDS(param, "output/R0101.param.rds")
} else {
  param <- readRDS("output/R0101.param.rds")
  n.gene = dim(DataTPM)[1] 
}

## 1.1 resulting figures ####
param %>% apply(2, mean, na.rm=T)
param %>% ggplot(aes(mu, theta, col=pi)) + geom_point()

# parameters without outliers
param.outlier = which(param$theta > 1e+2)
param [param.outlier,]
DataTPM[119100,-1] %>% as.numeric %>% round %>% hist
length(param.outlier) # 714 out of 2,386 genes have theta > 100

param %>% filter(theta<=100, mu<500) %>% ggplot(aes(theta, mu, col=pi)) + geom_point() + ggtitle("theta<100")
ggsave("figure/P0101-param-mtp.png")
saveRDS(param, "param.rds")

param %>% filter(theta<100) %>% ggplot(aes(theta)) + geom_density() 
# param %>% filter(theta<100) %>% ggplot(aes(theta)) + geom_density() + xlim(c(0,1)); 
# mode at theta=.1
ggsave("figure/P0101-param-t.png")

param %>% filter(theta <=100, mu<500) %>% ggplot(aes(mu)) + geom_density()
#param %>% filter(theta <=100) %>% ggplot(aes(mu)) + geom_density() + xlim(c(0,3));
ggsave("figure/P0101-param-m.png")
# mode at mu=1

param %>% filter(theta <=100) %>% ggplot(aes(pi)) + geom_density()
ggsave("figure/P0101-param-p.png")
# uniform [0,1] + uniform [0.75,1] with 50% chance


### 2. parameter estimates from real data (ZINB): delta (H-D group differences)
if (FALSE) {
  param.ECC = data.frame(gene.id = NA, 
                         mu0=NA, theta0=NA, pi0=NA, 
                         mu1=NA, theta1=NA, pi1=NA,
                         mu2=NA, theta2=NA, pi2=NA)
  tt(1); k=1
  for (i in 1:n.gene) {
    if (i %% 100) next  #only every 100 other genes
    if (!i %% 3000) cat(i," out of ",n.gene, "genes (", round(i/n.gene*100,2),"%), expected time = ",
                        (Sys.time()-time.tmp)/i*n.gene, "minutes.")
    param.ECC[k,] = c(i, ZINB.ML.time(DataTPM[i,H0]),
                      ZINB.ML.time(DataTPM[i,T1]), ZINB.ML.time(DataTPM[i,D2]))
    # filling in pi's with the proportion of zero counts, if the genes are not estimated.
    if (is.na(param.ECC[k,4])) {param.ECC[k,4] = mean(DataTPM[i,H0]==0)}
    if (is.na(param.ECC[k,7])) {param.ECC[k,7] = mean(DataTPM[i,T1]==0)}
    if (is.na(param.ECC[k,10])) {param.ECC[k,10] = mean(DataTPM[i,D2]==0)}
    k = k+1
  }
  tt(2)  # 26 mins for 2,386 genes
  saveRDS(param.ECC, "output/R0101.param.ECC.rds")
} else {
  param.ECC <- readRDS("output/R0101.param.ECC.rds")
}    


## 2.1 error analysis ####   
# NaN: estimation error
# NA: timeouts
apply(param.ECC, 2, function(x) mean(is.nan(x))) %>% round(2) -> tmp.NaN
apply(param.ECC, 2, function(x) mean(is.na(x))) %>% round(2) -> tmp.NA
tmp.NA = tmp.NA - tmp.NaN #NA includes NaN. Thus subtracted.
rbind(tmp.NaN, tmp.NA)[,c(2,5,8)]
#           H      T      D
# tmp.NaN   0.32   0.43   0.24
# tmp.NA    0.09   0.06   0.11

# when all (or all but a few) are zero counts, error
# when the few nonzero counts are large numbers, error (NaN)
# when the few nonzero counts are small numbers, time out (NA)
# 
# param.ECC[is.nan(param.ECC[,2]),1] %>% sapply(function(s) {mean(DataTPM[s,H0+1]==0)})
# param.ECC[is.nan(param.ECC[,5]),1] %>% sapply(function(s) {mean(DataTPM[s,T1+1]==0)})
# param.ECC[is.nan(param.ECC[,8]),1] %>% sapply(function(s) {mean(DataTPM[s,D2+1]==0)})
# Seen that most of the NaN are mostly zero counts
param.ECC[(!is.nan(param.ECC[,2])) & is.na(param.ECC[,2]) ,1] %>% sapply(function(s) {mean(DataTPM[s, H0]==0)})
DataTPM[1200,H0] %>% round %>% sort %>% as.numeric # example of timeout counts
DataTPM[2300,H0] %>% round %>% sort %>% as.numeric # example of timeout counts
DataTPM[3000,H0] %>% round %>% sort %>% as.numeric # example of timeout counts

# NaN and NA are very similar in that there are overwelming # of zero's.
# -> For those NaN and NA, pi's are filled in with the zero-proportion.

## 2.2 delta (ratio) #####
## delta.theta
data_frame(x=param.ECC[,2],y=param.ECC[,8], delta=pmax(y/x, x/y)) %>% 
  na.omit %>% arrange(delta) %>% mutate(index = 1:n()) %>% 
  mutate(facet=cut(index, c(0,900,1000,1200))) %>%
  mutate(facet=factor(facet, levels = levels(facet)[3:1])) -> tmp

log(tmp$delta) %>%  summary # 3Q: 0.72

tmp %>% ggplot(aes(index, delta)) + geom_point() + facet_grid(facet ~., scale="free") +
  ggtitle ("delta_theta of real data")
ggsave("figure/P0102-param-delta-t.png")
# suggests 5?

# distribution among reasonable parameter values (theta<=100)
tmp %<>% filter(x<=100, y<=100) 
tmp %>% ggplot(aes(x=delta)) + geom_histogram() + xlab("delta_theta") + ggtitle("theta<=100")
ggsave("figure/P0102-param-delta-t100.png")
tmp %>% filter(delta %btw% c(5,30)) %>% ggplot(aes(x,y, col=delta)) + geom_point() +
  ggtitle("theta estimate distribution where delta is between 5 and 30, theta<=100")+ xlab("theta(ECC0)") + ylab("theta(ECC2)")
ggsave("figure/P0102-param-delta-t100b.png")  


## delta.mu
data_frame(x=param.ECC[,3],y=param.ECC[,9], delta=pmax(y/x, x/y)) %>% 
  na.omit %>% arrange(delta) %>% mutate(index = 1:n()) %>% 
  mutate(facet=cut(index, c(0,800,1000,1200))) %>%
  mutate(facet=factor(facet, levels = levels(facet)[3:1])) -> tmp

log(tmp$delta) %>%  summary # 3Q: 3.66

tmp %>% ggplot(aes(index, delta)) + geom_point() + facet_grid(facet ~., scale="free") +
  ggtitle ("delta_mu of real data")
ggsave("figure/P0102-param-delta-m.png")
# suggests 20?

## distribution among reasonable parameter values (mu<=100)
tmp %<>% filter(x<=100, y<=100)
tmp %>% ggplot(aes(x=delta)) + geom_histogram() + xlab("delta_mu")  + ggtitle("mu<=100")
ggsave("figure/P0102-param-delta-m100.png")
tmp %>% filter(delta %btw% c(5,30)) %>% ggplot(aes(x, y, col=delta)) + geom_point() +
  ggtitle("mu estimate distribution where delta is between 5 and 30, mu<=100") + xlab("mu(ECC0)") + ylab("mu(ECC2)")
ggsave("figure/P0102-param-delta-m100b.png")



## delta.pi
data_frame(x=param.ECC[,4],y=param.ECC[,10], delta=pmax(y/x, x/y)) %>% 
  na.omit %>% arrange(delta) %>% mutate(index = 1:n()) %>% 
  mutate(facet=cut(index, c(0,2200,2300,3000))) %>%
  mutate(facet=factor(facet, levels = levels(facet)[3:1])) -> tmp
tmp %>% transmute(delta = qlogis(x) - qlogis(y)) %>% abs  %>% filter(is.finite(delta)) %>% summary
# 3Q = 0.99

tmp %>% ggplot(aes(index, delta)) + geom_point() + facet_grid(facet ~., scale="free") +
  ggtitle ("delta_pi of real data")
ggsave("figure/P0102-param-delta-p.png")
# suggests 3?

## distribution among reasonable parameter values (pi>1e-4)
tmp %<>% filter(x>1e-4, y>1e-4)
tmp %>% ggplot(aes(x=delta)) + geom_histogram() + xlab("delta_pi")  + ggtitle("pi>1e-4")
ggsave("figure/P0102-param-delta-p100.png")
tmp %>% ggplot(aes(x, y, col=delta)) + geom_point() +
  ggtitle("pi estimate distribution for all delta, pi>1e-4") + xlab("pi(ECC0)") + ylab("pi(ECC2)")
ggsave("figure/P0102-param-delta-p100b.png")



### 3. parameter estimates from re al data (ZINB): kappa (batch effects)    
if (FALSE) {
  param.bacth = data.frame(gene.id = NA, 
                           mu1=NA, theta1=NA, pi1=NA, 
                           mu2=NA, theta2=NA, pi2=NA)
  tt(1); k=1
  for (i in 1:n.gene) {
    if (i %% 100) next  #only every 100 other genes
    if (!i %% 3000) cat(i," out of ",n.gene, "genes (", round(i/n.gene*100,2),"%), expected time = ",
                        (Sys.time()-time.tmp)/i*n.gene, "minutes.")
    param.bacth[k,] = c(i, ZINB.ML.time(DataTPM[i,B1]), ZINB.ML.time(DataTPM[i,B2]))
    # filling in pi's with the proportion of zero counts, if the genes are not estimated.
    if (is.na(param.bacth[k,4])) {param.bacth[k,4] = mean(DataTPM[i,B1]==0)}
    if (is.na(param.bacth[k,7])) {param.bacth[k,7] = mean(DataTPM[i,B2]==0)}
    k = k+1
  }
  tt(2)  # 39 mins for 2,386 genes
  saveRDS(param.bacth, "output/R0101.param.batch.rds")
} else {
  param.bacth <- readRDS("output/R0101.param.batch.rds")
}

## 3.2 kappa (ratio) #####
## kappa.theta
data_frame(x=param.bacth[,2],y=param.bacth[,5], kappa=pmax(y/x, x/y)) %>% 
  na.omit %>% arrange(kappa) %>% mutate(index = 1:n()) %>% 
  mutate(facet=cut(index, c(0,800,1100,3000))) %>%
  mutate(facet=factor(facet, levels = levels(facet)[3:1])) -> tmp

log(tmp$kappa) %>%  summary # Med: 0.57

tmp %>% ggplot(aes(index, kappa)) + geom_point() + facet_grid(facet ~., scale="free") +
  ggtitle ("kappa_theta of real data")
ggsave("figure/P0102-param-kappa-t.png")

log(tmp$kappa) %>% mean %>% exp #1.8
# suggests 2

# distribution among reasonable parameter values (theta<=100)
tmp %<>% filter(x<=100, y<=100) 
tmp %>% ggplot(aes(x=kappa)) + geom_histogram() + xlab("kappa_theta") + ggtitle("theta<=100")
ggsave("figure/P0102-param-kappa-t100.png")
tmp %>% filter(kappa %btw% c(5,30)) %>% ggplot(aes(x,y, col=kappa)) + geom_point() +
  ggtitle("theta estimate distribution where kappa is between 5 and 30, theta<=100")+ xlab("theta(ECC0)") + ylab("theta(ECC2)")
ggsave("figure/P0102-param-kappa-t100b.png")  
log(tmp$kappa) %>% mean %>% exp #1.7
# suggests 2

## kappa.mu
data_frame(x=param.ECC[,3],y=param.ECC[,6], kappa=pmax(y/x, x/y)) %>% 
  na.omit %>% arrange(kappa) %>% mutate(index = 1:n()) %>% 
  mutate(facet=cut(index, c(0,650,850,1000))) %>%
  mutate(facet=factor(facet, levels = levels(facet)[3:1])) -> tmp

log(tmp$kappa) %>%  summary # Med: 1.07

tmp %>% ggplot(aes(index, kappa)) + geom_point() + facet_grid(facet ~., scale="free") +
  ggtitle ("kappa_mu of real data")
ggsave("figure/P0102-param-kappa-m.png")
log(tmp$kappa) %>% mean %>% exp #34
# suggests 30?

## distribution among reasonable parameter values (mu<=100)
tmp %<>% filter(x<=100, y<=100)
tmp %>% ggplot(aes(x=kappa)) + geom_histogram() + xlab("kappa_mu")  + ggtitle("mu<=100")
ggsave("figure/P0102-param-kappa-m100.png")
tmp %>% filter(kappa %btw% c(5,30)) %>% ggplot(aes(x, y, col=kappa)) + geom_point() +
  ggtitle("mu estimate distribution where kappa is between 5 and 30, mu<=100") + xlab("mu(ECC0)") + ylab("mu(ECC2)")
ggsave("figure/P0102-param-kappa-m100b.png")
log(tmp$kappa) %>% mean %>% exp #2.5


## kappa.pi
data_frame(x=param.ECC[,4],y=param.ECC[,7], kappa=pmax(y/x, x/y)) %>% 
  na.omit %>% arrange(kappa) %>% mutate(index = 1:n()) %>% 
  mutate(facet=cut(index, c(0,2200,2250,3000))) %>%
  mutate(facet=factor(facet, levels = levels(facet)[3:1])) -> tmp

tmp %>% transmute(kappa = qlogis(x) - qlogis(y)) %>% abs  %>% filter(is.finite(kappa)) %>% summary
# Median = 0.72

tmp %>% ggplot(aes(index, kappa)) + geom_point() + facet_grid(facet ~., scale="free") +
  ggtitle ("kappa_pi of real data")
ggsave("figure/P0102-param-kappa-p.png")
log(tmp$kappa) %>% mean %>% exp #Inf
# suggests 2?

## distribution among reasonable parameter values (pi>1e-4)
tmp %<>% filter(x>1e-4, y>1e-4)
tmp %>% ggplot(aes(x=kappa)) + geom_histogram() + xlab("kappa_pi")  + ggtitle("pi>1e-4")
ggsave("figure/P0102-param-kappa-p100.png")
tmp %>% ggplot(aes(x, y, col=kappa)) + geom_point() +
  ggtitle("pi estimate distribution for all kappa, pi>1e-4") + xlab("pi(ECC0)") + ylab("pi(ECC2)")
ggsave("figure/P0102-param-kappa-p100b.png")
log(tmp$kappa) %>% mean %>% exp #1.12



### 4. parameter estimates from real data (ZINB): grouped    
if (FALSE) {
  param.group = data.frame(gene.id = NA, 
                           mu_H0B1=NA, theta_H0B1=NA, pi_H0B1=NA, 
                           mu_D2B1=NA, theta_D2B1=NA, pi_D2B1=NA,
                           mu_H0B2=NA, theta_H0B2=NA, pi_H0B2=NA,
                           mu_D2B2=NA, theta_D2B2=NA, pi_D2B2=NA)
  tt(1); k=1
  for (i in 1:n.gene) {
    if (i %% 100) next  #only every 100 other genes
    if (!i %% 3000) cat(i," out of ",n.gene, "genes (", round(i/n.gene*100,2),"%), expected time = ",
                        (Sys.time()-time.tmp)/i*n.gene, "minutes.")
    param.group[k,] = c(i, ZINB.ML.time(DataTPM[i,H0B1]), ZINB.ML.time(DataTPM[i,D2B1]),
                        ZINB.ML.time(DataTPM[i,H0B2]), ZINB.ML.time(DataTPM[i,D2B2]))
    # filling in pi's with the proportion of zero counts, if the genes are not estimated.
    if (is.na(param.group[k,4])) {param.group[k,4] = mean(DataTPM[i,H0B1]==0)}
    if (is.na(param.group[k,7])) {param.group[k,7] = mean(DataTPM[i,D2B1]==0)}
    if (is.na(param.group[k,10])) {param.group[k,10] = mean(DataTPM[i,H0B2]==0)}
    if (is.na(param.group[k,13])) {param.group[k,13] = mean(DataTPM[i,D2B2]==0)}
    k = k+1
  }
  tt(2)  # 39 mins for 2,386 genes
  saveRDS(param.group, "output/R0101.param.group.rds")
} else {
  param.group <- readRDS("output/R0101.param.group.rds")
}

#################### p1 base parameter estimation ####################
param.group.H0B1 <- param.group %>% 
  mutate(disease = "H", batch = 1, group = "H1", mu = mu_H0B1, tht = theta_H0B1, pi = pi_H0B1) %>%
  dplyr::select(-(1:13)) %>%
  filter((pi > 0.27 & pi < 0.33) | (pi > 0.57 & pi < 0.63) | (pi > 0.87 & pi < 0.93)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((pi > 0.27 & pi < 0.33), 0.3, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.57 & pi < 0.63), 0.6, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.87 & pi < 0.93), 0.9, pi_id)) 

param.group.D2B1 <- param.group %>% 
  mutate(disease = "D", batch = 1, group = "D1", mu = mu_D2B1, tht = theta_D2B1, pi = pi_D2B1) %>%
  dplyr::select(-(1:13)) %>%
  filter((pi > 0.27 & pi < 0.33) | (pi > 0.57 & pi < 0.63) | (pi > 0.87 & pi < 0.93)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((pi > 0.27 & pi < 0.33), 0.3, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.57 & pi < 0.63), 0.6, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.87 & pi < 0.93), 0.9, pi_id)) 

param.group.H0B2 <- param.group %>% 
  mutate(disease = "H", batch = 2, group = "H2", mu = mu_H0B2, tht = theta_H0B2, pi = pi_H0B2) %>%
  dplyr::select(-(1:13)) %>%
  filter((pi > 0.27 & pi < 0.33) | (pi > 0.57 & pi < 0.63) | (pi > 0.87 & pi < 0.93)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((pi > 0.27 & pi < 0.33), 0.3, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.57 & pi < 0.63), 0.6, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.87 & pi < 0.93), 0.9, pi_id)) 

param.group.D2B2 <- param.group %>% 
  mutate(disease = "D", batch = 2, group = "D2", mu = mu_D2B2, tht = theta_D2B2, pi = pi_D2B2) %>%
  dplyr::select(-(1:13)) %>%
  filter((pi > 0.27 & pi < 0.33) | (pi > 0.57 & pi < 0.63) | (pi > 0.87 & pi < 0.93)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((pi > 0.27 & pi < 0.33), 0.3, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.57 & pi < 0.63), 0.6, pi_id)) %>%
  mutate(pi_id = ifelse((pi > 0.87 & pi < 0.93), 0.9, pi_id)) 

param.group.all = rbind(param.group.H0B1, param.group.D2B1, param.group.H0B2, param.group.D2B2)

param.group.all$pi_id_f = factor(param.group.all$pi_id,
                               levels = c(0.3, 0.6, 0.9),
                               labels = c(TeX("$\\pi \\approx 0.3$"), 
                                          TeX("$\\pi \\approx 0.6$"),
                                          TeX("$\\pi \\approx 0.9$")))

param.group.all %>%
  dplyr::filter(tht>0 & tht < 30 & mu < 100) %>%
  ggplot(aes(mu, tht, shape = group, col = group)) +
  geom_point() +
  scale_color_manual(name = "group",
                     values=c(D1 = "skyblue", 
                              D2 = "deepskyblue4",  
                              H1 = "tomato3", 
                              H2 = "tomato3")) +
  scale_shape_manual(name = "group",values=c(17, 15, 2, 0)) +
  xlab(TeX('$\\mu$')) + 
  ylab(TeX('$\\theta')) + 
  ggtitle("baseline parameter estimates") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> para_zinb_1


#################### p2 delta estimation (healthy/diseased) ####################

param.group.delta.B1 <- param.group %>% 
  mutate(delta_mu = pmax(mu_H0B1/mu_D2B1, mu_D2B1/mu_H0B1),
         delta_tht = pmax(theta_H0B1/theta_D2B1, theta_D2B1/theta_H0B1),
         delta_pi = pmax(pi_H0B1/pi_D2B1, pi_D2B1/pi_H0B1),
         batch = '1') %>% 
  dplyr::select(-(1:13)) %>% na.omit() %>%
  filter((delta_pi < 1.02) | 
           (delta_pi > 3 & delta_pi < 5) | 
           (delta_pi > 6 & delta_pi < 8)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((delta_pi < 1.02), 1, pi_id)) %>%
  mutate(pi_id = ifelse((delta_pi > 3 & delta_pi < 5), 4, pi_id)) %>%
  mutate(pi_id = ifelse((delta_pi > 6 & delta_pi < 8), 7, pi_id)) 

param.group.delta.B2 <- param.group %>% 
  mutate(delta_mu = pmax(mu_H0B2/mu_D2B2, mu_D2B2/mu_H0B2),
         delta_tht = pmax(theta_H0B2/theta_D2B2, theta_D2B2/theta_H0B2),
         delta_pi = pmax(pi_H0B2/pi_D2B2, pi_D2B2/pi_H0B2),
         batch = '2') %>% 
  dplyr::select(-(1:13)) %>% na.omit() %>%
  filter((delta_pi < 1.02) | 
           (delta_pi > 3 & delta_pi < 5) | 
           (delta_pi > 6 & delta_pi < 8)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((delta_pi < 1.02), 1, pi_id)) %>%
  mutate(pi_id = ifelse((delta_pi > 3 & delta_pi < 5), 4, pi_id)) %>%
  mutate(pi_id = ifelse((delta_pi > 6 & delta_pi < 8), 7, pi_id)) 

param.group.delta.all = rbind(param.group.delta.B1, param.group.delta.B2)

param.group.delta.all$pi_id_f = factor(param.group.delta.all$pi_id,
                                 levels = c(1, 4, 7),
                                 labels = c(TeX("$\\pi \\approx 1$"), 
                                            TeX("$\\pi \\approx 4$"),
                                            TeX("$\\pi \\approx 7$")))

param.group.delta.all %>%
  dplyr::filter(delta_tht > 0 & delta_tht < 30 & delta_mu < 10) %>%
  ggplot(aes(delta_mu, delta_tht, col = batch, shape = batch)) +
  geom_point() +
  scale_shape_manual(name = "batch",values=c(17, 15)) +
  xlab(TeX('$\\delta_\\mu$')) + 
  ylab(TeX('$\\delta_\\theta')) + 
  ggtitle(TeX("(|$\\delta_{\\mu}$|, |$\\delta_{\\theta}$|, |$\\delta_{\\pi}$|) estimates")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> para_zinb_2



#################### p3 kappa estimation (b1/b2) ####################
param.group.kappa.healthy <- param.group %>% 
  mutate(kappa_mu = pmax(mu_H0B1/mu_H0B2, mu_H0B2/mu_H0B1),
         kappa_tht = pmax(theta_H0B1/theta_H0B2, theta_H0B2/theta_H0B1),
         kappa_pi = pmax(pi_H0B1/pi_H0B2, pi_H0B2/pi_H0B1),
         disease = 'healthy') %>% 
  dplyr::select(-(1:13)) %>% na.omit() %>%
  filter((kappa_pi < 1.02) | 
           (kappa_pi > 3 & kappa_pi < 5) | 
           (kappa_pi > 6 & kappa_pi < 8)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((kappa_pi < 1.02), 1, pi_id)) %>%
  mutate(pi_id = ifelse((kappa_pi > 3 & kappa_pi < 5), 4, pi_id)) %>%
  mutate(pi_id = ifelse((kappa_pi > 6 & kappa_pi < 8), 7, pi_id)) 

param.group.kappa.diseased <- param.group %>% 
  mutate(kappa_mu = pmax(mu_D2B1/mu_D2B2, mu_D2B2/mu_D2B1),
         kappa_tht = pmax(theta_D2B1/theta_D2B2, theta_D2B2/theta_D2B1),
         kappa_pi = pmax(pi_D2B1/pi_D2B2, pi_D2B2/pi_D2B1),
         disease = 'diseased') %>% 
  dplyr::select(-(1:13)) %>% na.omit() %>%
  filter((kappa_pi < 1.02) | 
           (kappa_pi > 3 & kappa_pi < 5) | 
           (kappa_pi > 6 & kappa_pi < 8)) %>%
  mutate(pi_id = 0) %>%
  mutate(pi_id = ifelse((kappa_pi < 1.02), 1, pi_id)) %>%
  mutate(pi_id = ifelse((kappa_pi > 3 & kappa_pi < 5), 4, pi_id)) %>%
  mutate(pi_id = ifelse((kappa_pi > 6 & kappa_pi < 8), 7, pi_id)) 

param.group.kappa.all = rbind(param.group.kappa.healthy, param.group.kappa.diseased)

param.group.kappa.all$pi_id_f = factor(param.group.kappa.all$pi_id,
                                       levels = c(1, 4, 7),
                                       labels = c(TeX("$\\pi \\approx 1$"), 
                                                  TeX("$\\pi \\approx 4$"),
                                                  TeX("$\\pi \\approx 7$")))

param.group.kappa.all %>%
  dplyr::filter(kappa_tht > 0 & kappa_tht < 30 & kappa_mu < 5) %>%
  ggplot(aes(kappa_mu, kappa_tht, col = disease, shape = disease)) +
  geom_point() +
  scale_shape_manual(name = "disease",values=c(16,1)) +
  xlab(TeX('$\\kappa_\\mu$')) + 
  ylab(TeX('$\\kappa_\\theta')) + 
  ggtitle(TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> para_zinb_3

grid.arrange(
  grobs = list(para_zinb_1, para_zinb_2, para_zinb_3),
  ncol = 3, nrow = 1
) -> p_zinb

ggsave(("figure/para_selection_zinb.png"),  p_zinb,  width = 10, height = 10)#, dpi = 200)



######## exercise
if (FALSE) {
  DataTPM[5,-1] %>% as.numeric %>% round %>% apply(1, ZINB.ML)
  DataTPM[5,-1] %>% as.numeric %>% round %>% hist
  DataTPM[5,-1] %>% as.numeric %>% round %>% ZINB.ML.time (timeout=1)
  DataTPM[6,-1] %>% as.numeric %>% round %>% ZINB.ML.time (timeout=1)
  
  # running 10 genes for example
  param = DataTPM[1:10,-1] %>% round %>% apply(1, ZINB.ML.time) %>% t
  param = structure(param, dimnames=list(1:10, c("theta", "mu", "pi")))
  
}

