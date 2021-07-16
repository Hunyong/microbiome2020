### 0.1 library
library(dplyr); library(magrittr); library(ggplot2); library(gridExtra); library(cowplot)
library(MASS); library(boot); library(pscl); library(R.utils); library(latex2exp)
source("F00.00.generic.R")
source("F01.01.base.R")

type = "gene"
zoe = 0
### 0.2 Data
if (zoe == 1) {
  gene.marginal.RPK.DRNA <- readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE1.rds")
} else if (zoe == 2) {
  gene.marginal.RPK.DRNA <- readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE2.rds")
} else if (zoe == 0) {
  gene.marginal.RPK.DRNA <- readRDS("../MicrobiomePaper2020/Nature2019data/data.geneRPK.marginal.DRNA.IBD.rds")
  gene.marginal.RPK.DRNA$meta$id <- gene.marginal.RPK.DRNA$meta$External.ID
}
zoe.nm = if (zoe %in% 1:2) paste0("_zoe", zoe) else "_NEWDATA"

excluded.subject <- gene.marginal.RPK.DRNA$meta$id %in% c(352, 420, 10083, 11210, 11259, 11790, 12623)
DataMeta = gene.marginal.RPK.DRNA$meta[!excluded.subject,]
n.test   = 1000

# Raw data of 118 subjects
RNA     = gene.marginal.RPK.DRNA$otu[,, 2]
DNA     = gene.marginal.RPK.DRNA$otu[,, 1]
DataRPKRNA  = RNA[,colnames(RNA) %in% DataMeta$id]
DataRPKDNA  = DNA[,colnames(DNA) %in% DataMeta$id]
ST_RNA    = apply(DataRPKRNA, 2, sum, na.rm = TRUE)
ST_DNA    = apply(DataRPKDNA, 2, sum, na.rm = TRUE)
mean(ST_RNA, na.rm = TRUE) # 5,551,718
DataTPM <- t(t(DataRPKRNA)/ST_RNA) * 5E+6
DataTPMDNA = t(t(DataRPKDNA)/ST_DNA) * 5E+6
rm(gene.marginal.RPK.DRNA, excluded.subject, RNA, DNA)

############################################################################################
# basic statistics and distribution
############################################################################################
dim(DataTPM) # 342,006 genes
DataRPKRNA %>% apply(2, sum, na.rm = TRUE) %>% mean(na.rm = TRUE) # 5,551,718 average 'total RPKs per sample'

zero.proportion =
  DataRPKRNA %>% 
  apply(1, function(x) mean(x == 0, na.rm = TRUE)) 
zero.proportion %>% mean(na.rm = TRUE) # avg zero proportion = 0.879 0.804 / 0.963

zero.proportion.DNA =
  DataRPKDNA %>% 
  apply(1, function(x) mean(x == 0, na.rm = TRUE)) 
zero.proportion.DNA %>% mean(na.rm = TRUE) # avg zero proportion = 0.684 0.750 0.878

mean(zero.proportion.DNA <= 0.2, na.rm = TRUE) # 16% 11% / IBD 1.1%
mean(zero.proportion <= 0.2, na.rm = TRUE)     #  3%  6% / IBD 0.2%
mean(zero.proportion.DNA >= 0.9, na.rm = TRUE) # 43% 53.6% / IBD 68.5%
mean(zero.proportion >= 0.9, na.rm = TRUE)     # 71% 58.9% / IBD 90.7%
mean(zero.proportion.DNA >= 0.95, na.rm = TRUE) # xx% 44.7% / IBD 55.8%
mean(zero.proportion >= 0.95, na.rm = TRUE)     # xx% 50.2% / IBD 85.5%


mean.var = 
  tibble(nz.mean.DNA = DataRPKDNA %>% apply(1, function(x) mean(x[x>0], na.rm = TRUE)),
         nz.var.DNA =  DataRPKDNA %>% apply(1, function(x)  var(x[x>0], na.rm = TRUE)),
         overdispersion.DNA = nz.var.DNA/nz.mean.DNA^2,
         nz.mean.RNA = DataRPKRNA %>% apply(1, function(x) mean(x[x>0], na.rm = TRUE)),
         nz.var.RNA =  DataRPKRNA %>% apply(1, function(x)  var(x[x>0], na.rm = TRUE)),
         overdispersion.RNA = nz.var.RNA/nz.mean.RNA^2)

tmp.prop = with(mean.var, mean(nz.var.RNA/nz.mean.RNA^2 > nz.var.DNA/nz.mean.DNA^2, na.rm = TRUE))
set.seed(1)
mean.var %>% 
  na.omit %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.var.RNA/nz.mean.RNA^2, nz.var.DNA/nz.mean.DNA^2)) + geom_point() + theme_classic() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  xlab(TeX("$var (RNA| +) / \\mu (RNA| +)^2$")) +
  ylab(TeX("$var (DNA| +) / \\mu (DNA| +)^2$")) +
  annotate("text", x = 10, y = 5, label = paste0(round(tmp.prop*100, 1), "%"), col = "#F8766D") +
  annotate("text", y = 10, x = 5, label = paste0(round(100 - tmp.prop*100, 1), "%"), col = "#619CFF")
ggsave(sprintf("figure/C0101overdispersion1_%s%s.png", type, zoe.nm))

tmp.prop = with(mean.var, mean(nz.var.RNA/nz.mean.RNA > nz.var.DNA/nz.mean.DNA, na.rm = TRUE))
set.seed(1)
mean.var %>% 
  na.omit %>% 
  sample_frac(0.01) %>% 
  ggplot(aes(nz.var.RNA/nz.mean.RNA, nz.var.DNA/nz.mean.DNA)) + geom_point() + theme_classic() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  xlim(c(0, 1000)) + ylim(c(0,1000)) +
  xlab(TeX("$var (RNA| +) / \\mu (RNA| +)$")) +
  ylab(TeX("$var (DNA| +) / \\mu (DNA| +)$")) +
  annotate("text", x = 900, y = 500, label = paste0(round(tmp.prop*100, 1), "%"), col = "#F8766D") +
  annotate("text", y = 900, x = 500, label = paste0(round(100 - tmp.prop*100, 1), "%"), col = "#619CFF")
ggsave(sprintf("figure/C0101overdispersion2_%s%s.png", type, zoe.nm))

mean.var %>% 
  na.omit %>% 
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

const = ifelse(zoe == 0, 5, 1)
const2 = ifelse(zoe == 0, 2, 1)

zp.DNA <-
  zero.proportion.DNA %>%
  data.frame(x = .) %>%
  ggplot(aes(x)) + geom_histogram() + theme_classic() +
  xlab("proportion of zeros") + ylab("frequency") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_y_continuous(limits = c(0, 1.6e+5 * const)) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5,vjust=0.5)) +
  geom_vline(xintercept = 0.9, color = "red", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.95, y = .5e+5 * 1.2 * const2, col = "red",
           label = paste0(round(mean(zero.proportion.DNA >= 0.9, na.rm = TRUE) * 100, 0), "%")) +
  geom_vline(xintercept = 0.8, color = "orange", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.85, y = .5e+5 * .75 * const2, col = "orange",
           label = paste0(round(mean(zero.proportion.DNA >= 0.8, na.rm = TRUE) * 100, 0), "%")) +
  geom_vline(xintercept = 0.7, color = "cadetblue4", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.75, y = .5e+5 * .5 * const2, col = "cadetblue4",
           label = paste0(round(mean(zero.proportion.DNA >= 0.7, na.rm = TRUE) * 100, 0), "%")) +
  geom_vline(xintercept = 0.6, color = "darkslategray", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.65, y = .5e+5 * .35 * const2, col = "darkslategray",
           label = paste0(round(mean(zero.proportion.DNA >= 0.6, na.rm = TRUE) * 100, 0), "%"))
ggsave(sprintf("figure/C0101zeroproportion_DNA_%s%s.png", type, zoe.nm), zp.DNA)

zp.RNA <-
  zero.proportion %>% 
  data.frame(x = .) %>% 
  ggplot(aes(x)) + geom_histogram() + theme_classic() + 
  xlab("proportion of zeros") + ylab(" ") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  # scale_y_continuous(limits = c(0, 1.6e+5 * const)) +
  theme(axis.text.y=element_blank()) +
  geom_vline(xintercept = 0.9, color = "red", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.95, y = 1e+5 * 1.2, col = "red",
           label = paste0(round(mean(zero.proportion >= 0.9, na.rm = TRUE) * 100, 0), "%")) +
  geom_vline(xintercept = 0.8, color = "orange", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.85, y = 1e+5 * .75, col = "orange",
           label = paste0(round(mean(zero.proportion >= 0.8, na.rm = TRUE) * 100, 0), "%")) +
  geom_vline(xintercept = 0.7, color = "cadetblue4", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.75, y = 1e+5 * .5, col = "cadetblue4",
           label = paste0(round(mean(zero.proportion >= 0.7, na.rm = TRUE) * 100, 0), "%")) +
  geom_vline(xintercept = 0.6, color = "darkslategray", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.65, y = 1e+5 * .35, col = "darkslategray",
           label = paste0(round(mean(zero.proportion >= 0.6, na.rm = TRUE) * 100, 0), "%"))
ggsave(sprintf("figure/C0101zeroproportion_RNA_%s%s.png", type, zoe.nm), zp.RNA)

zp <-
  plot_grid(zp.DNA, zp.RNA + ylab(NULL), labels = "AUTO", rel_widths = c(1.12, 1), label_x = c(0.2, 0.1))
save_plot(sprintf("figure/C0101zeroproportion_%s%s.png", type, zoe.nm), plot = zp, base_width = 9, base_height = 5)

saveRDS(zp.DNA, sprintf("output/C0101zeroproportion_%s%s_DNA.rds", type, zoe.nm))
saveRDS(zp.RNA, sprintf("output/C0101zeroproportion_%s%s_RNA.rds", type, zoe.nm))




zero.proportion %>% mean(na.rm = TRUE)
zero.proportion %>% 
  cut(breaks = c(-0.1, 0.5, 0.8, 0.9, 0.95, 0.99, Inf)) %>% 
  table %>% {./length(zero.proportion)}


if (FALSE) {
    zp.DNA.1 <-   readRDS("output/C0101zeroproportion_gene_zoe1_DNA.rds")
    zp.RNA.1 <-   readRDS("output/C0101zeroproportion_gene_zoe1_RNA.rds")
    zp.DNA.2 <-   readRDS("output/C0101zeroproportion_gene_zoe2_DNA.rds")
    zp.RNA.2 <-   readRDS("output/C0101zeroproportion_gene_zoe2_RNA.rds")
    zp.DNA.IBD <- readRDS("output/C0101zeroproportion_gene_NEWDATA_DNA.rds")
    zp.RNA.IBD <- readRDS("output/C0101zeroproportion_gene_NEWDATA_RNA.rds")
    zp <-
      plot_grid(zp.DNA.2 + xlab(NULL) + ylab(" "), zp.RNA.2 + xlab(NULL) + ylab(NULL), 
                zp.DNA.1 + xlab(NULL), zp.RNA.1 + xlab(NULL) + ylab(NULL), 
                zp.DNA.IBD + ylab(" "), zp.RNA.IBD + ylab(NULL), 
                labels = "AUTO", nrow = 3, ncol = 2, rel_widths = c(1.12, 1), label_x = c(0.2, 0.1))
    save_plot("figure/C0101zeroproportion_gene_all.png", plot = zp, base_width = 9, base_height = 10)
    
}

# zp.mean = 
#   data.frame(zp.DNA = zero.proportion.DNA, mean.DNA = DataRPKDNA %>% apply(1, mean, na.rm = TRUE),
#              zp.RNA = zero.proportion, mean.RNA = DataRPKRNA %>% apply(1, mean, na.rm = TRUE))
# 
# ggplot(zp.mean) +
#   geom_point(aes(zp.DNA, mean.DNA)) +
#   scale_y_continuous(trans = "log")
# ggplot(zp.mean) +
#   geom_point(aes(zp.RNA, mean.RNA)) +
#   scale_y_continuous(trans = "log")
# ggplot(zp.mean) +
#   geom_point(aes(zp.DNA, zp.RNA))

############################################################################################
# variances of DNA and RNA
############################################################################################
# variances = data.frame(var.DNA = apply(DataRPKDNA, 1, var), var.RNA = apply(DataRPKRNA, 1, var),
#                        mean.DNA = apply(DataRPKDNA, 1, mean), mean.RNA = apply(DataRPKRNA, 1, mean))
# variances =
#   variances %>% 
#   mutate(disp.DNA = ifelse(mean.DNA < 1e-3, NA, var.DNA/mean.DNA),
#          disp.RNA = ifelse(mean.RNA < 1e-3, NA, var.RNA/mean.RNA),
#          RNAtoDNA = disp.RNA/ disp.DNA)
# variances %>% 
#   filter(var.DNA < 10000, var.RNA < 10000) %>% 
#   ggplot(aes(var.DNA, var.RNA)) + geom_point()
# variances %>% 
#   sample_frac(0.01) %>% 
#   ggplot(aes(disp.DNA, disp.RNA)) + geom_point() + xlim(c(0, 2000)) + ylim(c(0, 2000))
# variances$RNAtoDNA %>% {function(s) mean(s > 1, na.rm=TRUE)}(.)
# 
# var.nz <- function(x) var(x[x>0])
# variances.nz = data.frame(var.DNA = apply(DataRPKDNA, 1, var.nz), var.RNA = apply(DataRPKRNA, 1, var.nz))
# variances.nz %>% 
#   filter(var.DNA < 10000, var.RNA < 10000) %>% 
#   ggplot(aes(var.DNA, var.RNA)) + geom_point()
# 

