library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
library(latex2exp)
source("F00.00.generic.R")
source("F01.01.base.R")

### 0.2 Data
# Raw data of 118 subjects
gene.marginal.RPK.DRNA <- readRDS("Data/data.geneRPK.marginal.DRNA.180.rds")
outcome <- readRDS("Data/data.outcome.180.rds") 
# gene.marginal.RPK.DRNA$otu[, , 2] is the RNA infomation

# Get only batch.DNA == 170421, and remove subject 253, 420
outcome.new = outcome[outcome$batch.DNA == 170421,]
outcome.new = outcome.new[outcome.new$id != 352 & outcome.new$id != 420,]
DataMeta116 = outcome.new

RNA = gene.marginal.RPK.DRNA$otu[, , 2]
RNA.new = RNA[,colnames(RNA) %in% DataMeta116$id]
DataRPK116 = RNA.new

rm(gene.marginal.RPK.DRNA, outcome, outcome.new, RNA, RNA.new)

# disease groups (0:Healthy, 1:Treated. 2:Diseased)
HD = which(DataMeta116$ECC %in% c(0,2))
H0 = which(DataMeta116$ECC == 0)
T1 = which(DataMeta116$ECC == 1)
D2 = which(DataMeta116$ECC == 2)

# batch
B1 = which(DataMeta116$RNA.date == "170628")
B2 = which(DataMeta116$RNA.date == "170718")

DataMeta116 <-
  DataMeta116 %>% 
  mutate(group = paste0(ifelse(cariesfree == 1, "H", "D"), ifelse(batch.RNA == "170628", 1, 2)),
         group_batch =  ifelse(batch.RNA == "170628", 1, 2),
         group_disease = ifelse(cariesfree == 1, "H", "D"))



#grp <- c("H1", "H2", "D1", "D2")
grp <- unique(DataMeta116$group)

zero.prob <- function (vec) {mean(vec == 0)}

### marginal sample (not considering batches and disease groups)
zero.proportion <- apply(DataRPK116, 1, zero.prob)
genes.regular.index <- which(zero.proportion <= 0.95) # 35%

if(F){
  set.seed(1)
  samp <- sample(genes.regular.index, 300)
  stat.ziln <- 
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, ]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi = pp,
        mu = mean(log(y.nz)),
        sig = sd(log(y.nz)))
    }) %>% t %>% as.data.frame
  
  ggplot(stat.ziln, aes(mu, sig, col = pi)) +
    geom_point() 
  
  set.seed(1)
  samp <- sample(genes.regular.index, 300)
  stat.ziln2 <- 
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, ]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi = pp,
        mu = mean(log(y.nz)),
        tht = var(log(y.nz))/mean(log(y.nz)))
    }) %>% t %>% as.data.frame
  
  ggplot(stat.ziln2, aes(mu, tht, col = pi)) +
    geom_point()  
  
  stat.ziln2 %>%
    dplyr::filter(0< tht & tht < 10) %>%
    ggplot(aes(mu, tht, col = pi)) +
    geom_point()
  
  
  ### conditional sample
  zero.proportion <- apply(DataRPK116, 1, zero.prob)
  genes.regular.index <- which(zero.proportion <= 0.95) # 35%
  
  set.seed(1)
  samp <- sample(genes.regular.index, 300)
  cond.ziln <- 
    lapply(grp, function(g) {
      sapply(samp, function(x) {
        yvec <- DataRPK116[x, DataMeta116$group == g]
        pp <- mean(yvec == 0)
        y.nz <- yvec[yvec != 0]
        
        c(pi = pp,
          mu = mean(log(y.nz)),
          sig = sd(log(y.nz)))
      }) %>% t %>% 
        as.data.frame %>% 
        mutate(group = g, disease = substr(group, 1, 1), batch = substr(group, 2, 2))
    }) 
  cond.ziln <- do.call(rbind, cond.ziln)
  
  ggplot(cond.ziln, aes(mu, sig, col = pi)) +
    geom_point() +
    facet_grid(disease ~ batch)
  cond.ziln %>%
    group_by(group) %>% 
    summarize(mu.mean = mean(mu, na.rm = TRUE),
              sig.mean = mean(sig, na.rm = TRUE),
              pi.mean = mean(pi, na.rm = TRUE))  
  
  
}
### Final para plot for main

######## p1

set.seed(1)
grp <- unique(DataMeta116$group)
samp <- sample(genes.regular.index, 300)
cond.ziln <- 
  lapply(grp, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi = pp,
        mu = mean(log(y.nz)),
        sig = sd(log(y.nz)),
        tht = var(log(y.nz))/mean(log(y.nz)))
    })%>% t %>% 
      as.data.frame %>% 
      mutate(group = g, disease = substr(group, 1, 1), batch = substr(group, 2, 2))
  }) 

cond.ziln <- do.call(rbind, cond.ziln)

cond.ziln.pi3 <- cond.ziln %>% filter (pi > 0.27 & pi < 0.33)
cond.ziln.pi5 <- cond.ziln %>% filter (pi > 0.47 & pi < 0.53)
cond.ziln.pi6 <- cond.ziln %>% filter (pi > 0.57 & pi < 0.63)
cond.ziln.pi7 <- cond.ziln %>% filter (pi > 0.67 & pi < 0.73)
cond.ziln.pi9 <- cond.ziln %>% filter (pi > 0.87 & pi < 0.93)
cond.ziln.pi3$pi_id = 0.3
cond.ziln.pi6$pi_id = 0.6
cond.ziln.pi9$pi_id = 0.9

cond.ziln.all = rbind(cond.ziln.pi3, cond.ziln.pi6, cond.ziln.pi9)
cond.ziln.all$group = substr(cond.ziln.all$group,1,2)

cond.ziln.all$pi_id_f = factor(cond.ziln.all$pi_id,
                               levels = c(0.3, 0.6, 0.9),
                               labels = c(TeX("$\\pi \\approx 0.3$"), 
                                          TeX("$\\pi \\approx 0.6$"),
                                          TeX("$\\pi \\approx 0.9$")))
# 1st col of parameter selection
cond.ziln.all %>%
  dplyr::filter(tht>0 & tht < 30) %>%
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
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p1

####### p2
set.seed(1)
samp <- sample(genes.regular.index, 300)
grp_health = c("H1","H2")
cond.ziln.delta.healthy <- 
  lapply(grp_health, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_h = pp,
        mu_h = mean(log(y.nz)),
        tht_h = var(log(y.nz))/mean(log(y.nz)))
    })%>% t %>% 
      as.data.frame %>%
      mutate( batch_h = substr(g, 2, 2))
  }) 
cond.ziln.delta.healthy <- do.call(rbind, cond.ziln.delta.healthy)

grp_diseased = c("D1", "D2")
cond.ziln.delta.diseased <- 
  lapply(grp_diseased, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_d = pp,
        mu_d = mean(log(y.nz)),
        tht_d = var(log(y.nz))/mean(log(y.nz)))
    })%>% t %>% 
      as.data.frame %>%
      mutate( batch_d = substr(g, 2, 2))
  }) 
cond.ziln.delta.diseased <- do.call(rbind, cond.ziln.delta.diseased)

cond.ziln.delta <- cbind(cond.ziln.delta.healthy, 
                         cond.ziln.delta.diseased) %>%
  mutate(delta_pi = pmax(pi_h/pi_d, pi_d/pi_h), 
         delta_mu = pmax(mu_h/mu_d, mu_d/mu_h), 
         delta_theta = pmax(tht_h/tht_d, tht_d/tht_h),
         batch = batch_h) %>%
  dplyr::select("batch", "delta_pi", "delta_mu", "delta_theta")

cond.ziln.delta.pi10 <- cond.ziln.delta %>% filter (delta_pi > 0.9 & delta_pi < 1.2)
cond.ziln.delta.pi15 <- cond.ziln.delta %>% filter (delta_pi > 1.3 & delta_pi < 1.7)
cond.ziln.delta.pi05 <- cond.ziln.delta %>% filter (delta_pi > 0.3 & delta_pi < 0.7)
cond.ziln.delta.pi10$pi_id = 1.0
cond.ziln.delta.pi15$pi_id = 1.5
cond.ziln.delta.pi05$pi_id = 0.5

cond.ziln.delta.all = rbind(cond.ziln.delta.pi05,
                            cond.ziln.delta.pi10, 
                            cond.ziln.delta.pi15)
rm(cond.ziln.delta.pi10, cond.ziln.delta.pi15, cond.ziln.delta.pi05)

cond.ziln.delta.all$pi_id_f = factor(cond.ziln.delta.all$pi_id,
                                     levels = c(0.5, 1.0, 1.5),
                                     labels = c(TeX("$\\delta_{\\pi} \\approx 0.5$"), 
                                                TeX("$\\delta_{\\pi} \\approx 1.0$"),
                                                TeX("$\\delta_{\\pi} \\approx 1.5$")))


cond.ziln.delta.all %>%
  dplyr::filter(delta_theta > 0 & delta_theta < 30 & delta_mu < 5) %>%
  ggplot(aes(delta_mu, delta_theta, col = batch, shape = batch)) +
  geom_point() +
  scale_shape_manual(name = "batch",values=c(17, 15)) +
  xlab(TeX('$\\delta_\\mu$')) + 
  ylab(TeX('$\\delta_\\theta')) + 
  ggtitle(TeX("(|$\\delta_{\\mu}$|, |$\\delta_{\\theta}$|, |$\\delta_{\\pi}$|) estimates")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p2



######## p3

set.seed(1)
samp <- sample(genes.regular.index, 300)
grp_batch1 = c("H1","D1")
cond.ziln.kappa.batch1 <- 
  lapply(grp_batch1, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_b1 = pp,
        mu_b1 = mean(log(y.nz)),
        tht_b1 = var(log(y.nz))/mean(log(y.nz)))
    })%>% t %>% 
      as.data.frame %>%
      mutate(disease1 = substr(g, 1, 1))
  }) 
cond.ziln.kappa.batch1 <- do.call(rbind, cond.ziln.kappa.batch1)

grp_batch2 = c("H2","D2")
cond.ziln.kappa.batch2 <- 
  lapply(grp_diseased, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_b2 = pp,
        mu_b2 = mean(log(y.nz)),
        tht_b2 = var(log(y.nz))/mean(log(y.nz)))
    })%>% t %>% 
      as.data.frame %>%
      mutate( disease2 = substr(g, 2, 2))
  }) 
cond.ziln.kappa.batch2 <- do.call(rbind, cond.ziln.kappa.batch2)


cond.ziln.kappa <- cbind(cond.ziln.kappa.batch1, 
                         cond.ziln.kappa.batch2) %>%
  mutate(kappa_pi = pmax(pi_b1/pi_b2, pi_b2/pi_b1), 
         kappa_mu = pmax(mu_b1/mu_b2, mu_b2/mu_b1), 
         kappa_theta = pmax(tht_b1/tht_b2, tht_b2/tht_b1),
         disease = ifelse(disease1 == "D", "dieased", "healthy")) %>%
  dplyr::select("disease", "kappa_pi", "kappa_mu", "kappa_theta")


cond.ziln.kappa.pi10 <- cond.ziln.kappa %>% filter (kappa_pi > 0.9 & kappa_pi < 1.2)
cond.ziln.kappa.pi15 <- cond.ziln.kappa %>% filter (kappa_pi > 1.3 & kappa_pi < 1.7)
cond.ziln.kappa.pi20 <- cond.ziln.kappa %>% filter (kappa_pi > 1.8 & kappa_pi < 2.2)
cond.ziln.kappa.pi10$pi_id = 1.0
cond.ziln.kappa.pi15$pi_id = 1.5
cond.ziln.kappa.pi20$pi_id = 2.0

cond.ziln.kappa.all = rbind(cond.ziln.kappa.pi10, cond.ziln.kappa.pi15, cond.ziln.kappa.pi20)

cond.ziln.kappa.all$pi_id_f = factor(cond.ziln.kappa.all$pi_id,
                                     levels = c(1.0, 1.5, 2.0),
                                     labels = c(TeX("$\\kappa_{\\pi} \\approx 1.0$"), 
                                                TeX("$\\kappa_{\\pi} \\approx 1.5$"),
                                                TeX("$\\kappa_{\\pi} \\approx 2.0$")))






cond.ziln.kappa.all %>%
  dplyr::filter(kappa_theta > 0 & kappa_theta < 30 & kappa_mu < 5) %>%
  ggplot(aes(kappa_mu, kappa_theta, col = disease, shape = disease)) +
  geom_point() +
  scale_shape_manual(name = "disease",values=c(16,1)) +
  xlab(TeX('$\\kappa_\\mu$')) + 
  ylab(TeX('$\\kappa_\\theta')) + 
  ggtitle(TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p3

grid.arrange(
  grobs = list(p1, p2, p3),
  ncol = 3, nrow = 1
) -> p

ggsave(("Document/plot/para_selection_ziln.png"),  p,  width = 10, height = 10)#, dpi = 200)
#dev.off()

cond.ziln.pi3 %>%
  group_by(group) %>% 
  summarize(mu.mean = mean(mu, na.rm = TRUE),
            sig.mean = mean(sig, na.rm = TRUE),
            pi.mean = mean(pi, na.rm = TRUE))  
