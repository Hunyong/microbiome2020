library(dplyr); library(magrittr); library(ggplot2); library(gridExtra)
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

zero.prob <- function (vec) {mean(vec == 0)}

# marginal sample (not considering batches and disease groups)
zero.proportion <- apply(DataRPK116, 1, zero.prob)
genes.regular.index <- which(zero.proportion <= 0.95) # 35%

set.seed(1)
samp <- sample(genes.regular.index, 300)
stat.zigamma <- 
  sapply(samp, function(x) {
    tmp <- DataRPK116[x, ]
    tmp.nz <- tmp[tmp>0]
    c(pi = zero.prob(tmp),
      mu = mean(tmp),             # = alpha * beta
      theta = var(tmp)/mean(tmp)) # = alpha * beta^2
  })
stat.zigamma <- as.data.frame(t(stat.zigamma))

ggplot(stat.zigamma, aes(mu, theta, col = pi)) +
  geom_point()


### Final para plot for main

######## p1
grp <- unique(DataMeta116$group)
set.seed(1)
samp <- sample(genes.regular.index, 300)
cond.zig <- 
  lapply(grp, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi = pp,
        mu = mean(y.nz),
        tht = var(y.nz)/mean(y.nz))
    })%>% t %>% 
      as.data.frame %>% 
      mutate(group = g, disease = substr(group, 1, 1), batch = substr(group, 2, 2))
  }) 

cond.zig <- do.call(rbind, cond.zig)

cond.zig.pi3 <- cond.zig %>% filter (pi > 0.27 & pi < 0.33)
cond.zig.pi5 <- cond.zig %>% filter (pi > 0.47 & pi < 0.53)
cond.zig.pi6 <- cond.zig %>% filter (pi > 0.57 & pi < 0.63)
cond.zig.pi7 <- cond.zig %>% filter (pi > 0.67 & pi < 0.73)
cond.zig.pi9 <- cond.zig %>% filter (pi > 0.87 & pi < 0.93)
cond.zig.pi3$pi_id = 0.3
cond.zig.pi6$pi_id = 0.6
cond.zig.pi9$pi_id = 0.9

cond.zig.all = rbind(cond.zig.pi3, cond.zig.pi6, cond.zig.pi9)
cond.zig.all$group = substr(cond.zig.all$group,1,2)

cond.zig.all$pi_id_f = factor(cond.zig.all$pi_id,
                               levels = c(0.3, 0.6, 0.9),
                               labels = c(TeX("$\\pi \\approx 0.3$"), 
                                          TeX("$\\pi \\approx 0.6$"),
                                          TeX("$\\pi \\approx 0.9$")))
# 1st col of parameter selection
cond.zig.all %>%
  #dplyr::filter(tht>0 & tht < 30 & mu < 10) %>%
  dplyr::filter(tht>0 & tht < 8 & mu < 7) %>%
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
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p_zig1

####### p2
set.seed(1)
samp <- sample(genes.regular.index, 300)
grp_health = c("H1","H2")
cond.zig.delta.healthy <- 
  lapply(grp_health, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_h = pp,
        mu_h = mean(y.nz),
        tht_h = var(y.nz)/mean(y.nz))
    })%>% t %>% 
      as.data.frame %>%
      mutate( batch_h = substr(g, 2, 2))
  }) 
cond.zig.delta.healthy <- do.call(rbind, cond.zig.delta.healthy)

grp_diseased = c("D1", "D2")
cond.zig.delta.diseased <- 
  lapply(grp_diseased, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_d = pp,
        mu_d = mean(y.nz),
        tht_d = var(y.nz)/mean(y.nz))
    })%>% t %>% 
      as.data.frame %>%
      mutate( batch_d = substr(g, 2, 2))
  }) 
cond.zig.delta.diseased <- do.call(rbind, cond.zig.delta.diseased)

cond.zig.delta <- cbind(cond.zig.delta.healthy, 
                        cond.zig.delta.diseased) %>%
  mutate(delta_pi = pmax(pi_h/pi_d, pi_d/pi_h), 
         delta_mu = pmax(mu_h/mu_d, mu_d/mu_h), 
         delta_theta = pmax(tht_h/tht_d, tht_d/tht_h),
         batch = batch_h) %>%
  dplyr::select("batch", "delta_pi", "delta_mu", "delta_theta")

cond.zig.delta.pi10 <- cond.zig.delta %>% filter (delta_pi > 0.9 & delta_pi < 1.2)
cond.zig.delta.pi15 <- cond.zig.delta %>% filter (delta_pi > 1.3 & delta_pi < 1.7)
cond.zig.delta.pi20 <- cond.zig.delta %>% filter (delta_pi > 1.8 & delta_pi < 2.2)
cond.zig.delta.pi10$pi_id = 1.0
cond.zig.delta.pi15$pi_id = 1.5
cond.zig.delta.pi20$pi_id = 2.0

cond.zig.delta.all = rbind(cond.zig.delta.pi10, cond.zig.delta.pi15, cond.zig.delta.pi20)

cond.zig.delta.all$pi_id_f = factor(cond.zig.delta.all$pi_id,
                                     levels = c(1.0, 1.5, 2.0),
                                     labels = c(TeX("$\\delta_{\\pi} \\approx 1.0$"), 
                                                TeX("$\\delta_{\\pi} \\approx 1.5$"),
                                                TeX("$\\delta_{\\pi} \\approx 2.0$")))


cond.zig.delta.all %>%
  #dplyr::filter(delta_theta > 0 & delta_theta < 10 & delta_mu < 10) %>%
  dplyr::filter(delta_theta > 0 & delta_theta < 30 & delta_mu < 4) %>%
  ggplot(aes(delta_mu, delta_theta, col = batch, shape = batch)) +
  geom_point() +
  scale_shape_manual(name = "batch",values=c(17, 15)) +
  xlab(TeX('$\\delta_\\mu$')) + 
  ylab(TeX('$\\delta_\\theta')) + 
  ggtitle(TeX("(|$\\delta_{\\mu}$|, |$\\delta_{\\theta}$|, |$\\delta_{\\pi}$|) estimates")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p_zig2



######## p3

set.seed(1)
samp <- sample(genes.regular.index, 300)
grp_batch1 = c("H1","D1")
cond.zig.kappa.batch1 <- 
  lapply(grp_batch1, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_b1 = pp,
        mu_b1 = mean(y.nz),
        tht_b1 = var(y.nz)/mean(y.nz))
    })%>% t %>% 
      as.data.frame %>%
      mutate(disease1 = substr(g, 1, 1))
  }) 
cond.zig.kappa.batch1 <- do.call(rbind, cond.zig.kappa.batch1)

grp_batch2 = c("H2","D2")
cond.zig.kappa.batch2 <- 
  lapply(grp_diseased, function(g) {
    sapply(samp, function(x) {
      yvec <- DataRPK116[x, DataMeta116$group == g]
      pp <- mean(yvec == 0)
      y.nz <- yvec[yvec != 0]
      
      c(pi_b2 = pp,
        mu_b2 = mean(y.nz),
        tht_b2 = var(y.nz)/mean(y.nz))
    })%>% t %>% 
      as.data.frame %>%
      mutate( disease2 = substr(g, 2, 2))
  }) 
cond.zig.kappa.batch2 <- do.call(rbind, cond.zig.kappa.batch2)


cond.zig.kappa <- cbind(cond.zig.kappa.batch1, 
                         cond.zig.kappa.batch2) %>%
  mutate(kappa_pi = pmax(pi_b1/pi_b2, pi_b2/pi_b1), 
         kappa_mu = pmax(mu_b1/mu_b2, mu_b2/mu_b1), 
         kappa_theta = pmax(tht_b1/tht_b2, tht_b2/tht_b1),
         disease = ifelse(disease1 == "D", "dieased", "healthy")) %>%
  dplyr::select("disease", "kappa_pi", "kappa_mu", "kappa_theta")


cond.zig.kappa.pi10 <- cond.zig.kappa %>% filter (kappa_pi > 0.9 & kappa_pi < 1.2)
cond.zig.kappa.pi15 <- cond.zig.kappa %>% filter (kappa_pi > 1.3 & kappa_pi < 1.7)
cond.zig.kappa.pi20 <- cond.zig.kappa %>% filter (kappa_pi > 1.8 & kappa_pi < 2.2)
cond.zig.kappa.pi10$pi_id = 1.0
cond.zig.kappa.pi15$pi_id = 1.5
cond.zig.kappa.pi20$pi_id = 2.0

cond.zig.kappa.all = rbind(cond.zig.kappa.pi10, cond.zig.kappa.pi15, cond.zig.kappa.pi20)

cond.zig.kappa.all$pi_id_f = factor(cond.zig.kappa.all$pi_id,
                                     levels = c(1.0, 1.5, 2.0),
                                     labels = c(TeX("$\\kappa_{\\pi} \\approx 1.0$"), 
                                                TeX("$\\kappa_{\\pi} \\approx 1.5$"),
                                                TeX("$\\kappa_{\\pi} \\approx 2.0$")))






cond.zig.kappa.all %>%
  #dplyr::filter(kappa_theta > 0 & kappa_mu < 10 & kappa_theta < 10) %>%
  dplyr::filter(kappa_theta > 0 & kappa_mu < 30 & kappa_theta < 3.5) %>%
  ggplot(aes(kappa_mu, kappa_theta, col = disease, shape = disease)) +
  geom_point() +
  scale_shape_manual(name = "disease",values=c(16,1)) +
  xlab(TeX('$\\kappa_\\mu$')) + 
  ylab(TeX('$\\kappa_\\theta')) + 
  ggtitle(TeX("(|$\\kappa_{\\mu}$|, |$\\kappa_{\\theta}$|, |$\\kappa_{\\pi}$|) estimates")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  facet_grid(rows = vars(pi_id_f), labeller = label_parsed) -> p_zig3

grid.arrange(
  grobs = list(p_zig1, p_zig2, p_zig3),
  ncol = 3, nrow = 1
) -> p_zig

ggsave(("Document/plot/para_selection_zig.png"),  p_zig,  width = 10, height = 10)#, dpi = 200)
#dev.off()


