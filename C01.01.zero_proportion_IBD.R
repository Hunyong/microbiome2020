library(dplyr); library(ggplot2); library(cowplot)
dat.new.DNA = readRDS("../MicrobiomePaper2020/Nature2019data/data.ecs_relab.geneProp.joint.DNA.rds")$otu
dat.new.RNA = readRDS("../MicrobiomePaper2020/Nature2019data/data.ecs_relab.geneProp.joint.RNA.rds")$otu

zero.proportion =
  dat.new.RNA %>% 
  apply(1, function(x) mean(x == 0)) 
zero.proportion %>% mean # avg zero proportion = 0.933

zero.proportion.DNA =
  dat.new.DNA %>% 
  apply(1, function(x) mean(x == 0)) 
zero.proportion.DNA %>% mean # avg zero proportion = 0.906

zp.DNA <-
  zero.proportion.DNA %>% 
  data.frame(x = .) %>% 
  ggplot(aes(x)) + geom_histogram() + theme_classic() + 
  xlab("proportion of zeros") + ylab("frequency") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_y_continuous(limits = c(0, 0.5e+5)) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5,vjust=0.5)) +
  geom_vline(xintercept = 0.9, color = "red", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.95, y = .3e+5 * 1.2, col = "red",
           label = paste0(round(mean(zero.proportion.DNA >= 0.9) * 100, 0), "%")) +
  geom_vline(xintercept = 0.8, color = "orange", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.85, y = .2e+5 * .75, col = "orange",
           label = paste0(round(mean(zero.proportion.DNA >= 0.8) * 100, 0), "%")) +
  geom_vline(xintercept = 0.7, color = "cadetblue4", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.75, y = .2e+5 * .5, col = "cadetblue4",
           label = paste0(round(mean(zero.proportion.DNA >= 0.7) * 100, 0), "%")) +
  geom_vline(xintercept = 0.6, color = "darkslategray", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.65, y = .2e+5 * .35, col = "darkslategray",
           label = paste0(round(mean(zero.proportion.DNA >= 0.6) * 100, 0), "%"))
ggsave(zp.DNA, "figure/C0101zeroproportionDNA_IBD.png")

zp.RNA <-
  zero.proportion %>% 
  data.frame(x = .) %>% 
  ggplot(aes(x)) + geom_histogram() + theme_classic() + 
  xlab("proportion of zeros") + ylab(" ") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_y_continuous(limits = c(0, 0.5e+5)) +
  theme(axis.text.y=element_blank()) +
  geom_vline(xintercept = 0.9, color = "red", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.95, y = 3e+4 * 1.2, col = "red",
           label = paste0(round(mean(zero.proportion >= 0.9) * 100, 0), "%")) +
  geom_vline(xintercept = 0.8, color = "orange", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.85, y = 1e+4 * .75, col = "orange",
           label = paste0(round(mean(zero.proportion >= 0.8) * 100, 0), "%")) +
  geom_vline(xintercept = 0.7, color = "cadetblue4", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.75, y = 1e+4 * .5, col = "cadetblue4",
           label = paste0(round(mean(zero.proportion >= 0.7) * 100, 0), "%")) +
  geom_vline(xintercept = 0.6, color = "darkslategray", size = 0.2, linetype = "dashed") +
  annotate(geom = "text", x =  0.65, y = 1e+4 * .35, col = "darkslategray",
           label = paste0(round(mean(zero.proportion >= 0.6) * 100, 0), "%"))
ggsave(zp.RNA, "figure/C0101zeroproportion_IBD.png")
zp <-
  plot_grid(zp.DNA, zp.RNA + ylab(NULL), labels = "AUTO", rel_widths = c(1.12, 1), label_x = c(0.2, 0.1))
save_plot("figure/C0101zeroproportion_IBD.png", plot = zp, base_width = 9, base_height = 5)


dat.new.RNA %>% apply(1, function(x) any(x>0, na.rm = T)) %>% sum
