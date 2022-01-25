source("F02.01.metrics.R")
source("F02.12.metrics-craft.R")
library(abind)

# a$ranks.TP

result = NULL
for (zoe in 1:3) {
  for(n.gene in c(1e+4, 1e+5)) {
    for(n.signal in c(100, 300, 1000)) {
      result.tmp =
        tab.metrics(zoe = zoe, type = "gene", n.signal = n.signal, n.gene = n.gene) 
      if (!is.null(result.tmp)) {
        result.tmp = result.tmp %>% 
          as_tibble(rownames = "method") %>% 
          tidyr::gather(value = "value", key = "metric", - method) %>% 
          mutate(n.signal = n.signal, n.gene = n.gene, zoe = zoe)
        result = rbind(result, result.tmp)
      }
    }
  }
}
result %>% 
  filter(n.signal == 300) %>% 
  ggplot(aes(method, value, col = method, group = method)) +
  geom_point() +
  facet_grid(metric ~ factor(n.signal)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))



