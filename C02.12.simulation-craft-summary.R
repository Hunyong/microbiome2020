### 0. library
  source("F02.01.metrics.R")
  source("F02.12.metrics-craft.R")
  library(abind)
  library(ggplot2)

### 1. Settings, naming, etc
  study.nm = c("ZOE 2.0 pilot", "ZOE 2.0", "IBD")
  
  methods.levels = c(paste0("LB",  c(".nonz", ".zero", ".glob", ".min")), 
                     paste0("MAST",  c(".nonz", ".zero", ".glob", ".min")),
                     paste0("Wg",  c(".nonz", ".zero", ".glob", ".min")),
                     "DS2", "DS2ZI", "MGS", "ANCOM.sz", "ANCOM", 
                     "LN", "KW", "KW-II", "WRS", "LFE", "ALDEX")
  methods.labels = 
    methods.levels %>% 
    gsub("\\.glob", "", .) %>% 
    gsub("^LFE$", "LEfSe", .)
  methods.to.include = c("LN", "LB", "MAST", "DS2", "DS2ZI", "MGS", "ANCOM", "KW", "KW-II", "LEfSe", "ALDEX")

  metric.levels = c("type1error", "FDR", "sensitivity", "accuracy", "AUC")
  metric.labels = c("Type 1 Error", "FDR", "Sensitivity", "Accuracy", "AUC")
  
### 2. Collection
  result = NULL
  for (zoe in 1:3) {
    for(n.gene in c(1e+4, 1e+5)) {
      for(n.signal in c(100, 300, 1000)) {
        for (BH in c(TRUE, FALSE)) {
          result.tmp =
            tab.metrics(zoe = zoe, type = "gene", n.signal = n.signal, n.gene = n.gene, BH.correction = BH) 
          
          if (!is.null(result.tmp)) {
            result.tmp = result.tmp %>% 
              as_tibble(rownames = "method") %>% 
              tidyr::gather(value = "value", key = "metric", - method) %>% 
              mutate(n.signal = n.signal, n.gene = n.gene, zoe = zoe) %>% 
              {if (!BH) {filter(., metric != "FDR")} else filter(., metric == "FDR")}
            result = rbind(result, result.tmp)
          }
        }
      }
    }
  }
  
### 3. plots
  result %>% 
    filter(n.gene == 10000) %>% 
    filter(n.signal %in% c(100, 1000)) %>%
    mutate(method = factor(method, levels = methods.levels, labels = methods.labels),
           study = factor(study.nm[zoe], levels = study.nm),
           metric = factor(metric, levels = metric.levels, labels = metric.labels)) %>% 
    filter(method %in% methods.to.include) %>% 
    filter(metric != "AUC") %>% 
    ggplot(aes(method, value, col = method, group = method)) +
    geom_point() + guides(col = FALSE) +
    facet_grid(metric ~ paste0("# signal genes = ", n.signal) + study, scales = "free_y") +
    geom_hline(data = data.frame(metric = "type1error"), col = "red", yintercept = 0.05, mapping = "metric") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsve()


