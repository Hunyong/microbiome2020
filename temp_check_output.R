# stat-n400-pert0.5-zig-9.5.3.rds
library(dplyr)

ptrn = "tmp_([0-9]*)_([0-9]*)_([0-9]*)_([0-9]*)_pert(.*)_n([0-9]*)_s([0-9]*)\\.txt"
tmp.files =
  list.files(".", pattern = "tmp_[0-9]+") %>% 
  data.frame(a = .) %>% 
  transmute(i = gsub(ptrn, "\\1", a) %>% as.numeric,
            j = gsub(ptrn, "\\2", a) %>% as.numeric,
            k = gsub(ptrn, "\\3", a) %>% as.numeric,
            model = gsub(ptrn, "\\4", a) %>% as.numeric %>% {c("zinb", "zig", "ziln")[.]},
            pert = gsub(ptrn, "\\5", a) %>% as.numeric %>% "/"(10),
            n = gsub(ptrn, "\\6", a) %>% as.numeric,
            s = gsub(ptrn, "\\7", a) %>% as.numeric)
tmp.files %T>%
  print %>% 
  write.csv("tmp_bookkeeping.csv")


ptrn2 = ".*stat-n([0-9]*)-pert(.*)-([a-z]*)-([0-9]*)\\.([0-9]*)\\.([0-9]*)\\.rds"
list.files("output", pattern = "stat.*\\.rds") %>% 
  data.frame(a = .) %>% 
  transmute(i = gsub(ptrn2, "\\4", a) %>% as.numeric,
            j = gsub(ptrn2, "\\5", a) %>% as.numeric,
            k = gsub(ptrn2, "\\6", a) %>% as.numeric,
            model = gsub(ptrn2, "\\3", a),
            pert = gsub(ptrn2, "\\2", a) %>% as.numeric,
            n = gsub(ptrn2, "\\1", a) %>% as.numeric) %>% 
  group_by(j, k, model, pert, n) %>% 
  summarize(n.i = n(), n.short = 10 - n.i) %>% 
  filter(n.short > 0) %>% 
  as.data.frame
# 