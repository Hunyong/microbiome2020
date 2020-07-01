x = (1:1000)/1000 # uniform numbers
x.a = asn(x) * pi/2 # arcsin and some normalization
x.s = sqrt(x)
plot(x, x.s - x.a)
plot(x, x.a, type = "l")
lines(x, x.s, col = "red")

RNA = readRDS("../Data-Processed/data.geneRPK.marginal.DRNA.ZOE1.rds")$otu[,, 2]
ST  = apply(RNA, 2, sum)
mean(ST) # 5,551,718 (ZOE1), 21M for ZOE2
tpm <- t(t(RNA)/ST)


# tpm
tpm1 <- tpm[19, ] ;   sum(tpm1 > 0) # 79
full.ln    <- gamlss(tpm1[tpm1 > 0] ~ 1, family = LOGNO(sigma.link = "log"), 
                     control = gamlss.control(n.cyc = 100, trace = FALSE))
a3 = param.ln(full.ln)
tmp.ln <- ks.test(tpm1[tpm1 > 0], "pLOGNO", mu = a3["mu"], sigma = a3["sig"])
tmp.ln$statistic %>% as.numeric
tmp.ln$p.value %>% as.numeric

# arcsin
asn1 <- asn(tpm1)
full.ln    <- gamlss(asn1[asn1 > 0] ~ 1, family = LOGNO(sigma.link = "log"), 
                     control = gamlss.control(n.cyc = 100, trace = FALSE))
a3 = param.ln(full.ln)
tmp.ln <- ks.test(asn1[asn1 > 0], "pLOGNO", mu = a3["mu"], sigma = a3["sig"])
tmp.ln$statistic %>% as.numeric
tmp.ln$p.value %>% as.numeric
