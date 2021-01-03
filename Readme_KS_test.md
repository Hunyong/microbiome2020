### the Kruskal-Wallis test (or the Lilliefors procedure)  
This can be implemented by the following example code.  
The code includes testing the Beta, log-normal, and gamma distributions. 
However, other distributions can be included by tailoring the code.  
The Kruskal-Wallis test works for continuous distributions.  
The distinction between the Kruskall-Wallis test and the Lilliefors procedure is that
the former tests the goodness of fit of data to a distribution with a fixed set of parameters,
while the latter tests that with an estimated distribution which is often the case in practice.  

'@ data A vector of the random variable.
'@ model The model against which the goodness of fit is measured. Either "beta", "gamma", or "ln" for Beta, Gamma, or log-normal, respectively.
'@ ks.pval p-value of the Kruskal-Wallis test. This is lower than the p-value of the Lilliefors procedure and should not be used when the model parameters are estimated.
'@ lilliefors Whether the Lilliefors procedure is used.
'@ n.lilliefors The number of replicates used to calculate the p-value of the Lilliefors procedure.

# example - Beta distribution

library(gamlss) # for beta distribution
source("F01.01.goodness_of_fit.R") # ks.empirical() function
set.seed(1)
dat = rbeta(100, shape1 = 1, shape2 = 1)
ks.empirical(dat, model = "beta", return.est = FALSE, ks.pval = FALSE, lilliefors = TRUE, n.lilliefors = 300)
# p-val = 0.787 # Failed to reject the null (null = lack of fit)
ks.empirical(dat, model = "gamma", return.est = FALSE, ks.pval = FALSE, lilliefors = TRUE, n.lilliefors = 300)
# p-val = 0.013 # The null is rejected at 5% significance level.
ks.empirical(dat, model = "ln", return.est = FALSE, ks.pval = FALSE, lilliefors = TRUE, n.lilliefors = 300)
# p-val = 0.000 # The null is rejected at 5% significance level.

