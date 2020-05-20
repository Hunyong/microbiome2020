
##### 1. parameter estimators
### 1.1 parameter estimator for ZINB ###
  library(pscl)
  
  ZINB.ML <- function(yvec, xvec = 1, formula = y ~ x, notation="mtp") {
    #notation (parametrization):
    # abp: alpha-beta-pi
    # mtp: mu-theta-pi
    # * linkage:: alpha = theta, beta = theta/mu
    if (!notation %in% c("abp", "mtp")) {stop("Notation is not correctly specified.")}
    
    d <- xvec %>% unique %>% length
    if (d==1) {formula = y~1} else if (d>2) {
      warning("factor x is forced to be numeric")
      d = 2
    }
    data = data.frame(y=as.numeric(yvec) %>% round, x = as.numeric(xvec))
    a <- pscl::zeroinfl(formula = formula, data = data, dist = "negbin", EM = TRUE)
    theta <- a$theta # theta
    bg <- a$coef # regression (beta, gamma)
    
    pp <- plogis(cumsum(bg[[2]]))
    mu.nonzero = exp(cumsum(bg[[1]])) # e^beta
    beta = theta/mu.nonzero
    
    if (notation=="abp") {
      result <- c(theta, beta, pp)
      names(result) = c("alpha", paste0("beta",0:(d-1)), paste0("pi",0:(d-1)))
    } else if (notation == "mtp") {
      result <- c(mu.nonzero, theta, pp)
      names(result) = c(paste0("mu",0:(d-1)), "theta", paste0("pi",0:(d-1)))
    }
    
    return(result)
  }
  
  
  timeout <- function(expression, timeout = 3, exception.out = rep(NA,3), error.out = rep(NaN,3)) {
    tryCatch(R.utils::withTimeout(expression, timeout = timeout), 
             TimeoutException = function(ex) {message("Timeout. Skipping.")
                                              return(exception.out)},
              error = function(ex) {message("Estimation error. Skipping.")
                                     return(error.out)}
             ) -> a
    return(a)
  }
  
  ZINB.ML.time <- function(yvec, xvec = 1, formula = y ~ x, notation="mtp",...) {
    timeout(expression = ZINB.ML(yvec, xvec = xvec, formula = formula, notation=notation), ...) 
  }
  
### 1.2 random number generator for ZINB ###
  rZINB <- function(n, param=NULL, m, t, p) {
    if (!is.null(param)) {m = param[1]; t = param[2]; p = param[3]}

    # gamma
    # rvec <- rgamma(n, shape = a, rate = 1/b)
    rvec <- rgamma(n, shape = m/t, rate = 1/t) #because (a = m/t) and (b = t)
    
    # poisson from gamma
    rvec <- rpois(n, rvec)
    
    # pi
    pvec <- rbinom(n, 1, 1-p)
    xvec <- rvec * pvec
    return(xvec)
  }
  
  rZINB.abp <- function(n, param=NULL, a, b, p) {
    if (!is.null(param)) {a = param[1]; b = param[2]; p = param[3]}
    
    # gamma
    rvec <- rgamma(n, shape = a, rate = 1/b)
    # poisson from gamma
    rvec <- rpois(n, rvec)
    
    # pi
    pvec <- rbinom(n, 1, 1-p)
    xvec <- rvec * pvec
    return(xvec)
  }
  
  which.mode <- function(x) {
    x <- x[complete.cases(x)]
    d <- density(x)
    return(d$x[which.max(d$y)])
  }
  
  # rZINB.sim is in F02.01.simulation.R
  
if (FALSE) { #example
  
  i=100
  tmp = data.frame(y=as.numeric(DataRPK116[i,-1]) %>% round, x = 1)
  ZINB.ML(tmp$y, tmp$x, notation="abp")
  ZINB.ML(tmp$y, tmp$x, notation="mtp")
  
  tmp = data.frame(y=as.numeric(DataRPK116[i,-1]) %>% round, x = DataMeta116$ECC)
  ZINB.ML(tmp$y, tmp$x)
  
  
  
  par(mfrow=c(1,2))
  as.numeric(gene.marginal.RPK.RNA[4,-1])%>% round %>% hist
  rZINB(100, 3.68, .7644, .236) %>% hist  # similar
  par(mfrow=c(1,1))
  
  gene.marginal.RPK.RNA[1,-1] %>% as.numeric %>% round -> a  # normal
  gene.marginal.RPK.RNA[5,-1] %>% as.numeric %>% round -> a  # time out
  gene.marginal.RPK.RNA[6,-1] %>% as.numeric %>% round -> a  # error
  
  timeout(ZINB.ML(a))
  
}
  
  
### 2.1 random number generator for ZI-Gamma ###
  rZIG <- function(n, param=NULL, m, t, p) {
    if (!is.null(param)) {m = param[1]; t = param[2]; p = param[3]}
    
    # gamma
    rvec <- rgamma(n, shape = m/t, rate = 1/t) #because (a = m/t) and (b = t)
    
    # pi
    pvec <- rbinom(n, 1, 1-p)
    xvec <- rvec * pvec
    return(xvec)
  }
  
  if (FALSE) {
    library(dplyr)
    rZIG(100, param = c(m= 1, t = 1, p = 0.5)) %>% {c(mean(.==0), mean(.[.>0]))}
  }
  
### 2.2 random number generator for ZI-Log-normal ###
  rZILN <- function(n, param=NULL, m, sig, p) {
    if (!is.null(param)) {m = param[1]; t = param[2]; p = param[3]}
    
    sig2 = log (t/m + 1)
    sig = sqrt(sig2)
    mm = log (m) - sig2/2
    
    # normal
    rvec <- exp(rnorm(n, mean = mm, sd = sig)) #because (a = m/t) and (b = t)
    
    # pi
    pvec <- rbinom(n, 1, 1 - p)
    xvec <- rvec * pvec
    return(xvec)
  }
  
  if (FALSE) {
    library(dplyr)
    #rZILN(100, param = c(m = 1, sig = 1, p = 0.5)) %>% {c(mean(.==0), mean(log(.[.>0])))}
    rZILN(100, param = c(m = 1, t = 1, p = 0.5)) %>% {c(mean(.==0), mean(log(.[.>0])))}
  }