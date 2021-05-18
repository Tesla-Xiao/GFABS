"
Created on Thur Nov 1 20:27:06 2018

simulations for varing coefficient spr+group fabs.

consider right censored response data in low dimension 
under the accelerated failure time model.
"

library(mvtnorm)
library(survival)
library(Matrix)
library(splines) 
library(grpreg)
source("standard.r") 
source("bic_grpFabs.r")
source("hvcspr.r")
source("hvcsurv.r")
#dyn.load("grpFabs.dll")
dyn.load("grpFabs.so")

# generate data =============================================================
assig = function(n_args) {
  # example: n_args = c(2, 3, 4)
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
  # expand.grid: Create a data frame from all combinations of the supplied vectors or factors
  t(expand.grid(cargs))
}


myAUC <- function(XBeta, y)
{
  AUC.k = cor( y, XBeta, method = "kendall")
  
  AUC.k
}

generator = function(n, p, rho, error, tran, censor.rate = NULL) {
  sigma = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y)))
  
  d=5
  x = rmvnorm(n, mean = rep(0,p), sigma = sigma)
  u = runif(n)
  b = cbind(5*sin(2*pi*u), 5*cos(2*pi*u), 5*exp(2*u-1)-5, 5, -5, matrix(0, nrow=n, ncol=p-d))
  
  # linear predictor
  eta = numeric(n)
  for (i in 1:n) eta[i] = x[i, 1:d] %*% b[i, 1:d]
  
  # error   ==================================================
  if(error == "norm") e = rnorm(n, 0, 1)
  #if(error =="contaminate") e = c(rnorm(0.7*n), rcauchy(0.3*n))
  if(error =="contaminate") e = c(0.7*rnorm(n)+0.3*rcauchy(n))
  if(substr(error, 1, 1) == "t") e = rt(n, as.numeric(substr(error,2,2)))
  # the natural log of a Weibull random time is an extreme value random observation.   
  if(error == "ev") e = log(rweibull(n, shape = 1, scale = 1)) - digamma(1)
  
  # transformation ===========================================
  #if(tran == "x^3")  g = function(x) x^3
  #if(tran == "lm") g = function(x) x
  if(tran == "log") g = function(x) exp(x)
  
  y0 = g(eta + e)
  
  # censor  ==================================================
  if (is.null(censor.rate) != TRUE){
    cens = quantile(y0, 1-censor.rate)
    y = pmin(y0, cens)
    status = 1 * (y0<=cens)
  }else{
    y = y0
    status = rep(1, n)
  }
  
  # In the transformation model, y may be infinity.
  y[which(y=="NaN"|y=="Inf")] =  max(y[which(y!="NaN"&y!="Inf")]) 
  list(x=x, u=u, y=y, beta=b, status=status)                        
}


main <- function(number) {
  number      = as.numeric(number)-1
  B           = 500
  tran        = "log"
  bs.df       = 7
  bs.degree   = 3
  eps         = c(1/10, 1/50, 1/100)
  model       = "cox"
  censor.rate = 0.2
  n           = 300
  p           = 500
  
  n_args      = c(B, 2, 3)
  jobs        = assig(n_args)
  
  neach = 20
  output <- vector("list", neach)
  for (loops in 1:neach) {
    njobs       = number*neach+loops
    id          = jobs[, njobs]
    rplc        = (1:B)[id[1]]
    rho         = c(0.3, 0.7)[id[2]]
    error       = c("ev","contaminate", "t2")[id[3]]
    set.seed(njobs*10000)
    
    # generate data ======================================================
    dat    = generator(n,  p,  rho, error, tran, censor.rate)
    x      = dat$x
    y      = dat$y
    status = dat$status
    u      = dat$u
    bet    = dat$beta
    #b.std  = bet/sqrt(sum(bet^2))
    
    testing_dat = generator(n,  p,  rho, error, tran, censor.rate)
    testing_basis <- bs(testing_dat$u, df = bs.df, degree = bs.degree)
    testing_phi   <- cbind(1, testing_basis)
    testing_W     <- as.matrix(t(KhatriRao(t(testing_dat$x), t(testing_phi))))
    
    # estimate ===========================================================
    # spr
    n_eps     <- length(eps)
    fit <- vector("list", n_eps+1)
    for (i in 1:n_eps) {
      fit[[i]] <- hvcspr(y, x, u, status, bs.df, bs.degree, eps[i], NULL, NULL, "spr", "bic")
      fit[[i]]$XBeta <- rowSums(x * fit[[i]]$Beta)
      fit[[i]]$AUC = myAUC(fit[[i]]$XBeta, y)
      fit[[i]]$testing_XBeta <- rowSums(testing_W %*% fit[[i]]$theta[,fit[[i]]$optimal])
      fit[[i]]$testing_AUC = myAUC(fit[[i]]$testing_XBeta, testing_dat$y)
    }
    # cox
    i = i+1
    fit[[i]] <- grpreg(y, x, u, status, bs.df, bs.degree)
    fit[[i]]$XBeta <- rowSums(x * fit[[i]]$Beta)
    fit[[i]]$AUC = myAUC(fit[[i]]$XBeta, y)
    fit[[i]]$testing_XBeta <- rowSums(testing_W %*% fit[[i]]$theta[,fit[[i]]$optimal])
    fit[[i]]$testing_AUC = myAUC(fit[[i]]$testing_XBeta, testing_dat$y)
    
    
    print(loops)
    
    # save result ========================================================
    output[[loops]] <- list(fit=fit, seed = njobs*10000)
  }
  saveRDS(output, file = paste(tran, n, p, rho, error, rplc-neach+1, rplc, ".rds", sep="-"))
}

args <- commandArgs(TRUE)
main(args[1])



