"
Created on Thur Nov 1 20:27:06 2018

simulations for varing coefficient spr+group fabs.

consider right censored response data in low dimension 
under the accelerated failure time model.
"

library(mvtnorm)
library(Matrix)
library(splines) 
library(GFabs) 
source("hvcreg.r")
library(grpreg)

# generate data =============================================================
assig = function(n_args) {
  # example: n_args = c(2, 3, 4)
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
  # expand.grid: Create a data frame from all combinations of the supplied vectors or factors
  t(expand.grid(cargs))
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
  if(error =="contaminate") e = c(0.7*rnorm(n)+0.3*rcauchy(n))
  if(substr(error, 1, 1) == "t") e = rt(n, as.numeric(substr(error,2,2)))
  
  # transformation ===========================================
  if(tran == "lm") g = function(x) x
  
  y0 = g(eta + e)
  list(x=x, u=u, y=y0, beta=b)                        
}

main <- function(number) {
  number      = as.numeric(number)-1
  B           = 500
  tran        = "lm"
  bs.df       = 7
  bs.degree   = 3
  eps         = c(1/10, 1/50, 1/100)
  model       = "ols"
  n           = 200
  p           = 500
  
  n_args      = c(B, 2, 3)
  jobs        = assig(n_args)
  
  output <- vector("list", 25)
  for (loops in 1:25) {
    njobs       = number*25+loops
    id          = jobs[, njobs]
    rplc        = (1:B)[id[1]]
    rho         = c(0.3, 0.7)[id[2]]
    error       = c("norm", "contaminate", "t2")[id[3]]
    set.seed(njobs*10000)
    
    # generate data ======================================================
    dat    = generator(n,  p,  rho, error, tran)
    x      = dat$x
    y      = dat$y
    u      = dat$u
    bet    = dat$beta
    #b.std  = bet/sqrt(sum(bet^2))
    
    testing_dat = generator(n,  p,  rho, error, tran)
    testing_basis <- bs(testing_dat$u, df = bs.df, degree = bs.degree)
    testing_phi   <- cbind(1, testing_basis)
    testing_W     <- as.matrix(t(KhatriRao(t(testing_dat$x), t(testing_phi))))
    
    # estimate ===========================================================
    # spr
    n_eps     <- length(eps)
    fit <- vector("list", n_eps+1)
    for (i in 1:n_eps) {
      fit[[i]] <- GFabs_vc(x, y, u, bs.df=bs.df, bs.degree=bs.degree, model="spr", eps=eps[i])
      fit[[i]]$XBeta <- rowSums(x * fit[[i]]$Beta)
      fit[[i]]$AUC = myAUC(fit[[i]]$XBeta, y)
      fit[[i]]$testing_XBeta <- rowSums(testing_W %*% fit[[i]]$theta[,fit[[i]]$optimal])
      fit[[i]]$testing_AUC = cor( testing_dat$y, fit[[i]]$testing_XBeta, method = "kendall")
    }
    # lm
    i = i+1
    fit[[i]] <- hvcreg(y, x, u, bs.df, bs.degree)
    fit[[i]]$XBeta <- rowSums(x * fit[[i]]$Beta)
    fit[[i]]$AUC = myAUC(fit[[i]]$XBeta, y)
    fit[[i]]$testing_XBeta <- rowSums(testing_W %*% fit[[i]]$theta[,fit[[i]]$optimal])
    fit[[i]]$testing_AUC = cor( testing_dat$y, fit[[i]]$testing_XBeta, method = "kendall")
    
    
    print(loops)
    
    # save result ========================================================
    output[[loops]] <- list(fit=fit, seed = njobs*10000)
  }
  saveRDS(output, file = paste(tran, n, p, rho, error, rplc-24, rplc, ".rds", sep="-"))
}

args <- commandArgs(TRUE)
main(args[1])



