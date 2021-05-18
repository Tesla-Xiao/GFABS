hvcspr = function(y, x, u, status = NULL, bs.df = 5, bs.degree = 3, eps = 0.01,
                    sigma=NULL, lambda.min=NULL, model="spr", criterion="bic")  
{
  basis <- bs(u, df = bs.df, degree = bs.degree)
  phi   <- cbind(1, basis)
  W     <- as.matrix(t(KhatriRao(t(x), t(phi))))
  
  n     <- nrow(x)
  p     <- ncol(x)
  group <- sort(c(2*(1:p)-1, rep(2*(1:p), each = bs.df)))
  
  if (is.null(status))         status = rep(1, n)
  if (is.null(sigma))           sigma = 1/sqrt(n)
  if (is.null(lambda.min)) lambda.min = {if (n > ncol(W)) 1e-4 else .02}
  
  if (criterion == "bic") 
    fit <- bicgrpFabs(W, y, group, status, sigma, NULL, model, "L2", TRUE, TRUE, 
                      eps, 10^-6, 10^4, lambda.min)

  # columns of Beta are the estimated beta_j(u)'s.
  Beta  <- matrix(0, n, p)
  q     <- ncol(phi)
  opt   <- fit$opt
  for (i in 1:p) 
    Beta[, i] <- fit$beta[((i-1)*q+1):(i*q), opt] %*% t(phi)
  
  list(theta     = fit$beta, 
       Beta      = Beta,
       lambda    = fit$lambda,
       BIC       = fit$BIC,
       optimal   = opt,
       direction = fit$direction, 
       active    = fit$active,
       iter      = fit$iter)
}



