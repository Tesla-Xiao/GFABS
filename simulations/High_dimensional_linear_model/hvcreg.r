hvcreg = function(y, x, u, bs.df, bs.degree, 
                  lambda.min=NULL, model="ols", criterion="bic") 
{
  basis <- bs(u, df =  bs.df, degree = bs.degree)
  phi   <- cbind(1, basis)
  W     <- as.matrix(t(KhatriRao(t(x), t(phi))))
  
  n     <- nrow(x)
  p     <- ncol(x)	
  group <- sort(c(2*(1:p)-1, rep(2*(1:p), each = bs.df)))
  if(is.null(lambda.min)) lambda.min = {if (n > ncol(W)) 1e-4 else .02}
  
  if (model == "ols"){
      if (criterion == "cv") {
      fit    <- cv.grpreg(W, y, group, lambda.min=lambda.min)
      theta  <- fit$fit$beta[-1, ]
      lambda <- fit$lambda
      opt    <- fit$min
      cve	   <- fit$cve
    } else {
      fit    <- grpreg::grpreg(W, y, group, lambda.min=lambda.min)
      theta  <- fit$beta[-1, ]
      lambda <- fit$lambda
      opt    <- which(lambda == select(fit)$lambda)
      cve    <- select(fit)$IC
    }	
  } 
  
  # columns of Beta are the estimated beta_j(u)'s.
  Beta  <- matrix(0, ncol=p, nrow=n)
  q     <- ncol(phi)
  for (i in 1:p)  
    Beta[, i] <- theta[((i-1)*q+1):(i*q), opt] %*% t(phi)
  
  list(theta     = theta, 
       Beta      = Beta,
       lambda    = lambda,
       BIC       = cve,
       optimal   = opt)
}



