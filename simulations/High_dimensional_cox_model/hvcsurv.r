grpreg <- function(y, x, u, status, bs.df, bs.degree, lambda.min = NULL) {
  basis <- bs(u, df = bs.df, degree = bs.degree)
  phi   <- cbind(1, basis)
  W     <- as.matrix(t(KhatriRao(t(x), t(phi))))
  
  n     <- nrow(x)
  p     <- ncol(x)
  group <- sort(c(2*(1:p)-1, rep(2*(1:p), each = bs.df)))
  if (missing(status)) status = rep(1, n)
  if(is.null(lambda.min)) lambda.min = {if (n > ncol(W)) 1e-4 else .02}
  
  fit <- grpreg::grpsurv(W, Surv(y, status), group, lambda.min=lambda.min)
  lambda <- fit$lambda
  df <- fit$df[1:length(lambda)]
  loss <- fit$loss
  bic = -loss/(n*(1-df/n)^2)
  opt = which.max(bic)
  
  beta <- -fit$beta[,1:length(lambda)]
  
  Beta <- matrix(0, n, p)
  q     <- ncol(phi)
  for (i in 1:p) {
      Beta[, i] <- phi %*% beta[((i-1)*q+1):(i*q), opt, drop=F]
  }
  
  
  val <- list(lambda     = lambda,
              theta      = beta,
              loss       = loss,
              df         = df,
              BIC        = bic,
              optimal    = opt,
              Beta       = Beta)
  return(val)
}