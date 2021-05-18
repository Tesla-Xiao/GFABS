bicgrpFabs = function(W, y, group, status=NULL, sigma=NULL, weight=NULL, 
                      model="spr", type = "L2", back=TRUE, stoping=TRUE, 
                      eps = 0.01, xi = 10^-6, iter=10^4, lambda.min = NULL)
{
  ## Reorder groups, if necessary
  gf    = as.factor(group)
  g     =  as.numeric(gf)
  G     = max(g)
  g.ord = order(g)
  g     = g[g.ord]
  W     = W[,g.ord]
  K     = as.numeric(table(g))
  K1    = cumsum(K)
  K0    = c(1, K1[-G]+1)
  
  ## standard
  std    = standard(W, group)
  W.std  = std[[1]]
  center = std[[2]]
  scale  = std[[3]]
  y.std  = y - mean(y)
  
  n = nrow(W.std)
  p = ncol(W.std)
  pg = length(K0)
  param = c(n, p, pg, 0)
  if (is.null(status))         status = rep(1, n)
  if (is.null(sigma))           sigma = 1/sqrt(n)
  if (is.null(lambda.min)) lambda.min = {if (n > p) 1e-4 else .02}
  if (is.null(weight))         weight = sqrt(K)
  
  if (model == "lm") model = 1
  if (model == "spr") model = 2
  if (type == "L2") type = 2
  
  fit <- .Call("BIC_grpFabs",
               as.numeric(t(W.std)),
               as.numeric(y.std),
               as.numeric(weight),
               as.integer(status),
               as.integer(K0-1),
               as.integer(K1-1),
               as.integer(model),
               as.numeric(sigma),
               as.numeric(eps),
               as.numeric(lambda.min),
               as.numeric(xi),
               as.integer(type),
               as.integer(back),
               as.integer(stoping),
               as.integer(iter),
               as.integer(param) )
  
  beta = matrix(fit$beta, nrow = p)
  beta = unstandard(beta, center, scale, group) 
  
  val = list(beta      = beta,
             lambda    = fit$lambda,
             direction = fit$direction,
             active    = fit$active,
             iter      = fit$iter,
             group     = group,
             BIC       = fit$bic,
             opt       = fit$opt)
}




