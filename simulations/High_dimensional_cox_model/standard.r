standard = function(W, group) 
{
  p = ncol(W)
  n = nrow(W)
  K1    = cumsum(tabulate(group))
  K0    = c(1, K1[-length(K1)]+1)
  
  center = colMeans(W)
  W.mean = W - matrix(rep(center, n), n, p, byrow=T)
  
  scale <- vector("list", length(K0))
  xx = matrix(0, n, p)
  for (i in 1:length(K0)) {
    idx = c(K0[i]:K1[i])
    xx[,idx] = W[,idx]
    scale[[i]] = solve(chol(t(xx[,idx]) %*% xx[,idx]/n))
    xx[,idx] = xx[,idx] %*% scale[[i]]
  }
  list(xx = xx, center = center, scale = scale)
}

unstandard = function(b, center, scale, group)
{
  K1    = cumsum(tabulate(group))
  K0    = c(1, K1[-length(K1)]+1)
  beta = matrix(0, nrow=nrow(b), ncol=ncol(b))
  for (i in 1:length(K0)) {
    idx = c(K0[i]:K1[i])
    beta[idx,] = scale[[i]] %*% b[idx,]
  }
  beta
}