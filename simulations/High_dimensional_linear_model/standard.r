standard = function(W, group) {
  p = ncol(W)
  n = nrow(W)
  center = colMeans(W)
  W.mean = W - matrix(rep(center, n), n, p, byrow=T)
  scale = sqrt(colSums(W.mean^2)/n)
  xx = t(t(W.mean)/scale)
  list(xx = xx, center = center, scale = scale)
}

unstandard = function(b, center, scale, group)
{
  b/scale
}