#' A Group Forward and Backward Stagewise (GFabs) algorithm for smoothed partial rank loss function (SPR) with the Group Lasso penalty.
#'
#' @param W The design matrix.
#' @param y The survival outcome.
#' @param group The grouping vector.
#' @param status The censoring indicator.
#' @param sigma The smoothing parameter in SPR.
#' @param weight The weight vector of groups.
#' @param model The loss function used.
#' @param type The norm of penalty used.
#' @param back The indicator of whether to take backward steps.
#' @param stoping The indicator of whether to stop iteration when lambda is less than lambda.min.
#' @param eps The step size for GFabs.
#' @param xi The threshhold for GFabs.
#' @param iter The maximum number of outer-loop iterations allowed.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max.
#'
#' @return A list.
#' \itemize{
#'   \item beta - estimation of covariates
#'   \item lambda - lambda sequence
#'   \item direction - direction of GFabs
#'   \item active - active sets
#'   \item iter - iterations
#'   \item BIC - bic criteria
#'   \item group - The grouping vector.
#'   \item opt - position of the optimal tuning based on BIC.
#' }
#' @export
#'
#' @examples
#' W = matrix(rnorm(80), 20, 4)
#' b = c(5, 5, -5, -5)
#' e = c(0.7*rnorm(20)+0.3*rcauchy(20))
#' y = W %*% b + e
#' group <- c(1, 1, 2, 2)
#' fit <- GFabs(W, y, group)

GFabs = function(W, y, group, status=NULL, sigma=NULL, weight=NULL,
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

  fit <- .Call(GFabs:::"BIC_grpFabs",
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

  # unstandardization
  K1    = cumsum(tabulate(group))
  K0    = c(1, K1[-length(K1)]+1)
  beta = matrix(0, nrow=nrow(beta), ncol=ncol(beta))
  for (i in 1:length(K0)) {
    idx = c(K0[i]:K1[i])
    beta[idx,] = scale[[i]] %*% beta[idx,]
  }

  val = list(beta      = beta,
             lambda    = fit$lambda,
             direction = fit$direction,
             active    = fit$active,
             iter      = fit$iter,
             group     = group,
             BIC       = fit$bic,
             opt       = fit$opt)
  return(val)
}




