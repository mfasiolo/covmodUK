#' Quantile function of the shash distribution
#'
#' @param p the probability level in (0, 1)
#' @param mu a matrix of shash parameters
#' @param wt a vector of weights
#' @param scale the scale parameter (ignored here)
#' @param logp if TRUE then p is log(p)
#'
#' @return a vector of quantiles
#' @export
#'
shash_qf <- function(p, mu, wt, scale, logp) {
  mu <- as.matrix(mu)
  if(ncol(mu)==1){ mu <- t(mu) }
  muE <- mu[ , 1, drop = TRUE]
  sigE <- exp(mu[ , 2, drop = TRUE])
  epsE <- mu[ , 3, drop = TRUE]
  delE <- exp(mu[ , 4, drop = TRUE])
  
  q <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(p, log.p = logp)) + (epsE/delE))
  
  return(q)
}