#' QF of shash distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @export
#'
#' @examples
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