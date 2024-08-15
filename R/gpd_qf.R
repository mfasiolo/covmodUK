#' QF of GPD GAMLSS
#' 
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#'
gpd_qf <- function(p, mu, logp) {
  mu <- as.matrix(mu)
  if( ncol(mu) == 1 ){ mu <- t(mu) }
  xi <- mu[ , 2, drop = TRUE]
  sig <- mu[ , 1, drop = TRUE] / (1 + xi)
  if(logp){
    p <- exp(p)
  }
  q <- sig * ((1-p)^(-xi) - 1) / xi
  return( q )
}