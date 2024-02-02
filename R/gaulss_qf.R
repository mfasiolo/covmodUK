#' QF of gaussian GAMLSS
#' 
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @export
#'
#' @examples
gaulss_qf <- function(p, mu, logp = FALSE) {
  qnorm(p, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp)
}