#' CDF of gaussian GAMLSS
#' 
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @export
#'
#' @examples
gaulss_cdf <- function(q, mu, wt, scale, logp = FALSE) {
  pnorm(q, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp)
}