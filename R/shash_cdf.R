#' CDF of shash distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#'
shash_cdf <- function(q, mu, wt, scale, logp) {

  return( shash()$cdf(q = q, mu = mu, wt = wt, scale = scale, logp = logp) )
  
}