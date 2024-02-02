#' CDF of shash distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @export
#'
#' @examples
shash_cdf <- function(q, mu, wt, scale, logp) {

  return( shash()$cdf(q = q, mu = mu, wt = wt, scale = scale, logp = logp) )
  
}