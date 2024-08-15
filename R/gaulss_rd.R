#' RD of gaulss distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#'
gaulss_rd <- function(mu, wt, scale)  {
  
  return( gaulss()$rd(mu, wt, scale)  )

}