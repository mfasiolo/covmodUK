#' RD of GPD distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#'
gpd_rd <- function(mu, wt, scale)  {
  
  return( bundle_GPD$rd(mu, wt, scale)  )
  
}