#' RD of shash distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#'
shash_rd <- function(mu, wt, scale)  {
  
  return( shash()$rd(mu, wt, scale)  )

}