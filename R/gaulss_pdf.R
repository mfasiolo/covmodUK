#' PDF of gaussian GAMLSS
#' 
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#' 
gaulss_pdf <- function(y, eta, log) {
  dnorm(y, mean = eta[ , 1], sd = 1/eta[ , 2], log = log)
}