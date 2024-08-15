##########################
#' Density of Gaussian copula
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @export
#' @noRd
#'
d_gau_copula <- function(logu, R, log = TRUE){
  logu <- as.matrix(logu)
  n <- nrow(logu)
  out <- numeric(n)
  for(ii in 1:n){
    z <- qnorm(logu[ii, , drop = TRUE], log.p = TRUE)
    logdetR <- 2 * sum(log(diag(R)))
    Riz <- forwardsolve(R, z)
    out[ii] <- - logdetR / 2 - (crossprod(Riz) - crossprod(z))/2
  }
  return(out)
}