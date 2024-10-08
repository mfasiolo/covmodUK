#' The pdf of the GPD distribution
#'
#' @param y the data on which the pdf is evaluated
#' @param eta matrix of linear predictors
#' @param log if TRUE the log-pdf is returned
#' 
#' @return the pdf corresponding to each y
#' @export
#'
gpd_pdf <- function(y, eta, log){
  if( is.vector(eta) ) { eta <- matrix(eta, nrow = 1) }
  phi <- eta[ , 1, drop = TRUE]
  xi <- eta[ , 2, drop = TRUE]
  tol <- 1e-7
  if( any(abs(xi) < tol) ){ 
    xi[xi < tol & xi > 0] <- tol
    xi[xi > -tol & xi <= 0] <- - tol
  }
  Cxi <- xi * (1 + xi) * y
  l <- log1p(xi) - log(phi) - (1/xi+1) * log1p(Cxi / phi)
  return(l)
}