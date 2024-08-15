#' pdf of shash distribution
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @noRd
#' @export
#'
shash_pdf <- function(y, eta, log){
  .log1pexp <- function(x)
  {
    indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), T)
    kk <- which(indx==1)
    if( length(kk) ){  x[kk] <- exp(x[kk])  }
    kk <- which(indx==2)
    if( length(kk) ){  x[kk] <- log1p( exp(x[kk]) ) }
    kk <- which(indx==3)
    if( length(kk) ){  x[kk] <- x[kk] + exp(-x[kk]) }
    return(x)
  }
  mu <- eta[ , 1]
  sig <- exp(eta[ , 2])
  eps <- eta[ , 3]
  del <- exp(eta[ , 4])
  tau <- log(sig)
  phi <- log(del)
  z <- (y - mu) / (sig*del)
  dTasMe <- del*asinh(z) - eps
  g <- -dTasMe
  CC <- cosh( dTasMe )
  SS <- sinh( dTasMe )
  l <- - tau - 0.5*log(2*pi) + log(CC) - 0.5*.log1pexp(2*log(abs(z))) - 0.5*SS^2
  if( !log ){ l <- exp(l) }
  return(l)
}
