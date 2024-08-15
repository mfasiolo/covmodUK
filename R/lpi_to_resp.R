##########################
#' Convert linear predictors (on MCD scale) to covariance matrix
#' 
#' @param lpi a matrix of linear predictors
#' @param d the dimension of the response vector
#' @return a matrix with means and covariances on each row
#' @export
#' 
lpi_to_resp <- function(lpi, d){
  n <- nrow(lpi)
  mean_resp<- lpi[,1:d,drop=FALSE]
  cov_resp <- matrix(0, n, d*(d+1)/2)
  pred_resp <- matrix(0, n, d + d*(d+1)/2)
  logD2 <- matrix(0,d,d)
  T  <- matrix(0,d,d)
  Sigma <- matrix(0,d,d)
  for(i in 1:n){
    diag(logD2) <- lpi[i,(d+1):(2*d),drop=FALSE]
    diag(T) <- 1
    count <- 2*d+1
    for(j in 2:d){
      for(k in 1:(j-1)){
        T[j,k] <- lpi[i,count,drop=FALSE]
        count <- count +1
      }
    }
    Sigma <- solve(t(T)%*% diag(exp(-diag(logD2)))%*%T)
    cov_resp[i, 1:d] <- diag(Sigma)
    cov_resp[i, (d+1):(d*(d+1)/2)] <- t(Sigma)[upper.tri(t(Sigma), diag=FALSE)]
    pred_resp[i, ] <- c(mean_resp[i,], cov_resp[i,])
  }
  return(pred_resp)
}