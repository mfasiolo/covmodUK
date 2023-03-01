##########################
#' Get performance scores
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @export
#' @importFrom Matrix rankMatrix
#' @details This function is meant for internal use only.
#' 
get_scores <- function(x, d, nsim, ncores, A){
  
  out <- mclapply(x, function(dat){
    y <- as.matrix(dat$y)
    mu <- dat$mu_hat
    n <- nrow(y)
    defic <- rankMatrix(A) < nrow(A)
    scores <- matrix(NA, n, 4)
    for(ii in 1:n){
      yii <- drop(A %*% y[ii, , drop = TRUE])
      mean_hat <- A %*% mu[ii, 1:d, drop = TRUE]
      sig_hat <- chol(A %*% Sigma_mat(mu[ii, -(1:d)]) %*% t(A), pivot = TRUE)
      traj_emp <- mvnfast::rmvn(n=nsim, mu = mean_hat[attr(sig_hat, "pivot")], 
                                sigma = sig_hat, isChol = TRUE)[ , order(attr(sig_hat, "pivot"))]
      if(!defic){
        scores[ii, 1] <- - mvnfast::dmvn(X = yii[attr(sig_hat, "pivot")], mu = mean_hat[attr(sig_hat, "pivot")], 
                                         sigma = sig_hat, log=TRUE, isChol = TRUE)
      }
      scores[ii, 2] <- es_sample(y = yii, dat = t(traj_emp))
      scores[ii, 3] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 0.5)
      scores[ii, 4] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 1)
    }
    return( scores )
  }, mc.cores = ncores)
  
  return(out)
  
}

