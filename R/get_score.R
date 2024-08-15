##########################
#' Get performance scores for MVN GAMs on test set
#' 
#' @param x a matrix of predictions from the MVN GAM
#' @param d the dimension of the response vector
#' @param nsim the number of simulation use to compute (e.g.) the CRPS loss
#' @param ncores the number of cores
#' @param A a matrix used to specify the linear combination of regions we are looking at
#' @return A list of performance scores.
#' @export
#' 
get_scores <- function(x, d, nsim, ncores, A){
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("x", "d", "nsim", "A"), 
                envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  my_fun_afbias <- function(dat){
    y <- as.matrix(dat$y)
    mu <- dat$mu_hat
    n <- nrow(y)
    #defic <- rankMatrix(A) < nrow(A)
    #if(defic) stop("A is rank deficient")
    scores <- matrix(NA, n, 9)
    for(ii in 1:n){
      yii <- drop(A %*% y[ii, , drop = TRUE])
      mean_hat <- A %*% mu[ii, 1:d, drop = TRUE]
      sig_hat <- A %*% Sigma_mat(mu[ii, -(1:d)]) %*% t(A)
      sig_chol <- chol(sig_hat)
      traj_emp <- mvnfast::rmvn(n=nsim, mu = mean_hat, 
                                sigma = sig_chol, isChol = TRUE)
      d_A <- nrow(A)
      #if(!defic){
        scores[ii, 1] <- - mvnfast::dmvn(X = yii, mu = mean_hat, sigma = diag(diag(sig_hat), ncol = nrow(A)), log=TRUE)
        scores[ii, 2] <- - mvnfast::dmvn(X = yii, mu = mean_hat, sigma = sig_chol, log=TRUE, isChol = TRUE) - scores[ii, 1]
        scores[ii, 3] <- sum(sapply(1:d_A, function(kk) crps_sample(y = yii[kk], dat = traj_emp[ , kk])))
        scores[ii, 4] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = qnorm(0.001, 0, sqrt(sig_hat[kk, kk])), qu = 0.001)))
        scores[ii, 5] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = qnorm(0.01, 0, sqrt(sig_hat[kk, kk])), qu = 0.01)))
        scores[ii, 6] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = qnorm(0.99, 0, sqrt(sig_hat[kk, kk])), qu = 0.99)))
        scores[ii, 7] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = qnorm(0.999, 0, sqrt(sig_hat[kk, kk])), qu = 0.999)))
      #}
      scores[ii, 8] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 0.5)
      scores[ii, 9] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 1)
    }

    return( scores )
  }
  environment(my_fun_afbias) <- .GlobalEnv
  
  out <- parLapply(NULL, x, my_fun_afbias)
  
  stopCluster(cl)
  rm(cl)
  
  return(out)
  
}

