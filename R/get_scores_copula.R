##########################
#' Get performance scores
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @export
#' @details This function is meant for internal use only.
#' 
get_scores_copula <- function(fit_lss, nsim, ncores, A, indep = FALSE){
  
  d <- length(fit_lss[[1]])
  
  fam <- attr(fit_lss, "nam")
  qf_fun <- get(paste0(fam, "_qf"))
  pdf_fun <- get(paste0(fam, "_pdf"))
  
  d_A <- nrow(A)
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("qf_fun", "fam", "pdf_fun", "indep", "d", "d_A", "nsim", "A", "d_gau_copula"), 
                envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  my_fun_aofhs <- function(my_list){
    
    n <- length(my_list[[1]]$resid_test)
    scores <- matrix(0, n, 9)
    
    # Static copula: empirical correlation of z-transformed residuals used
    if( !indep ){
      S <- cor(sapply(my_list, function(x) qnorm(x$resid_train, log.p = TRUE)))
    } else {
      S <- diag(1, 14)
    }
    
    y <- sapply(my_list, "[[", "y_test")
    mu_train <-  lapply(my_list, "[[", "mu_train")
    mu_test <-  lapply(my_list, "[[", "mu_test")
    
    if(d_A == d){
      scores[ , 1] <- - rowSums(sapply(my_list, "[[", "llk_test")) 
      scores[ , 2] <- - d_gau_copula(sapply(my_list, function(x) x$resid_test), R = t(chol(S)), log = TRUE) 
      my_q_fun <- function(p, mu, traj){
        qf_fun(p, mu, logp = FALSE)
      }
    } else {
      my_q_fun <- function(p, mu, traj){
        quantile(traj, p)
      }
      if(fam == "gaulss"){
        for(ii in 1:n){
          yii <- drop(A %*% y[ii, , drop = TRUE])
          mu_A <- A %*% sapply(mu_test, function(x) x[ii, 1])
          sd_mar <- 1 / sapply(mu_test, function(x) x[ii, 2])
          Sig_A <- A %*% (diag(sd_mar, ncol=d) %*% S %*% diag(sd_mar, ncol=d)) %*% t(A)
          scores[ii, 1] <- - sum(dnorm(yii, mu_A, sqrt(diag(Sig_A)), log = TRUE))
          scores[ii, 2] <- - (mvnfast::dmvn(yii, mu_A, Sig_A, log = TRUE) + scores[ii, 1])
        }
      }
    }
    
    sig_hat <- chol(S)

    for(ii in 1:n){
      yii <- drop(A %*% y[ii, , drop = TRUE])
      traj_emp <- mvnfast::rmvn(n = nsim, mu = rep(0, 14), sigma = sig_hat, isChol = TRUE)
      for(kk in 1:14){
        mu_ii_kk <- matrix(mu_test[[kk]][ii, ], ncol = ncol(mu_train[[1]]), nrow = nsim, byrow = TRUE) 
        traj_emp[ , kk] <- qf_fun(pnorm(traj_emp[ , kk], log.p = TRUE), mu = mu_ii_kk, logp = TRUE)
        #scores[ii, 10] <- scores[ii, 10] + pdf_fun(yii[kk], mu_ii_kk[1, , drop = FALSE], log = TRUE)
      }
      traj_emp <- t(A %*% t(traj_emp)) 
      scores[ii, 3] <- sum(sapply(1:d_A, function(kk) crps_sample(y = yii[kk], dat = traj_emp[ , kk])))
      scores[ii, 4] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = my_q_fun(0.001, mu_test[[kk]][ii, , drop = FALSE], traj_emp[ , kk]), qu = 0.001)))
      scores[ii, 5] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = my_q_fun(0.01, mu_test[[kk]][ii, , drop = FALSE], traj_emp[ , kk]), qu = 0.01)))
      scores[ii, 6] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = my_q_fun(0.99, mu_test[[kk]][ii, , drop = FALSE], traj_emp[ , kk]), qu = 0.99)))
      scores[ii, 7] <- sum(sapply(1:d_A, function(kk) pinLoss(y = yii[kk], mu = my_q_fun(0.999, mu_test[[kk]][ii, , drop = FALSE], traj_emp[ , kk]), qu = 0.999)))
      scores[ii, 8] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 0.5)
      scores[ii, 9] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 1)
    } 
    return(scores)
  }
  environment(my_fun_aofhs) <- .GlobalEnv
  
  out <- parLapply(NULL, fit_lss, my_fun_aofhs)
  
  stopCluster(cl)
  rm(cl)
  
  return(out)
}
