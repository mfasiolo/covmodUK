##########################
#' Get performance scores from shash+gpd+copula model
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @export
#' @details This function is meant for internal use only.
#' 
get_scores_copula_gpd <- function(fit_shash, fit_gpd, sets, nsim, ncores, A, dyn_copula = NULL, indep = FALSE){
  
  d <- length(fit_shash[[1]])
  d_A <- nrow(A)
  wrapper <- postprocess_gpd_fits(fit_shash, fit_gpd)
  resid <- t(sapply(1:tail(sets, 1), wrapper$resid))
  llk <- t(sapply(1:tail(sets, 1), wrapper$llk))
  y <- sapply(1:d, function(ii) c(fit_shash[[1]][[ii]]$y_train, do.call("c", lapply(fit_shash, function(o) o[[ii]]$y_test))))
  
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("sets", "d", "d_A", "y", "llk", "wrapper", "resid", "llk", "nsim", "A", "d_gau_copula", "indep", "dyn_copula"), 
                envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  if( !is.null(dyn_copula) ){
    mu_cop <- matrix(NA, nrow(y), d * (d + 1) / 2)
    tmp <- do.call("rbind", lapply(dyn_copula, "[[", "mu_hat"))[ , -(1:14)]
    mu_cop[(sets[1]+1):tail(sets, 1), ] <- tmp
    clusterExport(NULL, c("mu_cop"), envir = environment())
  }
  
  my_fun_fafaf <- function(ii){
    train <- 1:sets[ii]
    test <- (sets[ii]+1):sets[ii+1]
    out <- matrix(NA, length(test), 9)
    
    # Static copula: empirical correlation of z-transformed residuals used
    if( !indep ){
      S <- cor(qnorm(resid[train, ], log.p = TRUE))
    } else {
      S <- diag(1, 14)
    }
    
    if(d_A == d){
      my_q_fun <- function(index, p, traj, ll){
        wrapper$qf(index = kk, p = p, logp = FALSE)[ll]
      }
    } else {
      my_q_fun <- function(index, p, traj, ll){
        quantile(traj[ll, ], p)
      }
    }
    
    ytest <- y[test, ]
    sig_hat <- chol(S)
    jj <- 1
    for(kk in test){
      if(!is.null(dyn_copula)){
        S <- cov2cor(Sigma_mat(mu_cop[kk, ]))
        sig_hat <- chol(S)
      }
      ykk <- drop(A %*% y[kk, , drop = TRUE])
      psim <- apply(mvnfast::rmvn(n = nsim, mu = rep(0, 14), sigma = sig_hat, isChol = TRUE), 2, pnorm, log.p = TRUE)
      traj_emp <- A %*% apply(psim, 1, function(p) wrapper$qf(index = kk, p = p, logp = TRUE))
      if(d_A == d){
      out[jj, 1] <- - sum(llk[kk, ])  
      out[jj, 2] <- - d_gau_copula(resid[kk, , drop = FALSE], R = t(sig_hat), log = TRUE)
      }
      out[jj, 3] <- sum(sapply(1:d_A, function(ll) crps_sample(y = ykk[ll], dat = traj_emp[ll, ])))
      out[jj, 4] <- sum(sapply(1:d_A, function(ll) pinLoss(y = ykk[ll], mu = my_q_fun(index = kk, p = 0.001, traj_emp, ll), qu = 0.001)))
      out[jj, 5] <- sum(sapply(1:d_A, function(ll) pinLoss(y = ykk[ll], mu = my_q_fun(index = kk, p = 0.01, traj_emp, ll), qu = 0.01)))
      out[jj, 6] <- sum(sapply(1:d_A, function(ll) pinLoss(y = ykk[ll], mu = my_q_fun(index = kk, p = 0.99, traj_emp, ll), qu = 0.99)))
      out[jj, 7] <- sum(sapply(1:d_A, function(ll) pinLoss(y = ykk[ll], mu = my_q_fun(index = kk, p = 0.999, traj_emp, ll), qu = 0.999)))
      out[jj, 8] <- vs_sample(y = ykk, dat = traj_emp, w = NULL, p = 0.5)
      out[jj, 9] <- vs_sample(y = ykk, dat = traj_emp, w = NULL, p = 1)
      jj <- jj + 1
    }
    return(out)
  }
  environment(my_fun_fafaf) <- .GlobalEnv
  
  scores <- parLapply(NULL, 1:(length(sets)-1), my_fun_fafaf)
  
  stopCluster(cl)
  rm(cl)
  
  return(scores)
}
