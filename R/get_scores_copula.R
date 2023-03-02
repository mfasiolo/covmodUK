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
get_scores_copula <- function(x, cop_param, d, nsim, ncores, A){
  
  out <- lapply(x, function(dat){
    test <- dat$test
    train <- dat$train 
    modpred <- dat[1:14]
    n <- nrow(dat$test)
    scores <- matrix(NA, n, 4)
    
    fam_nam <- modpred[[1]]$model$family$family
    if(fam_nam == "shash"){
      my_pdf <- shash_pdf
      cdf <- modpred[[1]]$model$family$cdf
      qf <- function(p, mu, wt, scale, logp) {
        mu <- as.matrix(mu)
        if(ncol(mu)==1){ mu <- t(mu) }
        muE <- mu[ , 1, drop = TRUE]
        sigE <- exp(mu[ , 2, drop = TRUE])
        epsE <- mu[ , 3, drop = TRUE]
        delE <- exp(mu[ , 4, drop = TRUE])
        
        q <- muE + (delE * sigE) * sinh((1/delE) * asinh(qnorm(p, log.p = logp)) + (epsE/delE))
        
        return( q)
      }
    } else {
      if(fam_nam == "gaulss"){
        my_pdf <- function(y, eta, log) {
          dnorm(y, mean = eta[ , 1], sd = 1/eta[ , 2], log = log)
        }
        cdf <- function(q, mu, logp = FALSE) {
          pnorm(q, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp)
        }
        qf  <- function(p, mu, logp = FALSE) {
          qnorm(p, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp)
        }
      } else {
        stop("unknown family")
      }
    }
    
    # For each marginal: get z by inverting shash CDF and compute log-likelihood
    my_list <- lapply(modpred, function(o){
    
      y_nam <- as.character(o$model$formula[[1]][2])
      y_nam <- substr(y_nam,nchar(y_nam)-1,nchar(y_nam))
      ytrain <- o$model$y
      ytest <- test[ , paste0("res_", y_nam)]
      ztrain <- cdf(ytrain, o$model$fitted.values, logp = TRUE)
      ztest <- cdf(ytest, o$eta_hat, logp = TRUE)
      
      return(list("ytest" = ytest, "ztrain" = ztrain, "ztest" = ztest, "dtest" = my_pdf(ytest, o$eta_hat, log = TRUE), 
                  "eta_train" = o$model$fitted.values, "eta_test" = o$eta_hat))
    })
    
    if( is.null(cop_param) ){ # Static copula case: empirical correlation of z-transformed residuals used
      S <- cor(sapply(my_list, function(x) qnorm(x$ztrain, log.p = TRUE)))
      cop <- normalCopula(P2p(S), dim = 14, dispstr = "un")
      scores[ , 1] <- -dCopula(sapply(my_list, function(x) pmin(exp(x$ztest), 1-1e-8)), cop, log = TRUE) - rowSums(sapply(my_list, "[[", "dtest"))
      
      y <- sapply(my_list, "[[", "ytest")
      sig_hat <- chol(S)
      eta_train <-  lapply(my_list, "[[", "eta_train")
      eta_test <-  lapply(my_list, "[[", "eta_test")
      for(ii in 1:n){
        yii <- y[ii, , drop = TRUE]
        traj_emp <- mvnfast::rmvn(n = nsim, mu = rep(0, 14), sigma = sig_hat, isChol = TRUE)
        for(kk in 1:14){
          mu_ii_kk <- matrix(eta_test[[kk]][ii, ], ncol = ncol(eta_train[[1]]), nrow = nsim, byrow = TRUE) 
          traj_emp[ , kk] <- qf(pnorm(traj_emp[ , kk], log.p = TRUE), mu = mu_ii_kk, logp = TRUE)
        }
        scores[ii, 2] <- es_sample(y = yii, dat = t(traj_emp))
        scores[ii, 3] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 0.5)
        scores[ii, 4] <- vs_sample(y = yii, dat = t(traj_emp), w = NULL, p = 1)
      }
      
    } else { # Dynamical copula case: we fitting the MCD model the z-transformed residuals
      
    }
    
    return(scores)
  })
  
  return(out)
}
