##########################
#' Get normalized residuals from GAMLSS model
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @noRd
#' @export
#' @importFrom Matrix rankMatrix
#' 
get_norm_resid <- function(x, d){
  
  out <- lapply(x, function(dat){
    test <- dat$test
    train <- dat$train 
    modpred <- dat[1:d]
    
    fam_nam <- modpred[[1]]$model$family$family
    cdf <- modpred[[1]]$model$family$cdf
    if(is.null(cdf)){
      if(fam_nam == "gaulss"){
        cdf <- function(q, mu, wt, scale, logp = FALSE) {
          pnorm(q, mean = mu[ , 1], sd = 1/mu[ , 2], log.p = logp)
        }
      }
    }
 
    # For each marginal: get z by evaluating model CDF and inverting Gaussian CDF
    my_list <- lapply(modpred, function(o){
      
      y_nam <- as.character(o$model$formula[[1]][2])
      y_nam <- substr(y_nam,nchar(y_nam)-1,nchar(y_nam))
      ytrain <- o$model$y
      ytest <- test[ , paste0("res", y_nam)]
      ztrain <- qnorm(cdf(ytrain, o$model$fitted.values, logp = TRUE), log.p = TRUE)
      ztest <- qnorm(cdf(ytest, o$eta_hat, logp = TRUE), log.p = TRUE)
      
      return(list("ztrain" = ztrain, "ztest" = ztest))
    })
  
    return(my_list)
  })
  
  Z_train <- do.call("cbind", lapply(out[[1]], function(x) x$ztrain))
  Z_test <- do.call("rbind", lapply(out, function(x) do.call("cbind", lapply(x, function(x1) x1$ztest))))

  Z <- rbind(Z_train, Z_test)
  
  return(Z)
}


