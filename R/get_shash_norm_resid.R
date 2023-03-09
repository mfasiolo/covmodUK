##########################
#' Get normalized residuals from shash model
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @export
#' @importFrom Matrix rankMatrix
#' @details This function is meant for internal use only.
#' 
get_shash_norm_resid <- function(x, d, ncores){
  
  out <- lapply(x, function(dat){
    test <- dat$test
    train <- dat$train 
    modpred <- dat[1:14]
    
    fam_nam <- modpred[[1]]$model$family$family
      my_pdf <- shash_pdf
      cdf <- modpred[[1]]$model$family$cdf
 
    # For each marginal: get z by evaluating shash CDF and inverting Gaussian CDF
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
  
  #out <- rbind(out[[1]]$ztrain, do.call("rbind", lapply(out, "[[", "ztest"))) 
  
  Z <- rbind(Z_train, Z_test)
  
  return(Z)
}


