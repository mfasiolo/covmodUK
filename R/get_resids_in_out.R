##########################
#' Get uniform residuals from GAMLSS models on test and training sets
#'  
#' @param mod The output of the fitted GAMLSS models
#' @param final Same as above but for the models fitted to the whole data
#' @param logp if TRUE the uniform residuals are on log scale
#' @return A list of in-sample and out-of-sample residuals
#' @export
#' 
get_resids_in_out <- function(mod, final, logp){
  # In sample residuals (train set)
  d <- length(final)
  cdf <- get(paste0(final[[1]]$family$family, "_cdf"))
  resid_in <- sapply(final,
                     function(o){
                       return( cdf(o$y, o$fitted.values, logp = logp) )
                     })
  
  # Out sample residuals (test set)
  resid_out <- lapply(mod, function(o){
    sapply(o, "[[", "resid_test")
  })
  resid_out <- do.call("rbind", resid_out)
  if(!logp){ 
    resid_out <- exp(resid_out)
  }
  return( list("in" = resid_in, "out" = resid_out) )
}


