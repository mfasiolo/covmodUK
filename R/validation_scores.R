##########################
#' Get validation scores
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
get_validation_scores <- function(nstop_cross){
  
  ngrid <- length(nstop_cross[[1]])
  
  # Rearrange predictions
  eta_hat <- lapply(1:ngrid,
                    function(ii){ do.call("rbind", lapply(nstop_cross, function(x) x[[ii]]$eta_hat)) } )
  
  y <- as.matrix( do.call("rbind", lapply(nstop_cross, function(x) x[[1]]$y)) )
  
  # Get log-likelihood at a function of number of effects
  loglik <- sapply(1:ngrid, function(ii) return( -nll_mcd(eta_hat[[ii]], y) ))
  
  return( loglik )
  
}