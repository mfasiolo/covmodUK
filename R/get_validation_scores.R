##########################
#' Get validation scores to choose number of effects
#' 
#' @param nstop_cross The output of the validation_step() function
#' @return A vector of out-of-sample log-likelihood values
#' @export
#' 
get_validation_scores <- function(nstop_cross){
  
  ngrid <- length(nstop_cross[[1]])
  
  # Rearrange predictions
  eta_hat <- lapply(1:ngrid,
                    function(ii){ do.call("rbind", lapply(nstop_cross, function(x) x[[ii]]$eta_hat)) } )
  
  # Sigmas <- lapply(eta_hat, function(x) apply(lpi_to_resp(x, 14)[ , -(1:14)], 1, Sigma_mat, simplify = FALSE))
  # COR <- lapply(Sigmas, function(x) lapply(x, cov2cor))
  # loglik <- sapply(1:ngrid, function(ii) -sum(sapply(1:nrow(y), function(kk) dmvn(y[kk, ], rep(0, 14), COR[[ii]][[kk]], log = TRUE))))
  
  y <- as.matrix( do.call("rbind", lapply(nstop_cross, function(x) x[[1]]$y)) )
  
  # Get log-likelihood at a function of number of effects
  loglik <- sapply(1:ngrid, function(ii) return( -nll_mcd(eta_hat[[ii]], y) ))
  
  return( loglik )
  
}