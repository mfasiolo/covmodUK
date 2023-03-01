##########################
#' Boosting ordering
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
boosting_step <- function(data, data_v, ord_idx, effects, 
                          constraints, nstep, rate, diagonal, ncores, ncv = NULL, eff_seq = NULL, verbose = FALSE){
  
  if( is.null(eff_seq) ){
  out <- boost_order(y = as.matrix(data[ , ord_idx]), 
                     X =  data[ , -c(1:14)],
                     y_v = as.matrix(data_v[ , ord_idx]), 
                     X_v =  data_v[ , -c(1:14)],
                     effects = effects,
                     constraints = constraints, nstep = nstep, rate = rate,
                     diagonal = diagonal,
                     ncores = ncores, eff_seq = eff_seq, verbose = verbose)
  } else {
    
    ndat <- nrow(data_v)
    sets <- floor(seq(1, ndat, length.out = ncv+1) )
 
    out <- mclapply(1:ncv, function(ii){
                     y <- as.matrix(data[ , ord_idx])
                     X <- data[ , -c(1:14)]
                     if(ii > 1){
                       y <- rbind(y, as.matrix(data_v[sets[1]:sets[ii], ord_idx]))
                       X <- rbind(X, data_v[sets[1]:sets[ii], -c(1:14)])
                     }
                      res <- boost_order(
                       y = y, 
                       X = X,
                       y_v = as.matrix(data_v[(sets[ii]+1):sets[ii+1], ord_idx]), 
                       X_v =  data_v[(sets[ii]+1):sets[ii+1], -c(1:14)],
                       effects = effects,
                       constraints = constraints, nstep = nstep, rate = rate,
                       diagonal = diagonal,
                       ncores = 1, eff_seq = eff_seq, verbose = FALSE)$ll_v
                      return(unlist(res))
      }, mc.cores = ncores)
    
  }
  
  return( out )
  
}
