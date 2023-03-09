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
boosting_step <- function(data, data_v, effects, 
                          constraints, nstep, rate, diagonal, ncores, ncv = NULL, verbose = FALSE){
  
  out <- boost_order(y = as.matrix(data[ , 1:14]), 
                     X =  data,
                     y_v = as.matrix(data_v[ , 1:14]), 
                     X_v =  data_v,
                     effects = effects,
                     constraints = constraints, nstep = nstep, rate = rate,
                     diagonal = diagonal,
                     ncores = ncores, verbose = verbose)
  
  ndat <- nrow(data_v)
  sets <- floor(seq(1, ndat, length.out = ncv+1) )
  
  eff_seq <- list("index" = sapply(out$idx, function(x) x[1,1]))
  eff_seq$name <- out$effects[sapply(out$idx, function(x) x[1,2])]
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("data", "data_v", "effects", "sets", "eff_seq", "constraints", "nstep", 
                        "rate", "diagonal"), envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  ll_v <- parLapply(NULL, 1:ncv, function(ii){
    y <- as.matrix(data[ , 1:14])
    X <- data
    if(ii > 1){
      y <- rbind(y, as.matrix(data_v[sets[1]:sets[ii], 1:14]))
      X <- rbind(X, data_v[sets[1]:sets[ii], ])
    }
    res <- boost_order(
      y = y, 
      X = X,
      y_v = as.matrix(data_v[(sets[ii]+1):sets[ii+1], 1:14]), 
      X_v =  data_v[(sets[ii]+1):sets[ii+1], ],
      effects = effects,
      constraints = constraints, nstep = nstep, rate = rate,
      diagonal = diagonal,
      ncores = 1, eff_seq = eff_seq, verbose = FALSE)$ll_v
    return(unlist(res))
  })
  
  stopCluster(NULL)
  rm(cl)
  
  out$ll_v <- Reduce("+", ll_v)
  
  return( out )
  
}
