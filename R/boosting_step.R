##########################
#' Order effects by gradient boosting
#' 
#' @param data training data
#' @param data_v validation data
#' @param effect list of effects
#' @param constraints list of constraints
#' @param nstep number of boosting steps
#' @param rate learning rate
#' @param diagonal if TRUE only diagonal of MCD is modelled
#' @param ncores number of parallel thread to use
#' @param ncv number of fold for cross-validation
#' @param old_fit can be use to start boosting from an old boosting run
#' @param verbose print information as we boost?
#' @return A list containing information on the effects and their importance
#' @export
#' 
boosting_step <- function(data, data_v, effects, 
                          constraints, nstep, rate, diagonal, ncores, ncv, old_fit = NULL, verbose = FALSE){
  
  out <- boost_order(y = as.matrix(data[ , 1:14]), 
                     X =  data,
                     y_v = as.matrix(data_v[ , 1:14]), 
                     X_v =  data_v,
                     effects = effects,
                     constraints = constraints, nstep = nstep, rate = rate,
                     diagonal = diagonal,
                     old_fit = old_fit,
                     ncores = ncores, verbose = verbose)
  
  ndat <- nrow(data_v)
  sets <- floor(seq(1, ndat, length.out = ncv+1) )
  
  eff_seq <- list("index" = sapply(out$idx, function(x) x[1,1]))
  eff_seq$name <- out$effects[sapply(out$idx, function(x) x[1,2])]
  
  cv_stuff <- NULL
  if( !is.null(old_fit) ){
    cv_stuff <- list("eta" = old_fit$cv_eta, "eta_v" = old_fit$cv_eta_v)
  }
  
  ncores <- min(c(ncv, ncores))
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("data", "data_v", "effects", "sets", "eff_seq", "constraints", "nstep", 
                        "rate", "diagonal", "cv_stuff"), envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  cv_runs <- parLapply(NULL, 1:ncv, function(ii){
    y <- as.matrix(data[ , 1:14])
    X <- data
    if(ii > 1){
      y <- rbind(y, as.matrix(data_v[sets[1]:sets[ii], 1:14]))
      X <- rbind(X, data_v[sets[1]:sets[ii], ])
    }
    my_old_fit <- NULL
    if(!is.null(cv_stuff)){
      my_old_fit <- list("eta" = cv_stuff$eta[[ii]], "eta_v" = cv_stuff$eta_v[[ii]])
    }
    res <- boost_order(
      y = y, 
      X = X,
      y_v = as.matrix(data_v[(sets[ii]+1):sets[ii+1], 1:14]), 
      X_v =  data_v[(sets[ii]+1):sets[ii+1], ],
      effects = effects,
      constraints = constraints, nstep = nstep, rate = rate,
      diagonal = diagonal,
      ncores = 1, eff_seq = eff_seq, old_fit = my_old_fit, verbose = FALSE)
    return(res)
  })
  
  stopCluster(cl)
  rm(cl)

  out$ll_v <- Reduce("+", lapply(cv_runs, function(x) unlist(x$ll_v)))
  
  out$cv_eta <- lapply(cv_runs, "[[", "eta")
  out$cv_eta_v <- lapply(cv_runs, "[[", "eta_v")
  
  if( !is.null(old_fit) ){
   out$ll_v <- c(old_fit$ll_v, out$ll_v)
   out$idx <- c(old_fit$idx, out$idx)
   out$dll <- c(old_fit$dll, out$dll)
   out$nstep <- old_fit$nstep + nstep
  }

  return( out )
  
}
