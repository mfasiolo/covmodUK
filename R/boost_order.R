#' Perform the gradient boosting
#' @description This function allows to obtain an ordering of importance of the effects on the unconstrained elements of the modified cholesky decomposition
#'
#' @param y outcome (residual)
#' @param X covariates
#' @param effects list of effects included in the boosting procedure
#' @param constraints list of constraints on the effects
#' @param nstep numbers of step of the boosting procedure
#' @param verbose to print the iteration summary of the procedure
#' @param rate learning rate
#' @param ncores number of cores for parallel computation
#'
#' @return A list of summaries
#' @export
#' @importFrom parallel mclapply
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom stats as.formula
#' @importFrom stats cov
#'
#' @examples
#'
boost_order <- function(y, X,
                        effects, constraints,
                        y_v, X_v,
                        nstep, 
                        rate, diagonal, ncores,
                        old_fit = NULL,
                        eff_seq = NULL, verbose = FALSE){
  d <- ncol(y)
  n <- nrow(y)
  n_v <-  nrow(y_v)
  neta <- d * (d + 1)/2

  if(diagonal){
    elem_seq <- 1:d
  } else {
    elem_seq <- 1:(d*(d + 1)/2)
  }
  
  if( is.null(old_fit) ){ # Initialise linear predictors OR....
    eta <- matrix(0, n, neta + d)
    eta_v <- matrix(0, n_v, neta + d)
    Svcov <- cov(y)
    MCD_dec <- mcd(Svcov)
    Theta_el <- c(diag(MCD_dec), MCD_dec[upper.tri(MCD_dec, diag=FALSE)])
    for(j in elem_seq){
      eta[ , j+d]  <- rep(Theta_el[j], n)
      eta_v[ , j+d] <- rep(Theta_el[j], n_v)
    }
  } else { # ... restart from latest values
    eta <- old_fit$eta
    eta_v <- old_fit$eta_v
  }
  
  score_init <- matrix(0, n, neta + d)
  
  l0 <- -nll_mcd(eta, y)
  
  baselearn <- lapply(effects, function(eff){
    # Build model matrix for effect and its QR decomposition
    form_gam <- as.formula(paste0("res_P ~ ", eff)) # Any response if fine here
    gam_obj <- gam(form_gam, data=X, fit = FALSE)
    Xmat <- gam_obj$X
    sm <- gam_obj$smooth
    if(!length(sm)){
      S <- diag(1, ncol(Xmat))
      cols <- 1:ncol(Xmat)
      ns <- 0
    } else {
      S <- gam_obj$smooth[[1]]$S[[1]]
      cols <- sm[[1]]$first.para:sm[[1]]$last.para
      ns <- sm[[1]]$null.space.dim
    }
    out <- blcon(X = Xmat, S = S, cols = cols, ns = ns, edf=4)
    out$Xval <- predict(gam(G = gam_obj), newdata = X_v, type = "lpmatrix")[ , cols]
    return( out )
    
  })
  names(baselearn) <- effects
  
  my_fit_base_learn <- function(eff){
    delta <- rep(0, neta)
    index <- constraints[[eff]]
    if( is.null(index) ) index <- 1:neta
    Q <- baselearn[[eff]]$Q
    # Fit effect to all relevant linear predictors
    for(jj in index){
      # Least squares fit to score vector
      score_hat <- Q %*% (t(Q) %*% score_init[ , jj + d])
      # Improve the following steps: avoid creating a copy
      eta1 <- eta
      eta1[ , jj+d] <- eta1[ , jj+d] + rate * score_hat
      delta[jj] <- -nll_mcd(eta1, y) - l0
    }
    return(delta)
  }
  environment(my_fit_base_learn) <- .GlobalEnv

  dll <- list()
  idx <- list()
  ll_v <- list()
  
  ncores <- min(c(ncores, length(effects)))
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("neta", "y", "rate", "baselearn", "d"), envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  for(ii in 1:nstep){
    d1_mcd(eta, y, score_init)
    
    if( !is.null(eff_seq) ){
      effects <- eff_seq$name[ii]
      constraints <- list(eff_seq$index[ii])
      names(constraints) <- effects
    }
    
    clusterExport(NULL, c("effects", "constraints", "score_init", "eta", "l0"), envir = environment())
    
    dll[[ii]] <- parLapply(NULL, effects, my_fit_base_learn)
    
    dll[[ii]] <- do.call("cbind", dll[[ii]])
    
    l0 <- l0 + max(dll[[ii]])
    
    idx[[ii]] <- which(dll[[ii]] == max(dll[[ii]]), arr.ind = TRUE)
    
    yyy <- score_init[, idx[[ii]][1, 1] + d]
    
    best_bl <- baselearn[[ effects[idx[[ii]][1, 2]] ]]
    beta_hat <- drop( backsolve(best_bl$R, t(best_bl$Q) %*% yyy) )
  
    yyy_hat <- drop( best_bl$Q %*% (t(best_bl$Q) %*% yyy) )
    
    eta[, d + idx[[ii]][1,1]] <- eta[, d + idx[[ii]][1,1]] + rate * yyy_hat
    
    eta_v[, d + idx[[ii]][1,1]] <- eta_v[, d + idx[[ii]][1,1]] + rate * drop(best_bl$Xval %*% beta_hat)
    
    ll_v[[ii]] <- -nll_mcd(eta_v, y_v)
    
    if(verbose == TRUE){
      print(paste("iter" = ii, "delta" = max(dll[[ii]]),
                  "id" = length(unique(sapply(idx, function(.x) paste0(.x[1,1], ".", .x[1,2])))),
                  "ll_v" = ll_v[[ii]]))
    }
    if(max(dll[[ii]]) < 0) break
  }
  
  stopCluster(cl)
  rm(cl)
  
  return(list(dll = dll, idx = idx, d = d, nstep = nstep, effects = effects, ll_v = ll_v, eta = eta, eta_v = eta_v))
}
