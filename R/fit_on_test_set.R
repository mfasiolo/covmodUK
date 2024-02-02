##########################
#' Fit MVN GAM on test set
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
fit_on_test_set <- function(data_boost, nstop, dat, start, ncv, ncores){
  
  mean_formula_int <- list(res_P ~ 1,res_N ~ 1,res_F ~ 1,res_M ~ 1,res_G ~ 1,res_D ~ 1,res_K ~ 1,
                           res_E ~ 1,res_B ~ 1,res_A ~ 1,res_C ~ 1,res_J ~ 1,res_H ~ 1,res_L ~ 1)
  d <- length(mean_formula_int)
  
  y_var_nam <- paste0("res", sapply(mean_formula_int, 
                                        function(form) {
                                          str <- as.character(form[2])
                                          substr(str,nchar(str)-1,nchar(str))
                                        }))

  if(nstop > 0){
  theta_formula <- formula_mcd(data_boost, stop_elem = nstop)
  theta_formula <- lapply(theta_formula,
                          function(x)
                            as.formula(paste(gsub(", sp = 0", '', x, fixed = TRUE), collapse = " ")))
  } else {
    theta_formula <- lapply(1:(d*(d+1)/2), function(nouse) ~ 1)
  }
  global_formula <- c(mean_formula_int,  theta_formula)
  
  ndat <- nrow(dat)
  sets <- floor( seq(start, ndat, length.out = ncv+1) )
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("dat", "sets", "global_formula", "y_var_nam"), 
                envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  my_fun_haffa <- function(ii){
    train <- dat[1:sets[ii], ]
    test <- dat[(sets[ii]+1):sets[ii+1], ]
    
    fit1 <- try(gam(global_formula,
                    family=mvn_mcd(d=14),
                    data = train, control=list(trace=TRUE)), TRUE)
    
    out <- list("eta_hat" = try(predict(fit1, newdata = test), TRUE), 
                "y" = test[ , y_var_nam])
    
    return(out)
  } 
  environment(my_fun_haffa) <- .GlobalEnv
  
  test_eval <-  parLapply(NULL, 1:(length(sets)-1), my_fun_haffa)
  
  stopCluster(cl)
  rm(cl)
  
  return( test_eval )
}