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
fit_on_test_set <- function(data_boost, nstop, data, start, ncv, ncores){
  
  mean_formula_int <- list(res_P ~ 1,res_N ~ 1,res_F ~ 1,res_M ~ 1,res_G ~ 1,res_D ~ 1,res_K ~ 1,
                           res_E ~ 1,res_B ~ 1,res_A ~ 1,res_C ~ 1,res_J ~ 1,res_H ~ 1,res_L ~ 1)
  d <- length(mean_formula_int)
  
  y_var_nam <- paste0("res_out", sapply(mean_formula_int, 
                                        function(form) {
                                          str <- as.character(form[2])
                                          substr(str,nchar(str)-1,nchar(str))
                                        }))
  
  res <- boost_eff(data_boost)
  
  if(nstop > 0){
  theta_formula <- formula_mcd(data_boost, res, stop_elem = nstop)
  theta_formula <- lapply(theta_formula,
                          function(x)
                            as.formula(paste(gsub(", sp = 0", '', x, fixed = TRUE), collapse = " ")))
  } else {
    theta_formula <- lapply(1:(d*(d+1)/2), function(nouse) ~ 1)
  }
  global_formula <- c(mean_formula_int,  theta_formula)
  
  ndat <- nrow(data)
  sets <- floor( seq(start, ndat, length.out = ncv+1) )
  
  test_eval <- mclapply(1:(length(sets)-1), function(ii){
    train <- residuals_data[1:sets[ii], ]
    test <- residuals_data[(sets[ii]+1):sets[ii+1], ]
    
    fit1 <- try(gam(global_formula,
                    family=mvn_mcd(d=14),
                    data = train, control=list(trace=TRUE)), TRUE)
    
    out <- list("eta_hat" = try(predict(fit1, newdata = test), TRUE), 
                "y" = test[ , y_var_nam])
    
    return(out)
  }, mc.cores = ncores)
  
  return( test_eval )
}