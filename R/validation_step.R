##########################
#' Validation for optimal number of effects
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
validation_step <- function(res, data, start, nstop, ncv, ncores){
  
  mean_formula_int <- list(res_P ~ 1,res_N ~ 1,res_F ~ 1,res_M ~ 1,res_G ~ 1,res_D ~ 1,res_K ~ 1,
                           res_E ~ 1,res_B ~ 1,res_A ~ 1,res_C ~ 1,res_J ~ 1,res_H ~ 1,res_L ~ 1)
  
  y_var_nam <- paste0("res", sapply(mean_formula_int, 
                                        function(form) {
                                          str <- as.character(form[2])
                                          substr(str,nchar(str)-1,nchar(str))
                                        }))
  
  d <- length(mean_formula_int)
  
  Theta_formula_int <- list()
  for(ii in 1:(d*(d + 1)/2)){
    Theta_formula_int[[ii]] <- as.formula("~ 1")
  }
  
  ndat <- nrow(data)
  sets <- floor(seq( start, ndat, length.out = ncv+1) )
  
  nstop <- pmin(nstop, nrow(res))
  
  nstop_cross <- mclapply(1:(length(sets)-1), function(ii){
    train <- residuals_data[1:sets[ii], ]
    test <- residuals_data[(sets[ii]+1):sets[ii+1], ]
    
    test_pred <- lapply(nstop, function(ns){
      
      theta_formula <- formula_mcd(data_boost, res, stop_elem = ns)
      theta_formula <- lapply(theta_formula,
                              function(x)
                                as.formula(paste(gsub(", sp = 0", '', x, fixed = TRUE), collapse = " ")))
      global_formula <- c(mean_formula_int,  theta_formula)
      
      fit1 <- try(gam(global_formula,
                      family=mvn_mcd(d=14),
                      data = train, control=list(trace=FALSE)), TRUE)
      
      out <- list("eta_hat" = try(predict(fit1, newdata = test), TRUE), 
                  "y" = test[ , y_var_nam])
      
      return(out)
      
    })
    
    return( test_pred )
  }, mc.cores = ncores)
  
  
  return(nstop_cross)
}
