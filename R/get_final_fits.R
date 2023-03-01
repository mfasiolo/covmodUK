##########################
#' Fit all MVN GAM on all data
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
get_final_fits <- function(nstop_all, data, ncores){
  
  out <- mclapply(1:4, function(ii){
    boost_nam <- switch(as.character(ii), 
                        "1" = "data_boost_optimal_0.1_c.RData", 
                        "2" = "data_boost_optimal_0.1_c.RData", 
                        "3" = "data_boost_optimal_0.1_c+w.RData",
                        "4" = "data_boost_optimal_0.1_c+w+r.RData")
    
    load(file = boost_nam)
    
    mean_formula_int <- list(res_P ~ 1,res_N ~ 1,res_F ~ 1,res_M ~ 1,res_G ~ 1,res_D ~ 1,res_K ~ 1,
                             res_E ~ 1,res_B ~ 1,res_A ~ 1,res_C ~ 1,res_J ~ 1,res_H ~ 1,res_L ~ 1)
    d <- length(mean_formula_int)
    
    res <- boost_eff(data_boost)
    
    if(nstop_all[ii] > 0){
      theta_formula <- formula_mcd(data_boost, res, stop_elem = nstop_all[ii])
      theta_formula <- lapply(theta_formula,
                              function(x)
                                as.formula(paste(gsub(", sp = 0", '', x, fixed = TRUE), collapse = " ")))
    } else {
      theta_formula <- lapply(1:(d*(d+1)/2), function(nouse) ~ 1)
    }
    global_formula <- c(mean_formula_int,  theta_formula)
    
    fit1 <- try(gam(global_formula,
                    family=mvn_mcd(d=d),
                    data = data, control=list(trace=TRUE)), TRUE)
    
    return(fit1)
    
  }, mc.cores = ncores)
  
  return(out)
  
}