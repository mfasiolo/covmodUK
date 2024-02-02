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
get_final_fits <- function(nstop_all, data, ncores, files){
  
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("nstop_all", "data", "ncores", "files"), 
                envir = environment())
  clusterEvalQ(NULL, {
    library(covmodUK)
  })
  
  my_fun_ahfuosab <- function(ii){
    boost_nam <- files[ii]
    
    load(file = boost_nam)
    
    mean_formula_int <- list(res_P ~ 1,res_N ~ 1,res_F ~ 1,res_M ~ 1,res_G ~ 1,res_D ~ 1,res_K ~ 1,
                             res_E ~ 1,res_B ~ 1,res_A ~ 1,res_C ~ 1,res_J ~ 1,res_H ~ 1,res_L ~ 1)
    d <- length(mean_formula_int)
    
    if(nstop_all[ii] > 0){
      theta_formula <- formula_mcd(data_boost, stop_elem = nstop_all[ii])
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
    
  }
  environment(my_fun_ahfuosab) <- .GlobalEnv
  
  out <- parLapply(NULL, 1:length(nstop_all), my_fun_ahfuosab)
  
  stopCluster(cl)
  rm(cl)
  
  return(out)
  
}