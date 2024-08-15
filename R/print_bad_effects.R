#' Prints bad effects 
#' @return C
#' @noRd
#' @export
#'
print_bad_effects <- function(data_boost, n){
  
  bad_diff <- diff(data_boost$ll_v)
  bad_eff <- order(bad_diff) + 1
  
  for(ii in 1:n){
    
    aaa <- data_boost
    aaa$dll <- aaa$dll[1:(bad_eff[ii]-1)]
    aaa$ll_v <- aaa$ll_v[1:(bad_eff[ii]-1)]
    aaa$idx <- aaa$idx[1:(bad_eff[ii]-1)]
    aaa$nstep <- (bad_eff[ii]-1)
    aaa$eta <- NULL
    aaa$eta_v <- NULL
    aaa <- boost_eff(aaa)
    
    bbb <- data_boost
    bbb$dll <- bbb$dll[1:bad_eff[ii]]
    bbb$ll_v <- bbb$ll_v[1:bad_eff[ii]]
    bbb$idx <- bbb$idx[1:bad_eff[ii]]
    bbb$nstep <- bad_eff[ii]
    bbb$eta <- NULL
    bbb$eta_v <- NULL
    
    bbb <- boost_eff(bbb)
    
    print(paste(setdiff(paste0(bbb$name_eff, bbb$uni_label), paste0(aaa$name_eff, aaa$uni_label)), ", delta =",  bad_diff[bad_eff[ii]-1]))
  }
  
  return(NULL)
  
}