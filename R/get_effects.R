##########################
#' Creates list of effects for gradient boosting 
#' 
#' @param version what type of model are we using "c", "c+w" etc
#' @param allNs vector of names of GSP groups
#' @param diagonal if TRUE the effects will be used to model only diagonal of MCD
#' @return List of effects to be passed to boosting_step() function
#' @export
#' 
get_effects <- function(version, allNs, diagonal){
  
  d <- 14
  
  G <- mat2vec(d)
  
  effects <- c("dow", "t + I(t^2)", "s(clock_hour, sp = 0)", "s(doy, sp = 0)")
  constraints <- list()
  
  if( grepl("w", version, fixed = TRUE) ){
    wind_eff <- sapply(1:d,
                       function(ii) paste0("s(WindSpd100_weighted.mean_cell_",
                                           allNs[ii], ", by = EMBEDDED_WIND_CAPACITY, k = 5, sp = 0)" ))
    solar_eff <- sapply(1:d,
                        function(ii) paste0("s(SSRD_mean_2_Cap_",
                                            allNs[ii], ", k = 5, sp = 0)" ))
    effects <- c(effects, wind_eff, solar_eff)
    tmp <- rep(lapply(1:d, function(ii) G[ii, 1:ii]), 2)
    names(tmp) <- c(wind_eff, solar_eff)
    constraints <- c(constraints, tmp)
  }
  
  if( grepl("r", version, fixed = TRUE) ){
    temp_eff <- sapply(1:d,
                       function(ii) paste0("s(x2T_weighted.mean_p_max_point_",
                                           allNs[ii], ",  k = 5, sp = 0)" ))
    prec_eff <- sapply(1:d,
                       function(ii) paste0("s(TP_weighted.mean_cell_",
                                           allNs[ii], ",  k = 5, sp = 0)" ))
    effects <- c(effects, temp_eff, prec_eff, "s(n2ex, k = 5, sp = 0)") 
    tmp <- rep(lapply(1:d, function(ii) G[ii, 1:ii]), 2)
    names(tmp) <- c(temp_eff, prec_eff)
    constraints <- c(constraints, tmp)
  }
  
  if(diagonal){
    for(ii in 1:length(effects))
      if(effects[ii] %in% names(constraints)){
        constraints[[effects[ii]]] <- constraints[[effects[ii]]][constraints[[effects[ii]]] <= d]
      } else {
        constraints[[effects[ii]]] <- 1:d
      }
  }
  
  return( list("effects" = effects, "constraints" = constraints) )
}
  
  