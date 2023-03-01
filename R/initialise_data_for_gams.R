##########################
#' Initialise data for fitting
#' 
#' @description A
#'  
#' @param theta B
#' @return C
#' @details This function is meant for internal use only.
#' @export
#' 
initialise_data_for_gams <- function(dat, d){
  out <- data.frame(matrix(0,nrow(dat[["A"]]),ncol=d),
                    dow = dat[["A"]]$dow,
                    dow_RpH_EW = dat[["A"]]$dow_RpH,
                    dow_RpH_Sc= dat[["P"]]$dow_RpH,
                    clock_hour = dat[["A"]]$clock_hour,
                    Date = dat[["A"]]$doy,
                    doy = dat[["A"]]$doy,
                    year = dat[["A"]]$year,
                    moy = dat[["A"]]$moy,
                    doy_s = dat[["A"]]$doy_s,
                    doy_c = dat[["A"]]$doy_c,
                    doy_s2 = dat[["A"]]$doy_s2,
                    doy_c2 = dat[["A"]]$doy_c2,
                    t =  dat[["A"]]$t,
                    School_Hol_EW = dat[["A"]]$School_Hol,
                    School_Hol_S = dat[["P"]]$School_Hol,
                    n2ex = dat[["A"]]$n2ex,
                    EMBEDDED_WIND_CAPACITY = dat[["A"]]$EMBEDDED_WIND_CAPACITY)
  
  colnames(out)[1:d] <-c("res_A", "res_B", "res_C", "res_D", "res_E",
                         "res_F", "res_G", "res_H", "res_J", "res_K",
                         "res_L", "res_M", "res_N", "res_P")
  
  count_out_and_fix <- ncol(out)
  
  set_variables<-c("node_n_sm_L1","x2T_weighted.mean_p_max_point", "x2Tsm_point",
                   "SSRD_mean_2_Cap", "WindSpd100_weighted.mean_cell",
                   "WindSpd10_weighted.mean_cell", "TP_weighted.mean_cell")
  
  out <- cbind(out, matrix(0,nrow(dat[["A"]]),
                           ncol=d*length(set_variables)))
  for(j in 1:length(set_variables)){
    for(k in 1:d){
      count_out_and_fix <- count_out_and_fix + 1
      colnames(out)[count_out_and_fix] <- paste0(set_variables[j],"_",allNs[k])
      idx_var <- which(colnames(dat[[allNs[k]]])%in%set_variables[j] )
      out[,colnames(out)[count_out_and_fix]] <- dat[[allNs[k]]][[idx_var]]
    }
  }
  return( out )
}
