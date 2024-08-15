#' Post process marginal GPD fits
#'
#' @param fit_s a list of fitted marginal shash models
#' @param fit_g a list of fitted marginl gam models
#' @param only_shash if TRUE the GPD distribution will not be used
#' @return A list containing the fit obtained by a shash+gpd model
#' @export
#'
postprocess_gpd_fits <- function(fit_s, fit_g, only_shash = FALSE){
  
  d <- length(fit_s[[1]])
  al <- attr(fit_g, "alpha_lev")
  mu_s <- lapply(1:d, function(ii) rbind(fit_s[[1]][[ii]]$mu_train, do.call("rbind", lapply(fit_s, function(o) o[[ii]]$mu_test))))
  llk_s <- sapply(1:d, function(ii) c(fit_s[[1]][[ii]]$llk_train, do.call("c", lapply(fit_s, function(o) o[[ii]]$llk_test))))
  r_s <- sapply(1:d, function(ii) c(fit_s[[1]][[ii]]$resid_train, do.call("c", lapply(fit_s, function(o) o[[ii]]$resid_test))))
  
  mu_gpd <- lapply(1:d, 
                   function(ii) {
                     list("up" = rbind(fit_g[[1]][[ii]][["up"]]$mu_train_all, do.call("rbind", lapply(fit_g, function(o) o[[ii]][["up"]]$mu_test_all))), 
                          "low" = rbind(fit_g[[1]][[ii]][["low"]]$mu_train_all, do.call("rbind", lapply(fit_g, function(o) o[[ii]][["low"]]$mu_test_all))))
                   })
  
  if( !only_shash ){
    
    gpd_dat_in <- lapply(fit_g[[1]], function(o){
      index <- c(o[["up"]]$index_train, o[["low"]]$index_train)
      r <- log(c(1 - al + al*exp(o[["up"]]$resid_train), al - al*exp(o[["low"]]$resid_train)))
      llk <- c(o[["up"]]$llk_train, o[["low"]]$llk_train) + log(al)
      q_s <- c(o[["up"]]$q_shash_train, o[["low"]]$q_shash_train)
      mu <- rbind(o[["up"]]$mu_train, o[["low"]]$mu_train)
      dir <- c(rep("up", length(o[["up"]]$index_train)), rep("low", length(o[["low"]]$index_train)))
      out <- data.frame("r" = r, "llk" = llk, "q_s" = q_s, "index" = index, "dir" = dir)
      out$mu <- mu
      if(!length(index)) stop("No residuals!")
      return( out )
    })
    
    gpd_dat_out <- lapply(1:d, function(ii) do.call("rbind", lapply(fit_g, function(o){
      index <- c(o[[ii]][["up"]]$index_test, o[[ii]][["low"]]$index_test)
      r <- log(c(1 - al + al*exp(o[[ii]][["up"]]$resid_test), al - al*exp(o[[ii]][["low"]]$resid_test)))
      llk <- c(o[[ii]][["up"]]$llk_test, o[[ii]][["low"]]$llk_test) + log(al)
      q_s <- c(o[[ii]][["up"]]$q_shash_test, o[[ii]][["low"]]$q_shash_test)
      mu <- rbind(o[[ii]][["up"]]$mu_test, o[[ii]][["low"]]$mu_test)
      dir <- c(rep("up", length(o[[ii]][["up"]]$index_test)), rep("low", length(o[[ii]][["low"]]$index_test)))
      out <- data.frame("r" = r, "llk" = llk, "q_s" = q_s, "index" = index, "dir" = dir)
      out$mu <- mu
      # if(!length(index)) stop("No residuals!")
      return( out )
    })))
    
    
    nas_in <- lapply(1:d, function(ii) which(!is.finite(gpd_dat_in[[ii]]$r) | gpd_dat_in[[ii]]$r == 0))
    if(sum(sapply(nas_in, length)) > 0) {
      stop(paste0("There are ", sum(sapply(nas_in, length)), " NAs in the in-sample GPD residuals"))
    }
    
    gpd_dat <- lapply(1:d, function(ii) rbind(gpd_dat_in[[ii]], gpd_dat_out[[ii]]))
    
    nas <- lapply(1:d, function(ii) which(!is.finite(gpd_dat[[ii]]$r) | gpd_dat[[ii]]$r == 0))
    print(paste0("There are ", sum(sapply(nas, length)), " NAs in the out-of-sample GPD residuals"))
    # for(ii in 1:d){
    #   for(kk in nas[ii]){
    #     gpd_dat[[ii]]$llk[kk] <- min(gpd_dat[[ii]]$llk, na.rm = TRUE)
    #     gpd_dat[[ii]]$r[kk] <- ifelse(gpd_dat[[ii]]$dir[kk] == "up", 
    #                                       max(gpd_dat[[ii]]$r, na.rm = TRUE), min(gpd_dat[[ii]]$r, na.rm = TRUE)) 
    #   }
    # }
  }
  
  resid_fun <- function(index){
    out <- r_s[index, ]
    if(!only_shash){
      for(ii in 1:d){
        info <- gpd_dat[[ii]]
        kk <- which(info$index == index)
        if(length(kk) && !(kk %in% nas[[ii]])){
          out[ii] <- info$r[kk]
        }
      }
    }
    return(out)
  }
  
  qf_fun <- function(index, p, logp){
    lp <- p
    if(!logp){
      lp <- log(p)
    }
    lal <- log(al)
    q <- shash_qf(lp, mu = t(sapply(mu_s, function(x) x[index, ])), logp = TRUE)
    if(!only_shash){
      ind <- which(lp < lal | lp > log1p(-al))
      for(kk in ind){
        info <- gpd_dat[[kk]]
        ee <- which(info$index == index)
        if( length(ee) && !(ee %in% nas[[kk]]) ){
          dir <- ifelse(lp[kk] < lal, "low", "up")
          lp0 <- lp[kk]
          lp0 <- ifelse(lp0 < lal, log1p(-exp(lp0-lal)), log1p(-(1-exp(lp0))/al))
          lp_sh <- ifelse(lp[kk] < lal, lal, log1p(-al))
          qs <- shash_qf(p = lp_sh, mu = mu_s[[kk]][index, ], logp = TRUE) 
          q[kk] <- qs + gpd_qf(lp0, mu_gpd[[kk]][[dir]][index, ], logp = TRUE) * ifelse(dir == "up", 1, -1)
        }
      }
    }
    return(q)
  }
  
  rd_fun <- function(index){
    p <- runif(d)
    y <- qf_fun(index = index, p = p, logp = FALSE)
    return(y)
  }
  
  llk_fun <- function(index){
    out <- llk_s[index, ]
    if(!only_shash){
      for(ii in 1:d){
        info <- gpd_dat[[ii]]
        kk <- which(info$index == index)
        if( length(kk) && !(kk %in% nas[[ii]]) ){
          out[ii] <- info$llk[kk]
        }
      }
    }
    return(out)
  } 
  
  return( list("rd" = rd_fun, "llk" = llk_fun, "resid" = resid_fun, "qf" = qf_fun) )
  
}
