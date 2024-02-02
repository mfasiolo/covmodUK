#' Fit a GAMLSS model to each margin
#'
#' @description XXX
#' @param X XXX
#'
#' @return XXX
#' @export
#'
#' @examples
fit_gamlss <- function(fam, form, sets, dat, cl){
  
  cdf <- get(paste0(fam, "_cdf"))
  pdf <- get(paste0(fam, "_pdf"))
  qf <- get(paste0(fam, "_qf"))
  rd <- get(paste0(fam, "_rd"))
  clusterExport(NULL, c("form", "sets", "dat", "fam", "cdf", "pdf", "rd", "qf"), envir = environment())
  
  my_fun_faff <- function(ii){
    train <- dat[1:sets[ii], ]
    test <- dat[(sets[ii]+1):sets[ii+1], ]
    out <- lapply(1:length(form), function(kk){
      my_form <- form[[kk]]
      fit1 <- gam(my_form, family = get(fam)(), data = train, optimizer = "efs")
      mu_test <- predict(fit1, newdata = test, type = "response")
      resid_train <- cdf(train[ , kk], fit1$fitted.values, logp = TRUE)
      resid_test <- cdf(test[ , kk], mu_test, logp = TRUE)
      llk_train <- pdf(train[ , kk], fit1$fitted.values, log = TRUE)
      llk_test <- pdf(test[ , kk], mu_test, log = TRUE)
      return( list("mu_train" = fit1$fitted.values, "mu_test" = mu_test, 
                   "resid_train" = resid_train, "resid_test" = resid_test, 
                   "llk_train" = llk_train, "llk_test" = llk_test, 
                   "y_train" = train[ , kk], "y_test" = test[ , kk]) )
    })
    out$test <- test
    out$train <- train
    return(out)
  }
  environment(my_fun_faff) <- .GlobalEnv
  
  fits <- parLapply(NULL, 1:(length(sets)-1), my_fun_faff)
  attr(fits, "nam") <- fam
  attr(fits, "train") <- lapply(fits, "[[", "train")
  attr(fits, "test") <- lapply(fits, "[[", "test")
  
  for(ii in 1:length(fits)){
    fits[[ii]]$train <- fits[[ii]]$test <- NULL
  }
  
  return(fits)
}
