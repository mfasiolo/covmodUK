#' Family definition (internal use)
#'
#' @description  Family function
#' @param d dimension of the outcome
#' @param nb number of chunkings
#'
#' @return family
#' @export
#'
#' @importFrom Rcpp evalCpp
#'
#' @examples

mvn_mcd <- function(d = 2, nb = 1){
  #Multivariate model with covariance modelling via modified Cholesky decomposition
  if ( d < 2 ) stop("mvn_scm requires to or more dimensional data")
  if ( nb < 1 ) stop("the number of observations' blocks must be greater than 0")

  no_eta <- d + d *(d + 1)/2
  param <- 1 #set to MCD parametrisation


  stats <- list()
  for (j in 1:no_eta) stats[[j]] <- make.link("identity")

  validmu <- function(mu) all(is.finite(mu))


  # z, w, G, t: auxiliary indices used in the derivative's computation (all of them used in the mcd and only z, w in the logm)
  z <- w <- t <- rep(0, (d * (d - 1)/2))
  Gm <- matrix(0, d - 1, d - 1)
  mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
  idx_zwGt(d, z, w, Gm, t)

  # indices of elements different from zero in the hessian matrix
  idx_jk <- idxHess_no0(no_eta, z, w, param)

  if ( param == 1 ) { # See the computational paper
    nHel <- d * (d^2 + 15 * d + 2)/6
    if ( d > 2 ) nHel <- nHel + d * (d - 1) * (d - 2)/3
  }

  # Quantities defined in the environment:
  # l1 and l1_l: matrix of the 1st derivatives
  getL1 <- function() get(".l1")   # only the first nb - 1 blocks of observations
  putL1 <- function(.l1) assign(".l1", .l1, envir = environment(sys.function()))
  getL1_l <- function() get(".l1_l")   # Last block of observations
  putL1_l <- function(.l1_l) assign(".l1_l", .l1_l, envir = environment(sys.function()))

  # l2 and l2_l: matrices of the 2nd derivatives excluding the intercepts' blocks
  getL2 <- function() get(".l2") # only the first nb - 1 blocks of observations
  putL2 <- function(.l2) assign(".l2", .l2, envir = environment(sys.function()))
  getL2_l <- function() get(".l2_l")   # Last block of observations
  putL2_l <- function(.l2_l) assign(".l2_l", .l2_l, envir = environment(sys.function()))

  # l2_v and l2_v_l: vectors of the (cumulated) second derivatives considering only the  intercepts' blocks
  getL2_v <- function() get(".l2_v")  # First nb - 1 blocks of observations
  putL2_v <- function(.l2_v) assign(".l2_v", .l2_v, envir = environment(sys.function()))
  getL2_v_l <- function() get(".l2_v_l")   # Last block of observations
  putL2_v_l <- function(.l2_v_l) assign(".l2_v_l", .l2_v_l, envir = environment(sys.function()))

  # Indices for observations' blocks
  getidx_b <- function() get(".idx_b")
  putidx_b <- function(.idx_b) assign(".idx_b", .idx_b, envir = environment(sys.function()))

  # List of auxiliary indices for hessian blocks  (intercept/partial/full taking into account the sparsity if it exists)
  getidx_aux <- function() get(".idx_aux")
  putidx_aux <- function(.idx_aux) assign(".idx_aux", .idx_aux, envir = environment(sys.function()))


  assign(".cflag", TRUE, envir = environment())
  getcflag <- function() get(".cflag")
  putcflag <- function(.cflag) assign(".cflag", .cflag, envir = environment(sys.function()))


  assign(".d", d, envir = environment())
  getd <- function() get(".d")
  putd <- function(.d) assign(".d", .d, envir = environment(sys.function()))

  assign(".no_eta", no_eta, envir = environment())
  getno_eta <- function() get(".no_eta")
  putno_eta <- function(.no_eta) assign(".no_eta", .no_eta, envir = environment(sys.function()))

  # Type of parametrisation
  assign(".param", param, envir = environment())
  getparam <- function() get(".param")
  putparam <- function(.param) assign(".param", .param, envir = environment(sys.function()))

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      d <- family$getd() # d <- ncol(y)
      no_eta <- family$getno_eta()
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      jj <- attr(x,"lpi")
      resid <-  matrix(0, nrow(y), d)
      if(is.list(x)){ # discrete case
        start <- rep(0,max(unlist(jj)))
        for(k in 1:d){
          yt1 <- y[,k]
          e1 <- E[ , jj[[k]], drop=FALSE] ## square root of total penalty
          R <- suppressWarnings(chol(XWXd(x$Xd, w = rep(1, length(yt1)),
                                          k = x$kd, ks = x$ks, ts = x$ts,
                                          dt = x$dt, v = x$v, qc = x$qc, nthreads = 1,
                                          drop = x$drop, lt = x$lpid[[k]]) +
                                       crossprod(E[ , jj[[k]]]), pivot = TRUE))
          Xty <- XWyd(x$Xd, rep(1, nrow(y)), yt1, x$kd,
                      x$ks, x$ts, x$dt, x$v, x$qc, x$drop, lt = x$lpid[[k]])
          piv <- attr(R, "pivot")
          rrank <- attr(R, "rank")
          startji <- rep(0, ncol(R))
          if (rrank < ncol(R)) {
            R <- R[1:rrank, 1:rrank]
            piv <- piv[1:rrank]
          }
          startMu <- numeric(ncol(R))
          startMu[piv] <- backsolve(R, forwardsolve(t(R), Xty[piv]))
          startMu[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startMu
          eta1 <- Xbd(x$Xd, start, k = x$kd, ks = x$ks,
                      ts = x$ts, dt = x$dt, v = x$v, qc = x$qc, drop = x$drop,
                      lt = x$lpid[[k]])
          resid[,k] <-  y[,k] -   eta1
        }
        Svcov <- cov(resid)
        MCD_dec <- mcd(Svcov)
        Theta_el<-c(diag(MCD_dec), MCD_dec[upper.tri(MCD_dec, diag=FALSE)])

      } else { #regular case
        start <- rep(0,ncol(x))
        for(k in 1:d){
          yt1 <- y[,k]
          x1 <- x[ , jj[[k]], drop=FALSE]
          e1 <- E[ , jj[[k]], drop=FALSE] ## square root of total penalty

          if(use.unscaled){
            x1 <- rbind(x1, e1)
            startMu <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
            startMu[ !is.finite(startMu) ] <- 0
          } else {
            startMu <- pen.reg(x1, e1, yt1)
          }
          if(!is.matrix(x[ , jj[[k]]])){
            resid[,k] <-  y[,k] - x[ , jj[[k]]]*startMu
          } else{
            resid[,k] <-  y[,k] - x[ , jj[[k]]]%*%startMu
          }
          start[jj[[k]]] <- startMu
        }

        Svcov <- cov(resid)
        MCD_dec <- mcd(Svcov)
        Theta_el<-c(diag(MCD_dec), MCD_dec[upper.tri(MCD_dec, diag=FALSE)])

        for(k in (d+1):no_eta){
          x1 <-  x[ , jj[[k]],drop=FALSE]
          startji <- qr.coef(qr(x1), c(rep(Theta_el[k-d],nrow(x1))))
          startji[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startji
        }
      }
      return(start)
    }
    if(is.null(start)){
      start <- my_init_fun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset)
    }
  }) ## initialize

  ##To do residuals and postproc
  residuals <- function(object, type=c("response","deviance")) {
    type <- match.arg(type)
    if(type == "deviance"){
      n <- dim(object$fitted.values)[1]
      res <- matrix(0, n, d)
      res_dev_mcd(object$fitted.values, object$y, res)
    } else {
      res <- object$y - object$fitted.values[,1:d]
    }
    res
  } ## residuals


  ll <- function(y, X, coef, wt, family,  offset = NULL,
                 deriv = 0, d1b = NULL, d2b = NULL, Hp = NULL,
                 rank = 0, fh = NULL, D = NULL) {
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index

    p <- ncol(X)
    n <- nrow(y)
    eta <- matrix(0, n, no_eta) #linear predictor matrix

    # number of observations in each block (nb - 1 blocks with the same observations; the last could differ)
    nlast <- n %% nb
    nset <- n %/% nb

    # Memory initialization of several quantities:
    l1 <- try(getL1(), TRUE) #First derivatives matrix (first nb - 1 blocks of observations)
    if ( "try-error" %in% class(l1) ) {
      l1 <- matrix(0, nset, no_eta)
      putL1(l1)
    }

    idx_aux <- try(getidx_aux(), TRUE) #list of auxiliary indices for hessian blocks (intercept/partial/full taking into account the sparsity if it exists)
    if ( "try-error" %in% class(idx_aux) ) {
      idx_aux <- aux_idx(jj, idx_jk, no_eta)
      putidx_aux(idx_aux)
    }

    if(nb > 1){
      l2 <- try(getL2(), TRUE)  #Second derivatives matrix (first nb - 1 blocks of observations and Hessian blocks not involving the intercepts)
      if ( "try-error" %in% class(l2) ) {
        l2 <- matrix(0, nset, nHel - idx_aux$llls)
        putL2(l2)
      }
    } else { # Such trick allows to pass the matrix l2 in the case nb = 1, avoiding cpp issues
      l2 <- try(getL2(), TRUE)  #Second derivatives matrix (first nb - 1 blocks of observations and Hessian blocks not involving the intercepts)
      if ( "try-error" %in% class(l2) ) {
        l2 <- matrix(0, 1, 1)
        putL2(l2)
      }
    }

    l2_v <- try(getL2_v(), TRUE) #Second derivatives vector (first nb - 1 blocks of observations and Hessian blocks involving the intercepts)
    if ( "try-error" %in% class(l2_v) ) {
      l2_v <- rep(0, idx_aux$llls)
      putL2_v(l2_v)
    }

    l2_v_l <- try(getL2_v_l(), TRUE)  #Second derivatives vector (last blocks of observations and Hessian blocks not involving the intercepts)
    if ( "try-error" %in% class(l2_v_l) ) {
      l2_v_l <- rep(0, idx_aux$llls)
      putL2_v_l(l2_v_l)
    }

    # Create the last block quantities
    if ( nlast == 0 ) {
      nobs_b <- nset
      idx_b_seq <- cumsum(rep(nset, nb - 1)) - 1
    } else {
      nobs_b <- nlast
      idx_b_seq <- cumsum(rep(nset, nb)) - 1
    }

    l1_l <- try(getL1_l(), TRUE)  #First derivatives matrix (Last block of observations)
    if ( "try-error" %in% class(l1_l) ) {
      l1_l <- matrix(0, nobs_b,  no_eta)
      putL1_l(l1_l)
    }

    l2_l <- try(getL2_l(), TRUE) #Second derivatives matrix (Last block of observations)
    if("try-error" %in% class(l2_l)){
      l2_l <- matrix(0, nobs_b, nHel - idx_aux$llls)
      putL2_l(l2_l)
    }

    idx_b <- try(getidx_b(), TRUE)  # Indices of observations' blocks
    if("try-error" %in% class(idx_b)){
      idx_b <- c(-1, idx_b_seq, n - 1)
      putidx_b(idx_b)
    }



    for(k in 1 : no_eta){
      eta[,k] <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[k]])
      else X[,jj[[k]],drop=FALSE]%*%coef[jj[[k]]]
    }

    #l1 <- matrix(0, n, no_eta) #initialization
    ## log-likelihood: eta is a matrix n*w and y is a matrix n*d

    #l <- ll_mcd(eta, y[,1:d]) - 0.5 * n * d * log(2 * pi)
    l <- ll_mcd(eta, y) - 0.5 * n * d * log(2 * pi)


    #if (deriv>0) {
      ## the first derivative: eta is a matrix n*w, y is a matrix n*d,
      ##                         l1 is a matrix n*w
      #d1_mcd(eta, y[,1:d], l1)
      #d1_mcd(eta, y, l1)

      ## the second derivatives
      #l2 <- matrix(0, n, no_eta * (no_eta + 1)/2) #initialization
      #d2_mcd(eta,y[,1:d], l2)
      #d2_mcd(eta,y, l2)

    #}

    l3 <- 0 ## defaults


    # if (deriv) {
    #   i2 <- family$tri$i2
    #   i3 <- family$tri$i3
    #
    #   ## get the gradient and Hessian...
    #   ret <- gamlss.gH(X, jj, l1, l2, i2, l3 = l3, i3 = i3,
    #                    d1b = d1b, deriv = deriv - 1, fh = fh, D = D)
    # } else ret <- list()

    if ( deriv ) {
      ret <- gamlss.gH_mcd(X, jj, eta, y, w, z, t, Gm,
                                      l1, l1_l, l2, l2_v, l2_l, l2_v_l,
                                      idx_b, idx_aux,
                                      d1b = d1b, deriv = deriv - 1, fh = fh, D = D)
    } else ret <- list()

    ret$l <- l
    ret
  } ## end ll mvn_mcd

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL) {
    ## optional function to give predicted values - idea is that
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information.
    ## if se = FALSE returns one item list containing matrix otherwise
    ## list of two matrices "fit" and "se.fit"...


    if (is.null(eta)) {
      discrete <- is.list(X)
      lpi <- attr(X,"lpi")
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      }
      K <- length(lpi) ## number of linear predictors
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,K)
      for (i in 1:K) {
        if (discrete) {
          eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]])
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        }
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
      }
    }

    out <- matrix(0, nrow(eta), ncol(eta))
    cor_flag <- getcflag()
    pred_mcd(eta, out,d, as.integer(cor_flag))
    list(fit = out)
  } ## mvncm predict

  jacobian <- function(eta, jj){
    #The following two lines could be unuseful
    #no_eta <- ncol(eta) # delete in SCM
    #d <- -3/2 + sqrt(9/4 + 2 * no_eta) # delete in SCM
    res <- matrix(0, nrow(eta), no_eta)
    G <- mat2vec(d) #put in cpp???

    cor_flag <- getcflag()

    if(jj <= d){
      res[,jj] <- 1
    } else {
      jj <- jj - d
      idx_jj <- which(G == jj, arr.ind = TRUE)
      S_row <- as.numeric(idx_jj[1, 1])
      S_col <- as.numeric(idx_jj[1, 2])

      rc_idx_s <- rep(NA, d * (d - 1)/2)
      rc_idx_t <- rep(NA, d * (d - 1)/2)

      count <- 1
      for(j in (d + 1):(d * (d + 1)/2)){
        rc_idx_s[count] <- which(G == j, arr.ind=TRUE)[1, 1]
        rc_idx_t[count] <- which(G == j, arr.ind=TRUE)[1, 2]
        count <- count + 1
      }

      jacobian_mcd(eta, res, d, S_row - 1, S_col - 1,
                   as.numeric(rc_idx_s) - 1, as.numeric(rc_idx_t) - 1, as.integer(cor_flag))
    }
    return(res)
  }

  structure(list(family = "Multivariate normal (MCD)", ll = ll, nlp = no_eta,
                 ##link=paste(link), ## unuseful?
                 tri = trind.generator(no_eta,ifunc=TRUE),
                 initialize = initialize,
                 getd = getd, putd = putd,
                 getno_eta = getno_eta, putno_eta = putno_eta,
                 getcflag = getcflag, put_cflag = putcflag,
                 getL1 = getL1, putL1 = putL1,
                 getL1_l = getL1_l, putL1_l = putL1_l,
                 getL2 = getL2, putL2 = putL2,
                 getL2_v = getL2_v, putL2_v = putL2_v,
                 getL2_l = getL2_l, putL2_l = putL2_l,
                 getL2_v_l = getL2_v_l, putL2_v_l = putL2_v_l,
                 getidx_b = getidx_b, putidx_b = putidx_b,
                 getidx_aux =  getidx_aux, putidx_aux =  putidx_aux,
                 getparam =  getparam, putparam =  putparam,
                 #postproc=postproc, ##to do
                 residuals=residuals,
                 predict = predict,
                 jacobian = jacobian,
                 linfo = stats, ## link information list
                 validmu = validmu,
                 #rd=rd,
                 d2link = 1, d3link = 1, d4link = 1, ## signals to fix.family.link that all done
                 ls = 1, ## signals that ls not needed here
                 available.derivs = 0, ## can use full Newton here
                 discrete.ok = TRUE
  ),class = c("general.family", "extended.family", "family"))
} ## end mvn_cmd
