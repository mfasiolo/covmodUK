#' Base learner construction
#' @description aaa
#'
#' @export
## Base learner constructor: takes a smooth and returns a base learner
## for it, seeking to make effective degrees of freedom edf + null.space.dim.
## INPUTS
##  - X: full model matrix from a gamObject
##  - S: penalty matrix
##  - col: columns of X that are relevant to the specific smooth
##  - ns: dimension of null space
##  - edf: desired degrees of freedom of the base learner
## OUTPUT 
##  - Qf: a matrix such that mu_hat = Q %*% t(Q) %*% y
##
## EXAMPLE
#
# library(mgcv)
# dat <- gamSim(1, n = 400)
# fit <- gam(y ~ x1 + s(x2), data = dat, fit = FALSE)
# Q <- blcon(X = fit$X, fit$sm[[1]]$S[[1]], cols = fit$sm[[1]]$first.para:fit$sm[[1]]$last.para, ns = 1, edf = 4)
# plot(dat$x2, Q %*% (t(Q) %*% dat$y))
# Based on initial code from Simon Wood
blcon <- function(X, S, cols, ns, edf) {
  X <- X[ , cols]
  rf <- function(rho,X,S,edf) {
    n <- nrow(X)
    X1 <- rbind(X,exp(rho)*t(mroot(S)))
    QR <- qr(X1)
    R <- qr.R(QR)
    Q <- qr.Q(QR)[1:n, ]
    res <- sum(Q^2) - edf
    attr(res,"Q") <- Q
    attr(res,"R") <- R
    res
  } ## rf
  ## search for s.p. giving target EDF...
  if((ncol(X)-ns) <= edf){
    QR <- qr(X)
    return( list("Q" = qr.Q(QR), "R" = qr.R(QR)) )
  } 
  er <- uniroot(rf,c(-5,5),X=X,S=S,edf=edf,tol=.01)
  out <- list("rho" = er$root, "Q" = attr(er$f.root,"Q"), "R" = attr(er$f.root,"R"), "S" = S, "col" = col)
  return(out)
} ## blcon
