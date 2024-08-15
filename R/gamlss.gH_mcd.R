#' New gamlss.gH_mcd
#'
#' @description  Routine for computing the derivatives
#' @param X model matrix
#' @param jj indices X for lpi
#' @param eta lpi
#' @param y outcome matrix
#' @param w indices vector
#' @param z indices vector
#' @param t indices vector
#' @param Gm indices matrix
#' @param l1 initialised score eta
#' @param l1_l initialised score eta (last block)
#' @param l2 initialised hessian eta
#' @param l2_v initialised hessian eta (cumulative)
#' @param l2_l initialised hessian eta (last block)
#' @param l2_v_l initialised hessian eta (cumulative - last block)
#' @param idx_b block indices
#' @param idx_aux list sparsity indices
#' @param l3 3rd der
#' @param i3 3rd der indices
#' @param l4 4th der
#' @param i4 4th der indices
#' @param d1b 1st der beta
#' @param d2b 2nd der beta
#' @param deriv integer for optimisation routine
#' @param fh Hessian inverse
#' @param D matrix ...
#'
#' @return derivatives
#' @noRd
#' @export
#'
#' @importFrom Rcpp evalCpp
#'
gamlss.gH_mcd <- function (X, jj, eta, y, w, z, t, Gm,
                           l1, l1_l, l2, l2_v, l2_l, l2_v_l, idx_b, idx_aux,
                           l3 = 0, i3 = 0,
                           l4 = 0, i4 = 0, d1b = 0, d2b = 0, deriv = 0, fh = NULL, D = NULL){
  K <- length(jj) # no_eta
  sparse <- discrete <- FALSE
  p <- ncol(X)
  n <- nrow(X)

  trHid2H <- d1H <- d2H <- NULL

  lpi <- lapply(jj, function(x) x - 1) # trick to allow the use of the c++functions
  #####################
  # Score w.r.t. beta #
  #####################
  lb <- rep(0, p)
  d1_beta(X, eta, y, lpi, K, lb, l1, l1_l, idx_b, z, w, Gm) #first derivatives w.r.t beta

  #######################
  # Hessian w.r.t. beta #
  #######################
  lbb <- if ( sparse )
    Matrix(0, p, p)
  else matrix(0, p, p) # dovremmo inizializzarla?????

  d2_beta(X, eta, y, lpi, K, lbb, l2, l2_v, l2_l, l2_v_l, idx_b, z, w, Gm, t, #second derivatives w.r.t beta
          idx_aux$b1_eta, idx_aux$b1, idx_aux$b2, idx_aux$b3, idx_aux$idx_b1_eta, idx_aux$idx_b3,
          idx_aux$l2_el, idx_aux$l2_el2)
  #######################################
  # Derivative of the Hessian w.r.t eta #
  #######################################

  if ( deriv > 0) {
    stop("error")
  }
  list(lb = lb, lbb = lbb, d1H = d1H, d2H = d2H, trHid2H = trHid2H)
}
