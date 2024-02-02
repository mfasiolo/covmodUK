#ifndef _SCM_H
#define _SCM_H



#include <RcppArmadillo.h>
using namespace Rcpp;

// mcd parametrisation
double d1_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w, arma::mat& G);
double d2_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v,
                  arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t,
                  Rcpp::List&  b1_eta, Rcpp::List&  b3, arma::vec& ib1_eta, arma::vec& ib3);

#endif
