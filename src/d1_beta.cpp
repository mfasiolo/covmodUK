#include "scm_der.h"

//' Score mcd for fitting
//'
//' @param X model matrix
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param jj list of indices to access X elements (by lpi)
//' @param K length ...
//' @param lb Score beta
//' @param l1 Memory initialization (first B-1 blocks)
//' @param l1_l Memory initialization (Last block)
//' @param ig indices for chunking
//' @param z indices vector
//' @param w indices vector
//' @param G indices matrix
//' @export

// [[Rcpp::export(name="d1_beta")]]
double d1_beta(arma::mat& X, arma::mat& eta,  arma::mat& y, Rcpp::List& jj, uint32_t& K, arma::vec& lb, arma::mat& l1, arma::mat& l1_l,  arma::uvec& ig, arma::vec& z, arma::vec& w, arma::mat& Gm){
  using namespace arma;
  uint32_t size_g = ig.n_elem - 1;
  uint32_t igf;
  uint32_t igl;

  uint32_t i;
  uint32_t j;

  Rcpp::IntegerVector ijjj;
  uint32_t l_jjj;
  uint32_t ijjjf;
  uint32_t ijjjl;


  for(i = 0; i < size_g; i++){
    igf = ig[i] + 1;
    igl = ig[i + 1];

    if ( i < (size_g - 1) ) {
        d1_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l1, z, w, Gm);
    } else {
        d1_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l1_l, z, w, Gm);
    }


      for(j = 0; j < K; j++) {
        ijjj = jj[j];
        l_jjj = ijjj.length() - 1;
        ijjjf = ijjj[0];
        ijjjl = ijjj[l_jjj];
        if ( i < (size_g - 1) ) {
          lb.subvec(ijjjf, ijjjl) += X(span(igf, igl), span(ijjjf, ijjjl)).t() * l1.col(j);
        } else {
          lb.subvec(ijjjf, ijjjl) += X(span(igf, igl), span(ijjjf, ijjjl)).t() * l1_l.col(j);
        }
      }
  }
  return(1.0);
}
