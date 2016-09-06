#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <limits>
#include <list>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

double kl_divergence(mat a, mat b, double lambda) {
  return accu(a % log(a / (b + lambda) + lambda) - a + b);
}

// [[Rcpp::export]]
List c_znmf(mat v, int k, int niter, double tol, int peek_interval, double lambda) {
  int n = v.n_rows;
  int m = v.n_cols;
  
  double alpha = sqrt(2 * mean(mean(v)) / k);
  
  mat w = randu(n, k) * alpha;
  mat h = randu(k, m) * alpha;
  
  list<double> errors;
  double last_error = numeric_limits<double>::infinity();
  for (int i=0; i<niter; i++) {
    h = (h % (w.t() * (v / (w * h + lambda))));
    h.each_col() /= sum(w, 0).t();  
    w = (w % ((v / (w * h + lambda)) * h.t()));
    w.each_row() /= sum(h, 1).t();
    if ((i+1) % peek_interval == 0) {
      double error = kl_divergence(v, w * h, lambda);
      printf("(%d) error: %f; diff: %f\n", i+1, error, last_error - error);
      errors.push_back(error);
      if (last_error - error < tol) {
        break;
      }
      last_error = error;
    }
  }
  
  return List::create(
    Named("H") = h,
    Named("W") = w,
    Named("errors") = errors
  );
}
