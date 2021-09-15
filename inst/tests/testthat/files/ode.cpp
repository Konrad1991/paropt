// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(paropt)]]
// [[Rcpp::plugins(cpp11)]]

typedef int (*OS)(double &t, std::vector<double> &params, std::vector<double> &states);

int ode_system(double &t, std::vector<double> &params, std::vector<double> & states) {

  // do not use any R-Code or R-Objects if the optimzation should run in parallel.
  // Users have to guarantee that the function can be called by several threads in parallel
  // define parameters (vector params contain the parameter in the order as defined in the corresponding textfiles)
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];

  // states have to be in the same order as specified in the textfile "states_LV.txt" (Otherwise the error-calculation does not work)
  double n1 = states[0];
  double n2 = states[1];

  states[0] = n1*c*n2 - n1*d;
  states[1] = n2*a - n2*b*n1;

  return 0;
}


// [[Rcpp::export]]
Rcpp::XPtr<OS> test_optimization() {
  Rcpp::XPtr<OS> xpfun = Rcpp::XPtr<OS>(new OS(&ode_system));

  return xpfun;
}
