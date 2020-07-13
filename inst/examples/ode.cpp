#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector ode(double t, std::vector<double> params, Rcpp::NumericVector y) {
  Rcpp::NumericVector ydot(y.length());

  // define parameters
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];

  //Rcpp::Rcerr << a << "\t" << b << "\t" << c << "\t" << d << std::endl;
  //Rcpp::Rcerr << y[0] << "\t" << y[1] << std::endl;
  // states
  double n1 = y[0];
  double n2 = y[1];

  ydot[0] = n1*c*n2 - n1*d;
  ydot[1] = n2*a - n2*b*n1;
  return ydot;
}
