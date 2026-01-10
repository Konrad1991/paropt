#include "header.hpp"

etr::Double default_error_fct(etr::Integer num_points, etr::Double a_inp, etr::Double b_inp) {
  const double a = get_val(a_inp);
  const double b = get_val(b_inp);
  return std::abs((a - b) / b) / static_cast<double>(get_val(num_points));
}

// [[Rcpp::export]]
Rcpp::XPtr<error_calc_fct> get_default_error_fct() {
  return Rcpp::XPtr<error_calc_fct>(new error_calc_fct(&default_error_fct));
}
