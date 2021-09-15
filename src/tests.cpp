/*
R package paropt
Copyright (C) 2021 Konrad Krämer

This file is part of R package paropt


paroptis free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with paropt
If not see: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html#SEC4
*/

#include <testthat.h>

#include "modify_dataframe.hpp"
#include "state.hpp"
#include "params.hpp"


context("Import parameter") {

  std::vector<int> params_cut_idx_vec;
  std::vector<double> params_time_combi_vec;
  std::vector<double> param_combi_lb;
  std::vector<double> param_combi_ub;
  std::vector<std::string> header_parameter;

  Rcpp::NumericVector time = {0, 2, 4, 6, 8, 10};
  Rcpp::NumericVector par1 = {1, 2, 3, 4, 5, 6};
  Rcpp::NumericVector par2 = {10, 20, 30, 40, 50, 60};

  Rcpp::DataFrame lower = Rcpp::DataFrame::create( Rcpp::Named("time") = time, Rcpp::Named("par1") = par1, Rcpp::Named("par2") = par2);

  Rcpp::NumericVector par3 = {10, 20, 30, 40, 50, 60};
  Rcpp::NumericVector par4 = {100, 200, 300, 400, 500, 600};
  Rcpp::DataFrame upper = Rcpp::DataFrame::create( Rcpp::Named("time") = time, Rcpp::Named("par1") = par3, Rcpp::Named("par2") = par4);

  ip (lower, upper, params_cut_idx_vec, params_time_combi_vec,
                   param_combi_lb, param_combi_ub, header_parameter);

  std::vector<int> cut_true{6, 6};

  test_that("Import parameter") {
    for(int i = 0; i < params_cut_idx_vec.size(); i++) {
        expect_true(params_cut_idx_vec[i] ==  cut_true[i]);
    }
  }
}
















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

Rcpp::XPtr<OS> test_solving() {
  Rcpp::XPtr<OS> xpfun = Rcpp::XPtr<OS>(new OS(&ode_system));

  return xpfun;
}
