#ifndef HEADERS
#define HEADERS

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

#include <cassert>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_math.h>
#include <cvode/cvode_diag.h> // for ADAMS
#include <arkode/arkode_erkstep.h> // for ERK
#include <arkode/arkode_arkstep.h> // for fully implicit systems
#include "modify_dataframe.hpp"
#include "param_interpolation.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <stdio.h>

#include <stdlib.h>
#include <vector>
#include <list>
#include <iterator>
#include <random>

#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cctype>

#include <thread>
#include <mutex>

#if _OPENMP
#include <omp.h>
#endif

#define NA std::nan("l")

struct settingsPSO {
  double err_tol;
  int pso_n_pop;
  int pso_n_gen;
  double pso_par_initial_w;
  double pso_par_w_max;
  double pso_par_w_min;
  double pso_par_w_damp;
};

struct settingsPSO_Rcpp_interface {
  double err_tol;
  int pso_n_pop;
  int pso_n_gen;
  double pso_par_initial_w;
  double pso_par_w_max;
  double pso_par_w_min;
  double pso_par_w_damp;
};

struct time_state_information {
  std::vector<double> init_state;
  std::vector<double> par_times;
  std::vector<int> param_idx_cuts;
  std::vector<double> state_measured;
  std::vector<double> state_times;
  std::vector<int> state_idx_cut;
  Rcpp::NumericVector integration_times;
  double reltol;
  Rcpp::NumericVector absolute_tolerances;
};

struct time_state_information_two_stage {
  std::vector<double> par_times;
  std::vector<int> param_idx_cuts;
  std::vector<double> state_measured;
  std::vector<double> states_derivative_measured;
  std::vector<double> state_times;
  std::vector<int> state_idx_cut;
  Rcpp::NumericVector integration_times;
};

struct time_state_information_Rcpp_interface {
  std::vector<double> init_state;
  std::vector<double> par_times;
  std::vector<int> param_idx_cuts;
  std::vector<double> state_measured;
  std::vector<double> state_times;
  std::vector<int> state_idx_cut;
  std::vector<double> integration_times;
  double reltol;
  std::vector<double> absolute_tolerances;
};

/*
print std vector
*/
template<typename T>
void pv(T const &s) {
  for(auto Data : s) {
    Rcpp::Rcout << Data << std::endl;
  }
}


#endif // HEADERS
