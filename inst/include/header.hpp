#ifndef HEAD
#define HEAD

// [[Rcpp::depends(ast2ast)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppThread.h>
#include "etr.hpp"

#include <cassert>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_math.h>
#include <cvode/cvode_diag.h> // for ADAMS

#include <limits>
#include <vector>

struct time_state_information {
  std::vector<double> init_state;
  std::vector<double> par_times;
  std::vector<int> param_idx_cuts;
  std::vector<double> state_measured;
  std::vector<int> state_idx_cut;
  std::vector<double> integration_times;
  double reltol;
  std::vector<double> absolute_tolerances;
};

typedef void (*OS)(double t, sexp& y, sexp& ydot, sexp& params);
typedef std::vector<double> vd;
typedef std::vector<int> vi;
typedef arma::vec av;
typedef arma::mat am;
typedef double (*solver_ptr)(std::vector<double> &param_combi_start, OS ode_system, time_state_information &solv_param_struc);
typedef double (*solver_ptr_save)(std::vector<double> &param_combi_start, OS ode_system, time_state_information solv_param_struc, Rcpp::NumericMatrix &DF);

#endif // HEAD
