#ifndef HEAD
#define HEAD

// [[Rcpp::depends(ast2ast)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp2a)]]
#include "etr.hpp"
#include <RcppThread.h>

#include <cassert>
#include <cvode/cvode.h>
#include <cvode/cvode_diag.h> // for ADAMS
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include <limits>
#include <vector>

typedef etr::Vec<double> (*error_calc_fct)(int num_points, double a,
                                           double b);
typedef etr::Vec<double> (*spline_fct)(double &t, etr::Vec<double> &time_vec,
                                       etr::Vec<double> &par_vec);
typedef void (*JAC)(double &t, etr::Vec<double, etr::Borrow<double>> &y,
                    etr::Vec<double, etr::Borrow<double>> &ydot,
                    etr::Mat<double, etr::Borrow<double>> &J,
                    etr::Vec<double, etr::Borrow<double>> &params);

struct time_state_information {
  std::vector<double> init_state;
  std::vector<double> par_times;
  std::vector<int> param_idx_cuts;
  std::vector<double> state_measured;
  std::vector<int> state_idx_cut;
  std::vector<double> integration_times;
  double reltol;
  std::vector<double> absolute_tolerances;
  error_calc_fct ecf;
  spline_fct sf;
  JAC jf;
};

typedef void (*OS)(double &t, etr::Vec<double, etr::Borrow<double>> &y,
                   etr::Vec<double, etr::Borrow<double>> &ydot,
                   etr::Vec<double, etr::Borrow<double>> &params);

typedef std::vector<double> vd;
typedef std::vector<int> vi;
typedef arma::vec av;
typedef arma::mat am;
typedef double (*solver_ptr)(std::vector<double> &param_combi_start,
                             OS ode_system,
                             time_state_information &solv_param_struc);
typedef double (*solver_ptr_save)(std::vector<double> &param_combi_start,
                                  OS ode_system,
                                  time_state_information solv_param_struc,
                                  Rcpp::NumericMatrix &DF);

#endif // HEAD
