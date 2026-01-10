/*
 BSD 3-Clause License
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution. The names of its
 contributors may not be used to endorse or promote products derived from this
 software without specific prior written permission. THIS SOFTWARE IS PROVIDED
 BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED
 WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "header.hpp"

typedef std::vector<double> VD;
typedef std::vector<int> VI;
typedef std::vector<std::vector<double>> MD;
typedef std::vector<std::vector<int>> MI;
typedef std::vector<std::string> VS;
typedef Rcpp::DataFrame DF;

using etr::Array;
using etr::Buffer;
using etr::get_val;
using etr::Double;

void params_sort(realtype &t, std::vector<double> &params,
                 std::vector<double> &par_vec, std::vector<double> &time_vec,
                 std::vector<int> &param_idx_cuts) {
  // Get number of parameter
  const std::size_t no_par = param_idx_cuts.size();

  params.resize(no_par);

  Array<Double, Buffer<Double>> tmp_time_vec;
  Array<Double, Buffer<Double>> tmp_par_vec;
  int tmp_no_vals;
  std::size_t idx_count = 0;

  for (std::size_t i = 0; i < no_par; ++i) {
    tmp_no_vals = static_cast<int>(param_idx_cuts[i]);
    if (tmp_no_vals == 1) {
      params[i] = par_vec[idx_count];
      ++idx_count;
    } else {
      tmp_par_vec = etr::numeric(etr::Integer(tmp_no_vals));
      tmp_time_vec = etr::numeric(etr::Integer(tmp_no_vals));
      for (int j = 0; j < tmp_no_vals; ++j) {
        tmp_par_vec.set(j, etr::Double(par_vec[idx_count]));
        tmp_time_vec.set(j, etr::Double(time_vec[idx_count]));
        ++idx_count;
      }
      params[i] = get_val(etr::cmr(etr::Double(t), tmp_time_vec, tmp_par_vec));
    }
  }
}

static int check_retval(void *returnvalue, const char *funcname, int opt) {
  int *retval;
  if (opt == 0 && returnvalue == NULL) {
    return (1);
  } else if (opt == 1) {
    retval = (int *)returnvalue;
    if (*retval < 0) {
      return (1);
    }
  } else if (opt == 2 && returnvalue == NULL) {
    return (1);
  }

  return (0);
}

// Rcpp::Rcerr is not thread safe!!!
void own_error_handler(int error_code, const char *module, const char *function,
                       char *msg, void *usr_data) {}

struct usr_data {
  OS ode_system;
  std::vector<double> parameter;
  std::vector<double> parameter_time;
  std::vector<int> parameter_cut_idx;
};

int wrapper_ode_system(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  // cast pointer to structure and store elements
  struct usr_data *my_ode_system = (struct usr_data *)user_data;
  OS odes;
  odes = (*my_ode_system).ode_system;
  std::vector<double> params = (*my_ode_system).parameter;
  std::vector<double> params_time = (*my_ode_system).parameter_time;
  std::vector<int> params_cut_idx = (*my_ode_system).parameter_cut_idx;

  // interpolate Parameter if necessary
  std::vector<double> parameter_input;
  params_sort(t, parameter_input, params, params_time, params_cut_idx);

  // extract time
  etr::Double time = t;

  etr::Array<etr::Double, etr::Borrow<etr::Double>> parameter(
    parameter_input.data(),
    parameter_input.size(),
    std::vector<std::size_t>{parameter_input.size()}
  );
  double* y_ptr = N_VGetArrayPointer(y);
  double* yd_ptr = N_VGetArrayPointer(ydot);
  const std::size_t n_y = static_cast<std::size_t>(NV_LENGTH_S(y));
  const std::size_t n_ydot = static_cast<std::size_t>(NV_LENGTH_S(ydot));

  if (!y_ptr || !yd_ptr) return -1;
  etr::Array<etr::Double, etr::Borrow<etr::Double>> y_(y_ptr, n_y, std::vector<std::size_t>{n_y});
  etr::Array<etr::Double, etr::Borrow<etr::Double>> ydot_(yd_ptr, n_ydot, std::vector<std::size_t>{n_ydot});
  odes(time, y_, ydot_, parameter);

  return 0;
}

std::mutex mtx2;

double solver_bdf(std::vector<double> &param_combi_start, OS ode_system,
                  time_state_information &solv_param_struc) {
  mtx2.lock();
  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec =
      solv_param_struc.state_measured;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;
  error_calc_fct ecf = solv_param_struc.ecf;
  mtx2.unlock();

  // Begin Solver
  SUNContext sunctx;
  int NEQ = hs_cut_idx_vec.size();
  realtype reltol, t;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;

  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  // Create the SUNDIALS context
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1))
    return (1);

  // initial conditions
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0))
    return (1);
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0))
    return (1);

  mtx2.lock();
  for (int i = 0; i < NEQ; ++i) {
    NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
    NV_Ith_S(y, i) = solv_param_struc.init_state[i];
  }

  reltol = solv_param_struc.reltol;
  mtx2.unlock();

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0))
    return (1);

  // set error handler to Rcpp::Rcerr
  // void* ptr_to_cvode_mem = &cvode_mem;
  void *ptr_to_nothing = &cvode_mem;
  retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler, ptr_to_nothing);

  double sum_of_least_squares = 0.;

  struct usr_data my_ode_system = {ode_system,
                                   param_combi_start,
                                   params_time_combi_vec,
                                   params_cut_idx_vec};
  void *ptr_to_my_ode_system = &my_ode_system;
  retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
  if (check_retval((void *)cvode_mem, "CVodeSetUserData", 0))
    return (1);

  retval = CVodeInit(cvode_mem, wrapper_ode_system, integration_times[0], y);
  if (check_retval(&retval, "CVodeInit", 1))
    return (1);

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1))
    return (1);

  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0))
    return (1);

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0))
    return (1);

  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1))
    return (1);

  double return_time;
  std::vector<double> temp_measured(init_state.size());

  for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
    return_time = integration_times[ti];
    retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
    for (int n = 0; n < NV_LENGTH_S(y); n++) {
      temp_measured[n] = hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
      if (std::isnan(temp_measured[n])) {
      } else {
        sum_of_least_squares += get_val(ecf(
          etr::Integer(static_cast<int>(integration_times.size())),
          etr::Double(NV_Ith_S(y, n)),
          etr::Double(temp_measured[n])));
      }
    }
    if (retval < 0) {
      sum_of_least_squares = std::numeric_limits<double>::max();
      break;
    }
  }

  N_VDestroy(y);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&sunctx);

  return sum_of_least_squares;
}

double solver_bdf_save(std::vector<double> &param_combi_start, OS ode_system,
                       time_state_information solv_param_struc,
                       Rcpp::NumericMatrix &DF) {

  mtx2.lock();
  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec =
      solv_param_struc.state_measured;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;
  error_calc_fct ecf = solv_param_struc.ecf;
  mtx2.unlock();

  // Begin Solver
  SUNContext sunctx;
  int NEQ = hs_cut_idx_vec.size();
  realtype reltol, t;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;

  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  // Create the SUNDIALS context
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1))
    return (1);

  // initial conditions
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0))
    return (1);
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0))
    return (1);

  mtx2.lock();
  for (int i = 0; i < NEQ; ++i) {
    NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
    NV_Ith_S(y, i) = solv_param_struc.init_state[i];
  }

  reltol = solv_param_struc.reltol;
  mtx2.unlock();

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0))
    return (1);

  // set error handler to Rcpp::Rcerr
  // void* ptr_to_cvode_mem = &cvode_mem;
  void *ptr_to_nothing = &cvode_mem;
  retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler, ptr_to_nothing);

  double sum_of_least_squares = 0.;

  struct usr_data my_ode_system = {ode_system,
                                   param_combi_start,
                                   params_time_combi_vec,
                                   params_cut_idx_vec};
  void *ptr_to_my_ode_system = &my_ode_system;
  retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
  if (check_retval((void *)cvode_mem, "CVodeSetUserData", 0))
    return (1);

  retval = CVodeInit(cvode_mem, wrapper_ode_system, integration_times[0], y);
  if (check_retval(&retval, "CVodeInit", 1))
    return (1);

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1))
    return (1);

  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0))
    return (1);

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0))
    return (1);

  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1))
    return (1);

  double return_time;
  std::vector<double> temp_measured(init_state.size());

  for (int i = 0; i < NV_LENGTH_S(y); i++) {
    DF(0, i) = NV_Ith_S(y, i);
  }

  int counter = 1;
  for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
    return_time = integration_times[ti];
    retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL); //
    for (int n = 0; n < NV_LENGTH_S(y); n++) {
      DF(counter, n) = NV_Ith_S(y, n);
      temp_measured[n] = hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
      if (std::isnan(temp_measured[n])) {
      } else {
        sum_of_least_squares += get_val(ecf(
          etr::Integer(static_cast<int>(integration_times.size())),
          etr::Double(NV_Ith_S(y, n)),
          etr::Double(temp_measured[n])));
      }
    }
    if (check_retval(&retval, "CVode", 1)) {
      break;
    }
    counter++;
  }

  N_VDestroy(y);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&sunctx);

  return sum_of_least_squares;
}

double solver_adams(std::vector<double> &param_combi_start, OS ode_system,
                    time_state_information &solv_param_struc) {

  mtx2.lock();
  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec =
      solv_param_struc.state_measured;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;
  error_calc_fct ecf = solv_param_struc.ecf;
  mtx2.unlock();

  // Begin Solver
  SUNContext sunctx;
  int NEQ = hs_cut_idx_vec.size();
  realtype reltol, t; // tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int retval; // retvalr, iout;

  // Create the SUNDIALS context
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1))
    return (1);

  y = abstol = NULL;
  cvode_mem = NULL;

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0))
    return (1);
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0))
    return (1);

  mtx2.lock();
  for (int i = 0; i < NEQ; ++i) {
    NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
    NV_Ith_S(y, i) = solv_param_struc.init_state[i];
  }

  reltol = solv_param_struc.reltol;
  mtx2.unlock();
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0))
    return (1);

  // set error handler to Rcpp::Rcerr
  // void* ptr_to_cvode_mem = &cvode_mem;
  void *ptr_to_nothing = &cvode_mem;
  retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler, ptr_to_nothing);

  double sum_of_least_squares = 0.;

  struct usr_data my_ode_system = {ode_system, param_combi_start,
                                   params_time_combi_vec, params_cut_idx_vec};
  void *ptr_to_my_ode_system = &my_ode_system;
  retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
  if (check_retval((void *)cvode_mem, "CVodeSetUserData", 0))
    return (1);

  retval = CVodeInit(cvode_mem, wrapper_ode_system, integration_times[0], y);
  if (check_retval(&retval, "CVodeInit", 1))
    return (1);

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1))
    return (1);

  retval = CVDiag(cvode_mem);
  if (check_retval(&retval, "CVDiag", 1))
    return (1);
  double return_time;

  std::vector<double> temp_measured(init_state.size());

  for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
    return_time = integration_times[ti];
    retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
    for (int n = 0; n < NV_LENGTH_S(y); n++) {
      temp_measured[n] = hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
      if (std::isnan(temp_measured[n])) {
      } else {
        sum_of_least_squares += get_val(ecf(
          etr::Integer(static_cast<int>(integration_times.size())),
          etr::Double(NV_Ith_S(y, n)),
          etr::Double(temp_measured[n])));
      }
    }
    if (retval < 0) {
      sum_of_least_squares =
          std::numeric_limits<double>::max(); // 1.79769e+308;//1000000.;
      break;
    }
  }

  N_VDestroy(y);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return sum_of_least_squares;
}

double solver_adams_save(std::vector<double> &param_combi_start, OS ode_system,
                         time_state_information solv_param_struc,
                         Rcpp::NumericMatrix &DF) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec =
      solv_param_struc.state_measured;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;
  error_calc_fct ecf = solv_param_struc.ecf;

  // Begin Solver
  SUNContext sunctx;
  int NEQ = hs_cut_idx_vec.size();
  realtype reltol, t; // tout;
  N_Vector y, abstol;
  void *cvode_mem;
  int retval; // retvalr, iout;

  // Create the SUNDIALS context
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1))
    return (1);

  y = abstol = NULL;
  cvode_mem = NULL;

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0))
    return (1);
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0))
    return (1);

  for (int i = 0; i < NEQ; ++i) {
    NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
    NV_Ith_S(y, i) = solv_param_struc.init_state[i];
  }

  reltol = solv_param_struc.reltol;

  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0))
    return (1);

  // set error handler to Rcpp::Rcerr
  // void* ptr_to_cvode_mem = &cvode_mem;
  void *ptr_to_nothing = &cvode_mem;
  retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler, ptr_to_nothing);

  double sum_of_least_squares = 0.;

  struct usr_data my_ode_system = {ode_system, param_combi_start,
                                   params_time_combi_vec, params_cut_idx_vec};
  void *ptr_to_my_ode_system = &my_ode_system;
  retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
  if (check_retval((void *)cvode_mem, "CVodeSetUserData", 0))
    return (1);

  retval = CVodeInit(cvode_mem, wrapper_ode_system, integration_times[0], y);
  if (check_retval(&retval, "CVodeInit", 1))
    return (1);

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1))
    return (1);

  retval = CVDiag(cvode_mem);
  if (check_retval(&retval, "CVDiag", 1))
    return (1);
  double return_time;

  std::vector<double> temp_measured(init_state.size());

  for (int i = 0; i < NV_LENGTH_S(y); i++) {
    DF(0, i) = NV_Ith_S(y, i);
  }

  int counter = 1;
  for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
    return_time = integration_times[ti];
    retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
    for (int n = 0; n < NV_LENGTH_S(y); n++) {
      DF(counter, n) = NV_Ith_S(y, n);
      temp_measured[n] = hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
      if (std::isnan(temp_measured[n])) {
      } else {
        sum_of_least_squares += get_val(ecf(
          etr::Integer(static_cast<int>(integration_times.size())),
          etr::Double(NV_Ith_S(y, n)),
          etr::Double(temp_measured[n])));
      }
    }
    if (check_retval(&retval, "CVode", 1)) {
      break;
    }
    counter++;
  }

  N_VDestroy(y);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return sum_of_least_squares;
}
