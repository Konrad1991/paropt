#include "solver.hpp"

// [[Rcpp::export]]
Rcpp::List wrapper_solver(vd& init_state,
                             vd& par_times,
                             vi& param_idx_cuts,
                             vd& parameter_vec,
                             vd& state_measured,
                             vi& state_idx_cuts,
                             vd& integration_times,
                             double reltol,
                             vd& absolute_tolerances,
                             Rcpp::XPtr<OS> fct,
                             int solvertype) {


  // add parameter to struct
  time_state_information tsi;

  tsi.init_state = init_state;
  tsi.par_times = par_times;
  tsi.param_idx_cuts = param_idx_cuts;
  tsi.state_measured = state_measured;
  tsi.state_idx_cut = state_idx_cuts;
  tsi.integration_times = integration_times;
  tsi.reltol = reltol;
  tsi.absolute_tolerances= absolute_tolerances;

  OS ode = *fct;

  // define solver
  solver_ptr_save save_fct;
  if(solvertype == 1) {
    save_fct = solver_bdf_save;
  } else if(solvertype == 2){
    save_fct = solver_adams_save;
  }
  Rcpp::NumericMatrix df(integration_times.size(), init_state.size());

  double error = save_fct(parameter_vec, ode, tsi, df);
  Rcpp::List L = Rcpp::List::create(error, df);
  return L;
}
