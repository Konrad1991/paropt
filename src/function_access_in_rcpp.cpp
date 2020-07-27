#include "basic_functions.hpp"
#include "optimizer.hpp"
#include "solver_Rcpp_interface.hpp"
#include "paropt_types.h"
using namespace Rcpp;


#include <thread>
#include <string> //cast int to string

//' @export
//' @useDynLib paropt, .registration = TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
Rcpp::List function_access(
  std::vector<double> integration_times,
  Rcpp::XPtr<OS> fctptr, double relative_tolerance,
  std::vector<double> absolute_tolerances,
  std::string start, std::string states,
  std::string where_to_save_output_states,std::string solvertype) {


  // extract parameters
  std::vector<int> params_cut_idx_vec;
  std::vector<double> params_time_combi_vec;
  std::vector<double> param_combi_start;
  std::vector<std::string> header_parameter;


  Import_start_parameter(start, params_cut_idx_vec, params_time_combi_vec,  param_combi_start, header_parameter);

  // extract states
  std::vector<int> hs_cut_idx_vec;
  std::vector<double> hs_time_combi_vec;
  std::vector<double> hs_harvest_state_combi_vec;
  std::vector<std::string> header_states;

  Import_states(states, hs_cut_idx_vec, hs_time_combi_vec, hs_harvest_state_combi_vec, header_states);


  // extract initial values
  int tmpcount=0;

  std::vector<double> init_state ( hs_cut_idx_vec.size() );
  for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
    init_state[i] = hs_harvest_state_combi_vec[tmpcount];
    tmpcount += hs_cut_idx_vec[i];
  }

  // check absolute_tolerances
  if(static_cast<int>(init_state.size()) > absolute_tolerances.size()) {
    Rcpp::stop("\nERROR: absolute tolerances not defined for each state");
    //exit (EXIT_FAILURE);
  }

  if(static_cast<int>(init_state.size()) < absolute_tolerances.size()) {
    Rcpp::stop("\nERROR: dimension error for absolute tolerances");
    //exit (EXIT_FAILURE);
  }

  // check time in parameters vs state time
  // ============================================================
  std::vector<double>::iterator max_time_param_vector = std::max_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
  std::vector<double>::iterator min_time_param_vector = std::min_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
  std::vector<double>::iterator max_time_harvest_vector = std::max_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());
  std::vector<double>::iterator min_time_harvest_vector = std::min_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());

  bool max_time_diff_zero = double_diff_Rcpp_interface(*max_time_param_vector, 0);
  if(!max_time_diff_zero) {
    if(*max_time_param_vector > *max_time_harvest_vector) { //check || or &&
      Rcpp::warning("\nERROR: Maximum of timevector of parameter larger then corresponding timepoint of state vector");
    } else if(*max_time_param_vector < *max_time_harvest_vector) {
      Rcpp::warning("\nERROR: Maximum of timevector of parameter smaller then corresponding timepoint of state vector");
    }
    if(*min_time_param_vector > *min_time_harvest_vector) {
      Rcpp::warning("\nERROR: Minimum of timevector of parameter larger then corresponding timepoint of state vector");
    } else if(*min_time_param_vector < *min_time_harvest_vector) {
      Rcpp::warning("\nERROR: Minimum of timevector of parameter smaller then corresponding timepoint of state vector");
    }
  } else {
  // all parameter are scalare
  }

  std::vector<double> integration_time_assumption(hs_cut_idx_vec[0]);
  for(int i = 0; i < hs_cut_idx_vec[0];i++) {
    integration_time_assumption[i] = hs_time_combi_vec[i];
  }
  if(integration_times.size() > static_cast<int>(integration_time_assumption.size())) {
    Rcpp::stop("\nERROR: integration_times must not be larger than time of state input");
    //exit (EXIT_FAILURE);
  }

  bool check_entries_time;
  for(int i = 0; i < integration_times.size(); i++) {
    check_entries_time = double_diff_Rcpp_interface(integration_times[i],integration_time_assumption[i]);
    if(!check_entries_time) {
      Rcpp::warning("\nERROR: integration_times has not the same entries as the time vector of state input");
      //exit (EXIT_FAILURE);
    }
  }
  // ============================================================

  // check size of parameters either constant length = 1 or length>4 => variable
  for(size_t i = 0; i < params_cut_idx_vec.size(); i++) {
    if(params_cut_idx_vec[i] == 1 || params_cut_idx_vec[i] >=4) {
      // everything is fine. 4 values needed for spline
    } else{
      Rcpp::stop("\nERROR: neither constant nor variable parameter. Variable parameters need at least four datapoints!");
    }
  }

  // ==================================================================

  // Integration
  time_state_information_Rcpp_interface param_model;

  param_model.init_state = init_state;
  param_model.par_times = params_time_combi_vec;
  param_model.param_idx_cuts = params_cut_idx_vec;
  param_model.state_measured = hs_harvest_state_combi_vec;
  param_model.state_times = hs_time_combi_vec;
  param_model.state_idx_cut = hs_cut_idx_vec;
  param_model.integration_times = integration_times;
  param_model.reltol = relative_tolerance;
  param_model.absolute_tolerances = absolute_tolerances;

  double smsq = 0.0;

  OS ode_system = *fctptr;

  if(solvertype == "bdf") {
    smsq = solver_bdf_save_Rcpp_interface(param_combi_start, fctptr, param_model, where_to_save_output_states, header_states);
  }
  else if(solvertype == "ADAMS") {
    smsq = solver_adams_save_Rcpp_interface(param_combi_start, fctptr, param_model, where_to_save_output_states, header_states);
  } else if(solvertype == "ERK") {
    smsq = solver_erk_save_Rcpp_interface(param_combi_start, fctptr, param_model, where_to_save_output_states, header_states);
  } else if(solvertype == "ARK") {
    smsq = solver_ark_save_Rcpp_interface(param_combi_start, fctptr, param_model, where_to_save_output_states, header_states);
  } else {
    Rcpp::stop("\nERROR: Unknown solvertyp");
  }

  // ================================================================================
  std::vector<std::thread> threads(2);
  std::vector<double> results(2);

  std::vector<std::vector<double> > par(2);

  par[0].resize(param_combi_start.size());
  par[1].resize(param_combi_start.size());

  for(int i = 0; i < par.size(); i++) {
    for(int j = 0; j < param_combi_start.size(); j++) {
        par[i][j] = param_combi_start[j]*static_cast<double>(i+j);
    }
  }

  for(int i = 0; i < 2; i++) {
    std::thread temp(solver_bdf_save_Rcpp_interface, std::ref(par[i]), fctptr, param_model, std::to_string(i), header_states);
    //where_to_save_output_states, header_states);
    threads[i] = std::move(temp);
  }
  threads[0].join();
  threads[1].join();
  // ================================================================================

  return Rcpp::List::create(Rcpp::Named("Error of input-parameters:") = smsq,
                     Rcpp::Named("Solver set by user:") = solvertype,
                     Rcpp::Named("relative tolerance:") = relative_tolerance,
                     Rcpp::Named("absolute tolerance(s):") = absolute_tolerances);
}
