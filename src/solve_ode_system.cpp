#include "header.hpp"
#include "optimizer.hpp"
#include "solver.hpp"


  // [[Rcpp::export]]
  Rcpp::List solve_ode_system(Rcpp::NumericVector integration_times,
                     SEXP ode_system, double relative_tolerance,
                     Rcpp::NumericVector absolute_tolerances,
                      std::string start, std::string states,
                     std::string where_to_save_output_states,std::string solvertype) {

    std::vector<int> params_cut_idx_vec;
    std::vector<double> params_time_combi_vec;
    std::vector<double> param_combi_start;
    std::vector<std::string> header_parameter;

    Import_start_parameter(start, params_cut_idx_vec, params_time_combi_vec,  param_combi_start, header_parameter);

    std::vector<int> hs_cut_idx_vec;
    std::vector<double> hs_time_combi_vec;
    std::vector<double> hs_harvest_state_combi_vec;
    std::vector<std::string> header_states;

    Import_states(states, hs_cut_idx_vec, hs_time_combi_vec, hs_harvest_state_combi_vec, header_states);

    int tmpcount=0;

    std::vector<double> init_state ( hs_cut_idx_vec.size() );
    for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
      init_state[i] = hs_harvest_state_combi_vec[tmpcount];
      tmpcount += hs_cut_idx_vec[i];
    }

    if(static_cast<int>(init_state.size()) > absolute_tolerances.length()) {
      Rcpp::stop("\nERROR: absolute tolerances not defined for each state");
      //exit (EXIT_FAILURE);
    }

    if(static_cast<int>(init_state.size()) < absolute_tolerances.length()) {
      Rcpp::stop("\nERROR: dimension error for absolute tolerances");
      //exit (EXIT_FAILURE);
    }

    std::vector<double>::iterator max_time_param_vector = std::max_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
    std::vector<double>::iterator min_time_param_vector = std::min_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
    std::vector<double>::iterator max_time_harvest_vector = std::max_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());
    std::vector<double>::iterator min_time_harvest_vector = std::min_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());

    bool max_time_diff_zero = double_diff(*max_time_param_vector, 0);
    //bool max_time_param_vs_max_time_harvest = double_diff(*max_time_param_vector, *max_time_harvest_vector);
    //bool min_time_param_vs_min_time_harvest = double_diff(*min_time_param_vector, *min_time_harvest_vector);
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
    if(integration_times.length() > static_cast<int>(integration_time_assumption.size())) {
      Rcpp::stop("\nERROR: integration_times must not be larger than time of state input");
      //exit (EXIT_FAILURE);
    }

    bool check_entries_time;
    for(int i = 0; i < integration_times.length(); i++) {
      check_entries_time = double_diff(integration_times[i],integration_time_assumption[i]);
      if(!check_entries_time) {
        Rcpp::warning("\nERROR: integration_times has not the same entries as the time vector of state input");
        //exit (EXIT_FAILURE);
      }
    }

    for(size_t i = 0; i < params_cut_idx_vec.size(); i++) {
      if(params_cut_idx_vec[i] == 1 || params_cut_idx_vec[i] >=4) {
        // everything is fine. 4 values needed for spline
      } else{
        Rcpp::stop("\nERROR: neither constant nor variable parameter. Variable parameters need at least four datapoints!");
      }
    }

    if(TYPEOF(ode_system) != CLOSXP) {
      Rcpp::stop("\nERROR: type of odesystem should be closure");
      //exit (EXIT_FAILURE);
    }

    // ==================================================================
    Rcpp::NumericVector new_states(init_state.size());
    realtype time_param_sort = params_time_combi_vec[0];
    std::vector<double> parameter_input;
    double time = params_time_combi_vec[0];

    params_sort(time_param_sort, parameter_input, param_combi_start, params_time_combi_vec, params_cut_idx_vec);

    Rcpp::NumericVector current_states(hs_cut_idx_vec.size());
    for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
      current_states[i] = init_state[i];
    }

    Rcpp::Function odes = ode_system;
    Rcpp::Function fo("formals");
    Rcpp::List test_num_arguments = fo(odes);
    if(test_num_arguments.length() != 3) {
      Rcpp::stop("\nERROR: odesystem should only accept three arguments");
      //exit (EXIT_FAILURE);
    }

    try {odes(time, parameter_input, current_states); } catch (...) {
    Rcpp::stop("\nERROR: odesystem cannot be called. May be wrong types of arguments (double, std::vector<double>, Rcpp::NumericVector)?");
    //exit (EXIT_FAILURE);
    }

    bool output_numeric = false;
    switch(TYPEOF(odes(time, parameter_input, current_states))) {
      case REALSXP: {
        //Rcpp::Rcerr << "correct output" << std::endl;
        output_numeric = true;
      }
    }

    if(output_numeric == true) {
      try {new_states = odes(time, parameter_input, current_states); } catch (...) {
        Rcpp::stop("\nERROR: output of odesystem is wrong! Has to be Rcpp::NumericVector");
        //exit (EXIT_FAILURE);
      }
    } else {
      Rcpp::stop("\nERROR: output of odesystem is wrong! Has to be Rcpp::NumericVector");
    }

    if(new_states.length() != static_cast<int>(init_state.size())) {
      Rcpp::stop("\nERROR: output of odesystem is wrong! Has to be same size as number of states");
      //exit (EXIT_FAILURE);
    }
    // ==================================================================

    time_state_information param_model;

    param_model.init_state = init_state;
    param_model.par_times = params_time_combi_vec;
    param_model.param_idx_cuts = params_cut_idx_vec;
    param_model.state_measured = hs_harvest_state_combi_vec;
    param_model.state_times = hs_time_combi_vec;
    param_model.state_idx_cut = hs_cut_idx_vec;
    param_model.integration_times = integration_times;
    param_model.reltol = relative_tolerance;
    param_model.absolute_tolerances = absolute_tolerances;

    double smsq;
    // test integration
    if(solvertype == "bdf") {
      smsq = solver_bdf_save(param_combi_start, ode_system, param_model, where_to_save_output_states, header_states);
    }
    else if(solvertype == "ADAMS") {
      smsq = solver_adams_save(param_combi_start, ode_system, param_model, where_to_save_output_states, header_states);
    } else if(solvertype == "ERK") {
      smsq = solver_erk_save(param_combi_start, ode_system, param_model, where_to_save_output_states, header_states);
    } else if(solvertype == "ARK") {
      smsq = solver_ark_save(param_combi_start, ode_system, param_model, where_to_save_output_states, header_states);
    } else {
      Rcpp::stop("\nERROR: Unknown solvertyp");
    }

    return Rcpp::List::create(Rcpp::Named("Error of input-parameters:") = smsq,
                       Rcpp::Named("Solver set by user:") = solvertype,
                       Rcpp::Named("relative tolerance:") = relative_tolerance,
                       Rcpp::Named("absolute tolerance(s):") = absolute_tolerances);
  }
