/* !!revision!!
improve error handling: try; catch

// better name of fct?

remove start values

Add feature to pass data.frame instead of string
*/

#include "header.hpp"
#include "modify_dataframe.hpp"
#include "optimizer_master.hpp"
#include "solver_master.hpp"
#include "paropt_types.h"

#define NA std::nan("l")

// [[Rcpp::export]]
Rcpp::List master(std::vector<double> integration_times,
                                Rcpp::XPtr<OS2> ode_sys, double relative_tolerance,
                                std::vector<double> absolute_tolerances,
                                Rcpp::DataFrame lower, Rcpp::DataFrame upper, Rcpp::DataFrame states,
                                int npop, int ngen, double error,
                                std::string solvertype) {


    // extract parameters
    //enum IMPORT_PARAMETER ret = IMPORT_PARAMETER::UNDEFINED;
    VI params_cut_idx_vec;
    VD params_time_combi_vec;
    VD param_combi_lb;
    VD param_combi_ub;
    VS header_parameter;

    ip (lower, upper, params_cut_idx_vec, params_time_combi_vec,
                     param_combi_lb, param_combi_ub, header_parameter);

    // create randomly param_combi_start vector
    VD param_combi_start(param_combi_lb.size());
    for(unsigned int i = 0; i < param_combi_lb.size(); i++) {
      double rd = arma::randu();
     param_combi_start[i] = param_combi_lb[i] + (param_combi_ub[i] - param_combi_lb[i])*rd;
    }
    // extract states
    VI hs_cut_idx_vec;
    VD hs_time_combi_vec;
    VD hs_harvest_state_combi_vec;
    VS header_states;
    Import_states(states, hs_cut_idx_vec, hs_time_combi_vec, hs_harvest_state_combi_vec, header_states);

    // check if ngen is pOS2itiv
    if(ngen <= 0) {
      Rcpp::stop("\nERROR: number of generations should be a pOS2itiv number");
    }

    // extract initial values
    int tmpcount=0;

    std::vector<double> init_state ( hs_cut_idx_vec.size() );
    for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
      init_state[i] = hs_harvest_state_combi_vec[tmpcount];
      tmpcount += hs_cut_idx_vec[i];
    }

    // check absolute_tolerances
    if(static_cast<unsigned int>(init_state.size()) > absolute_tolerances.size()) {
      Rcpp::stop("\nERROR: absolute tolerances not defined for each state");
      //exit (EXIT_FAILURE);
    }

    if(static_cast<unsigned int>(init_state.size()) < absolute_tolerances.size()) {
      Rcpp::stop("\nERROR: dimension error for absolute tolerances");
      //exit (EXIT_FAILURE);
    }

    // check time in parameters vs state time
    // ============================================================
    std::vector<double>::iterator max_time_param_vector = std::max_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
    std::vector<double>::iterator min_time_param_vector = std::min_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
    std::vector<double>::iterator max_time_harvest_vector = std::max_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());
    std::vector<double>::iterator min_time_harvest_vector = std::min_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());

    bool max_time_diff_zero = double_diff_a2a(*max_time_param_vector, 0);

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
    if(static_cast<unsigned int>(integration_times.size()) > integration_time_assumption.size()) {
      Rcpp::stop("\nERROR: integration_times must not be larger than time of state input");
      //exit (EXIT_FAILURE);
    }

    bool check_entries_time = false;
    for(unsigned int i = 0; i < integration_times.size(); i++) {
      check_entries_time = double_diff_a2a(integration_times[i],integration_time_assumption[i]);
      if(!check_entries_time) {
        Rcpp::warning("\nERROR: integration_times has not the same entries as the time vector of state input");
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

    // Test integration
    // ==================================================================
    time_state_information_a2a param_model;

    param_model.init_state = init_state;
    param_model.par_times = params_time_combi_vec;
    param_model.param_idx_cuts = params_cut_idx_vec;
    param_model.state_measured = hs_harvest_state_combi_vec;
    param_model.state_times = hs_time_combi_vec;
    param_model.state_idx_cut = hs_cut_idx_vec;
    param_model.integration_times = integration_times;
    param_model.reltol = relative_tolerance;
    param_model.absolute_tolerances = absolute_tolerances;

    // dereference XPtr
    OS2 ode_system = *ode_sys;

    // test integration
    Rcpp::NumericMatrix DF(integration_times.size(),init_state.size());
    if(solvertype == "bdf") {
      solver_bdf_save_a2a(param_combi_start, ode_system, param_model, DF);
    }
    else if(solvertype == "ADAMS") {
      solver_adams_save_a2a(param_combi_start, ode_system, param_model, DF);
    } else if(solvertype == "ERK") {
      solver_erk_save_a2a(param_combi_start, ode_system, param_model, DF);
    } else if(solvertype == "ARK") {
      solver_ark_save_a2a(param_combi_start, ode_system, param_model, DF);
    } else {
      Rcpp::stop("\nERROR: Unknown solvertyp");
    }
    // ==================================================================

    // Optimization
    // ==================================================================
    settingsPSO_a2a param_pso;

    param_pso.err_tol = error;
    param_pso.pso_n_pop = npop;
    param_pso.pso_n_gen = ngen;
    param_pso.pso_par_initial_w = 0.5; //1.0 0.1
    param_pso.pso_par_w_max = 0.90; //0.9
    param_pso.pso_par_w_min = 0.40; //0.4
    param_pso.pso_par_w_damp = 0.99;

    double (*fctptr_to_solver)(std::vector<double> &param_combi_start, OS2 ode_system, time_state_information_a2a &solv_param_struc);
    if(solvertype == "bdf") {
      fctptr_to_solver = solver_bdf_a2a;
    } else if(solvertype == "ADAMS") {
      fctptr_to_solver = solver_adams_a2a;
    } else if(solvertype == "ERK") {
      fctptr_to_solver = solver_erk_a2a;
    } else if(solvertype == "ARK") {
      fctptr_to_solver = solver_ark_a2a;
    } else {
      Rcpp::stop("\nERROR: Unknown solvertyp");
    }

    //Rcpp::List tempo = optimizer_pointer_new(param_combi_lb, param_combi_ub,
      //param_pso, param_model, fctptr_to_solver, ode_system);

    Optimizer_a2a optimizing(param_combi_lb,param_combi_ub,
                      param_pso, param_model, fctptr_to_solver ,ode_system);

    double smsq;
    smsq = optimizing.pso();
    std::vector<double> temp;
    optimizing.get_best_particle_param_values(temp);
    // ==================================================================

    // solve ODE-System with optimized parameters
    // ==================================================================
    if(solvertype == "bdf") {
    solver_bdf_save_a2a(temp, ode_system, param_model, DF);
    }
    else if(solvertype == "ADAMS") {
      solver_adams_save_a2a(temp, ode_system, param_model, DF);
    }
    else if(solvertype == "ERK") {
      solver_erk_save_a2a(temp, ode_system, param_model, DF);
    }
    else if(solvertype == "ARK") {
      solver_ark_save_a2a(temp, ode_system, param_model, DF);
    } else {
      Rcpp::stop("\nERROR: Unknown solvertyp");
    }
    Rcpp::CharacterVector CV(header_states.size()-1);
    for(unsigned int i = 1; i < header_states.size(); i++) {
      CV[i-1] = header_states[i];
    }
    colnames(DF) = CV;
    // ==================================================================

    // export parameter
    // ==================================================================
    int idxcount = 0;
    std::vector<std::vector<double> > params_export(params_cut_idx_vec.size());
    for(size_t i = 0; i < params_cut_idx_vec.size(); i++) {
      params_export[i].resize(params_cut_idx_vec[i]);
      for(size_t j = 0; j < params_export[i].size(); j++) {
        params_export[i][j] = temp[idxcount];
        idxcount++;
      }
    }

    std::vector<int>::iterator max_cut_vector = std::max_element(params_cut_idx_vec.begin(), params_cut_idx_vec.end());
    //int idx_largest_parameter = std::max_element(params_cut_idx_vec.begin(), params_cut_idx_vec.end()) - params_cut_idx_vec.begin();

    Rcpp::NumericMatrix PAROUT(*max_cut_vector, params_export.size()+1);
    int rowcounter = 0;
    while(rowcounter < *max_cut_vector) {
      for(size_t i = 0; i < params_export.size()+1; i++) {

        if(i == 0) {
          PAROUT(rowcounter, i) = integration_times[rowcounter]; //params_time_combi_vec[params_cut_idx_vec[idx_largest_parameter] + rowcounter];
        } else {
        if(rowcounter < params_cut_idx_vec[i-1]) { // <=
          PAROUT(rowcounter, i) = params_export[i-1][rowcounter];
        } else {
          PAROUT(rowcounter, i) = NA;
          }
        }
      }
      rowcounter++;
    }

    // ==================================================================

    return Rcpp::List::create(Rcpp::Named("Error of best parameters:") = smsq,
                       Rcpp::Named("Error set by user:") = error,
                       Rcpp::Named("Number of particles for PSO:") = npop,
                       Rcpp::Named("Number of generations for PSO:") = ngen,
                       Rcpp::Named("Solver set by user:") = solvertype,
                       Rcpp::Named("relative tolerance:") = relative_tolerance,
                       Rcpp::Named("absolute tolerance(s):") = absolute_tolerances,
                       Rcpp::Named("Parameter:") = PAROUT,
                     Rcpp::Named("States:") = DF);
  }
