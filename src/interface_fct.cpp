/* !!revision!!
improve error handling: try; catch
// better name of fct?
*/

#include "header.hpp"
#include "optimizer.hpp"
#include "solver.hpp"

#define NA std::nan("l")

typedef std::vector<double> VD;
typedef std::vector<int> VI;
typedef std::vector<std::vector<double> > MD;
typedef std::vector<std::vector<int> > MI;
typedef std::vector<std::string> VS;
typedef Rcpp::DataFrame DF;

//' Optimize parameters of ode-systems
//' @export
//' @useDynLib paropt, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @description Optimize parameters used in an ode equation in order to match values defined in the state-data.frame
//'
//' @param integration_times a vector containing the time course to solve the ode-system (see Details for more Information)
//'
//' @param ode_system the ode-system which will be integrated by the solver (see Details for more Information).
//'
//' @param relative_tolerance a number defining the relative tolerance used by the ode-solver.
//'
//' @param absolute_tolerances a vector containing the absolute tolerance(s) for each state used by the ode-solver.
//'
//' @param lb a data.frame containing the lower bounds for the parameters (see Details for more Information).
//'
//' @param ub a data.frame containing the upper bounds for the parameters (see Details for more Information).
//'
//' @param states a data.frame containing the predetermined course of the states (see Details for more Information).
//'
//' @param npop a number defining the number of particles used by the Particle Swarm Optimizer.
//'
//' @param ngen a number defining the number of generations the Particle Swarm Optimizer (PSO) should run.
//'
//' @param error a number defining a sufficient small error. When the PSO reach this value optimization is stopped.
//'
//' @param solvertype a string defining the type of solver which should be used (bdf, ADAMS, ERK or ARK. see Details for more Information).
//'
//' @details The vector containing the time course to solve the ode-system should contain
//' the same entries as the time vector in the state-data.frame (it can be also be a different variable instead of time).
//'
//' @details The ode system should be a Rcpp-function with a specific signature Rcpp::NumericVector ode(double time, std::vector<double> parameter, Rcpp::NumericVector states).
//' The first entry defines the time point when the function is called.
//' The second argument defines the parameter which should be optimized. There exist two different types of parameters.
//' Parameters can be either constant or variabel. In order to calculate a variable parameter at a specific timepoint the Catmull-Rom-Spline is used.
//' This vector contains the already interpolated parameters at the specific time-point, in the same order as defined in the data.frames containing the lower- and upper-boundaries.
//' The last argument is a vector containing the states in the same order as defined in the data.frame containing the state-information.
//' Thus, it is obligatory that the state-derivates in the ode-system are in the same order defined as in the data.frame.
//' Furthermore, it is mandatory that the function return a Rcpp::NumericVector with the same dimension as the input vector containing the states.
//' The resulting vector has to contain the right hand side of the ode-system.
//'
//' @details For constant parameters use only the first row (below the headers) if other parameters are variable use “NA“ in the following rows for the constant parameters.
//' @details For variable parameters at least four points are needed. If a variable parameter is not available at every time point use “NA“ instead. 
//'
//' @details The two data.frames containg lower and upper-boundaries need the parameter in the same order. 
//'
//' @details The data.frame containing the state information should hold the time course in the first column.
//' The header-name time is compulsory. The following columns contain the states. Take care that the states are in the same order defined in the ode system.
//' If a state is not available use “NA“. This is possible for every time points except the first one.
//' The ode solver need a start value for each state which is extracted from the first row of this file (below the headers).
//'
//' @details The error between the solver output and the measured states is the sum of the absolute differences divided by the number of time points.
//' It is crucial that the states are in the same order in the text file cointaining the state-information and in the ode-system to compare the states correctly!
//'
//' @details For solving the ode system the SUNDIALS Software is used (https://computing.llnl.gov/projects/sundials).
//' The last argument defines the solver-type which is used during optimization:
//' “bdf“,  “ADAMS“, “ERK“ or “ARK“. bdf = Backward Differentiation Formulas, ADAMS = Adams-Moulton, ERK = explicite Runge-Kutta and ARK = implicite Runge-Kutta.
//' All solvers are used in the NORMAL-Step method in a for-loop using the time-points defined in the text-file containing the states as output-points.
//' The bdf- and ARK-Solver use the SUNLinSol_Dense as linear solver. Notably here is that for the ARK-Solver the ode system is fully implicit solved (not only part of it).
//'
//' Examples can be found in the vignette.
// [[Rcpp::export]]
  Rcpp::List optimizer(Rcpp::NumericVector integration_times,
                     SEXP ode_system, double relative_tolerance,
                     Rcpp::NumericVector absolute_tolerances,
                       Rcpp::DataFrame lb, Rcpp::DataFrame ub, Rcpp::DataFrame states,
                     int npop, int ngen, double error,
                 std::string solvertype) {


     // extract parameters

    //enum IMPORT_PARAMETER ret = IMPORT_PARAMETER::UNDEFINED;
    VI params_cut_idx_vec;
    VD params_time_combi_vec;
    VD param_combi_lb;
    VD param_combi_ub;
    VS header_parameter;

    ip (lb, ub, params_cut_idx_vec, params_time_combi_vec,
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

    // check if ngen is positiv
    if(ngen <= 0) {
      Rcpp::stop("\nERROR: number of generations should be a positiv number");
    }

    // extract initial values
    int tmpcount=0;

    std::vector<double> init_state ( hs_cut_idx_vec.size() );
    for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
      init_state[i] = hs_harvest_state_combi_vec[tmpcount];
      tmpcount += hs_cut_idx_vec[i];
    }

    // check absolute_tolerances
    if(static_cast<int>(init_state.size()) > absolute_tolerances.length()) {
      Rcpp::stop("\nERROR: absolute tolerances not defined for each state");
      //exit (EXIT_FAILURE);
    }

    if(static_cast<int>(init_state.size()) < absolute_tolerances.length()) {
      Rcpp::stop("\nERROR: dimension error for absolute tolerances");
      //exit (EXIT_FAILURE);
    }

    // check time in parameters vs state time
    // ============================================================
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
    if(static_cast<unsigned int>(integration_times.length()) > integration_time_assumption.size()) {
      Rcpp::stop("\nERROR: integration_times must not be larger than time of state input");
      //exit (EXIT_FAILURE);
    }

    bool check_entries_time = false;
    for(int i = 0; i < integration_times.length(); i++) {
      check_entries_time = double_diff(integration_times[i],integration_time_assumption[i]);
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

    // checks of ODE-System
    // ==================================================================
    if(TYPEOF(ode_system) != CLOSXP) {
      Rcpp::stop("\nERROR: type of odesystem should be closure");
      //exit (EXIT_FAILURE);
    }

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

    // Test integration
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

    // test integration
    Rcpp::NumericMatrix DF(integration_times.size(),init_state.size());
    if(solvertype == "bdf") {
      solver_bdf_save(param_combi_start, ode_system, param_model, DF);
    }
    else if(solvertype == "ADAMS") {
      solver_adams_save(param_combi_start, ode_system, param_model, DF);
    } else if(solvertype == "ERK") {
      solver_erk_save(param_combi_start, ode_system, param_model, DF);
    } else if(solvertype == "ARK") {
      solver_ark_save(param_combi_start, ode_system, param_model, DF);
    } else {
      Rcpp::stop("\nERROR: Unknown solvertyp");
    }
    // ==================================================================

    // Optimization
    // ==================================================================
    settingsPSO param_pso;

    param_pso.err_tol = error;
    param_pso.pso_n_pop = npop;
    param_pso.pso_n_gen = ngen;
    param_pso.pso_par_initial_w = 0.5; //1.0 0.1
    param_pso.pso_par_w_max = 0.90; //0.9
    param_pso.pso_par_w_min = 0.40; //0.4
    param_pso.pso_par_w_damp = 0.99;

    double (*fctptr_to_solver)(std::vector<double> &param_vec, SEXP ode_system, time_state_information &model);
    if(solvertype == "bdf") {
      fctptr_to_solver = solver_bdf;
    } else if(solvertype == "ADAMS") {
      fctptr_to_solver = solver_adams;
    } else if(solvertype == "ERK") {
      fctptr_to_solver = solver_erk;
    } else if(solvertype == "ARK") {
      fctptr_to_solver = solver_ark;
    } else {
      Rcpp::stop("\nERROR: Unknown solvertyp");
    }

    Optimizer optimizing(param_combi_lb,param_combi_ub,param_pso,param_model, fctptr_to_solver ,ode_system);

    double smsq;
    smsq = optimizing.pso();
    std::vector<double> temp;
    optimizing.get_best_particle_param_values(temp);
    // ==================================================================

    // solve ODE-System with optimized parameters
    // ==================================================================
    if(solvertype == "bdf") {
    solver_bdf_save(temp, ode_system, param_model, DF);
    }
    else if(solvertype == "ADAMS") {
      solver_adams_save(temp, ode_system, param_model, DF);
    }
    else if(solvertype == "ERK") {
      solver_erk_save(temp, ode_system, param_model, DF);
    }
    else if(solvertype == "ARK") {
      solver_ark_save(temp, ode_system, param_model, DF);
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
    int idx_largest_parameter = std::max_element(params_cut_idx_vec.begin(), params_cut_idx_vec.end()) - params_cut_idx_vec.begin();

    Rcpp::NumericMatrix PAROUT(*max_cut_vector, params_export.size()+1);
    int rowcounter = 0;
    while(rowcounter < *max_cut_vector) {
      for(size_t i = 0; i < params_export.size()+1; i++) {

        if(i == 0) {
          PAROUT(rowcounter, i) = params_time_combi_vec[params_cut_idx_vec[idx_largest_parameter] + rowcounter];
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
