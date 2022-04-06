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

//' Optimize parameters of ode-systems
//' @export
//' @useDynLib paropt, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @description Optimize parameters used in an ode equation in order to match values defined in the state-data.frame
//'
//' @param integration_times a vector containing the time course to solve the ode-system (see Details for more Information)
//'
//' @param ode_sys the ode-system which will be integrated by the solver (see Details for more Information). Currently the ode system has to be converted using the function 'convert'.
//'
//' @param relative_tolerance a number defining the relative tolerance used by the ode-solver.
//'
//' @param absolute_tolerances a vector containing the absolute tolerance(s) for each state used by the ode-solver.
//'
//' @param lower a data.frame containing the lower bounds for the parameters (see Details for more Information).
//'
//' @param upper a data.frame containing the upper bounds for the parameters (see Details for more Information).
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
//' @details The ode system should be of type Rcpp::XPtr<OS>. The OS is predefined in the package.
//' The function should possess the following signature: int ode(double &time, std::vector<double> &parameter, std::vector<double> &states).
//' The first entry defines the time point when the function is called.
//' The second argument defines the parameter which should be optimized. There exist two different types of parameters.
//' Parameters can be either constant or variabel. In order to calculate a variable parameter at a specific timepoint the Catmull-Rom-Spline is used.
//' This vector contains the already interpolated parameters at the specific time-point, in the same order as defined in the data.frames containing the lower- and upper-boundaries.
//' The last argument is a vector containing the states in the same order as defined in the data.frame containing the state-information.
//' Thus, it is obligatory that the state-derivates in the ode-system are in the same order defined as in the data.frame.
//' Within the function the new states have to be saved in the states-vector.
//' Please be aware that when using the approach with the Rcpp::XPtr the optimization is run in parallel. Thus, the function has to be thread-safe (among other things don't use any R Code)!
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
//' @examples
//' library(paropt)
//' # slow
//' ode <- function(t, parameter, y, ydot) {
//'
//'   a = parameter[1]
//'   b = parameter[2]
//'   c = parameter[3]
//'   d = parameter[4]
//'
//'   predator = y[1]
//'   prey = y[2]
//'
//'   ydot[1] = predator*prey*c - predator*d
//'   ydot[2] = prey*a - prey*predator*b
//' }
//'
//' # fast (but use it carefully)
//' ode <- function(t, parameter, y, ydot) {
//'
//'   a_db = at(parameter, 1)
//'   b_db = at(parameter, 2)
//'   c_db = at(parameter, 3)
//'   d_db = at(parameter, 4)
//'
//'   predator_db = at(y,1)
//'   prey_db = at(y, 2)
//'
//'   ydot[1] = predator_db*prey_db*c_db - predator_db*d_db
//'   ydot[2] = prey_db*a_db - prey_db*predator_db*b_db
//' }
//'
//' # compile
//' r <- paropt::convert(ode, verbose = TRUE)
//'
//' path <- system.file("examples", package = "paropt")
//' states <- read.table(paste(path,"/states_LV.txt", sep = ""), header = TRUE)
//'
//' # parameter
//' lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
//' ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)
//'
//' # Optimizing
//' set.seed(1)
//' df <- paropt::po(integration_times = states$time, ode_sys = r(),
//'                      relative_tolerance = 1e-6, absolute_tolerances = c(1e-8, 1e-8),
//'                      lower = lb, upper = ub, states = states,
//'                      npop = 40, ngen = 1000, error = 0.0001, solvertype = "bdf")
//'
//' par(mfrow = c(2,1))
//' plot(states$time, df$States[,1], pch = 19, type = 'l', ylab = "predator", xlab = "time", ylim = c(0, 30))
//' points(states$time, states$n1, pch = 19, col = "darkred", type = 'p')
//' legend(80, 30, legend=c("in silico", "measured"),
//'        col=c("black", "darkred"), lty=1, cex=0.8)
//' plot(states$time, df$States[,2], pch = 19, type = 'l', ylab = "prey", xlab = "time", ylim = c(0, 65))
//' points(states$time, states$n2, pch = 19, col = "darkred", type = 'p')
//' legend(80, 60, legend=c("in silico", "measured"),
//'        col=c("black", "darkred"), lty=1, cex=0.8)
// [[Rcpp::export]]
Rcpp::List po(std::vector<double> integration_times,
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
