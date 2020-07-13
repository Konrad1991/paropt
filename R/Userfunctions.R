String_Mod <- function(file) {
  path <- getwd()
  # Test if user specify full path or not
  file_mod <- gsub(path,"", file)
  file_mod2 <- gsub("/", "", file) # test if other directory is used. 
  if(nchar(file) != nchar(file_mod) || nchar(file) != nchar(file_mod2)) {
    print("Use string from user. Do not modify")
    ret <- file
  } else if(nchar(file)) {
    print("Modify string. Using path from setwd")
    path <- paste(path, "/", sep = "")
    file <- paste(path, file, sep = "")
    ret <- file
  }
  return(ret)
}

Identify_datatyp <- function(file) {
  file_extension <- substr(file, nchar(file) - 3 + 1, nchar(file))
  return(file_extension)
}

#' Optimize parameters of ode-systems
#' @export
#' @useDynLib paropt, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @param integration_times a vector containing the time course to solve the ode-system (see Details for more Information)
#' @param ode_system the ode-system which will be integrated by the solver (see Details for more Information). 
#' @param relative_tolerance a number defining the relative tolerance used by the ode-solver.
#' @param absolute_tolerances a vector containing the absolute tolerance(s) for each state used by the ode-solver.
#' @param start a string-path to a tab seperated text-file containing startvalues for the parameters (see Details for more Information). 
#' @param lower a string-path to a tab seperated text-file containing the lower bounds for the parameters (see Details for more Information).
#' @param upper a string-path to a tab seperated text-file containing the upper bounds for the parameters (see Details for more Information).
#' @param states a string-path to a tab seperated text-file containing the measured states (see Details for more Information). 
#' @param npop a number defining the number of particles used by the Particle Swarm Optimizer.
#' @param ngen a number defining the number of generations the Particle Swarm Optimizer (PSO) should run.
#' @param error a number defining a sufficient small error. When the PSO reach this value optimization is stopped.
#' @param where_to_save_output_states a string-path defining a name for a textfile where the result of the integration of the states is saved. Using the previously optimized parameter for this integration.
#' @param where_to_save_output_parameter a string-path defining a name for a textfile where the optimized parameters are saved. 
#' @param solvertype a string defines the type of solver which should be used (bdf, ADAMS, ERK or ARK. see Details for more Information).
#' @description optimizer() finds parameters of an ode-system to match measured states.
#' @details The vector containing the time course to solve the ode-system should contain the same entries as the time vector in the text file containing the states (of course it can be also be a different variable instead of time). It is possible that the vector is shorter than the time vector defined in the state-file in order to optimize only a part of the problem.  
#' @details The ode system should be a Rcpp-function with a specific signature. The name of the function is free to choose. The following parameters have to be passed: a double t, a std::vector<double> params, and a Rcpp::NumericVector y.
#' @details The first entry defines the time point when the function is called.
#' @details The second argument defines the parameter which should be optimized. There exist two different types of parameters. Parameters can be either constant or variabel. In order to calculate a variable parameter at a specific timepoint the Catmull-Rom-Spline is used. This vector contains the already splined parameters, in the same order as defined in the text-files containing the start-values and the lower- and upper-boundaries. 
#' @details The last argument is a vector containing the states in the same order as defined in the text-file containing the state-information. Thus, it is obligatory that the state-derivates in the ode-system are in the same order defined as in the text-file. 
#' @details Furthermore, it is mandatory that the function return a Rcpp::NumericVector with the same dimension as the input vector containing the states. Naturally, the vector should contain the right hand side of the ode-system. 
#' @details The files containing the start values (used to test integration) for the parameter, the lower- and upper-boundaries must have the following layout. In the first column the time is defined. In the following columns the parameters are defined. Consider that the parameter order is the same as used in the ode-system. 
#' @details For constant parameters use only the first row (below the headers) if other parameters are variable use “NA“ in the following rows for the constant parameters. 
#' @details For variable parameters at least four points are needed. If a variable parameter is not available at every time point use “NA“ instead. . 
#' @details The three files start-values, lower and upper-boundaries need the parameter in the same order. The particles are randomly created within the lower and upper boundary. 
#' @details The file containing the state information should contain in the first column the time. The header-name time is compulsory. The following columns contain the states. Take care that the state order is the same as defined in the ode system. If a state is not available use “NA“. This is possible for every time points except the first one. The ode solver need a start value for each state which is extracted from the first row of this file (below the headers).
#' @details The error between the solver output and the measured states is the sum of the absolute differences divided by the number of time points. It is crucial that the states are in the same order in the text file cointaining the state-information and in the ode-system to compare the states correctly!
#' @details For solving the ode system the SUNDIALS Software is used (https://computing.llnl.gov/projects/sundials). The last argument defines the solver-type which is used during optimization: “bdf“,  “ADAMS“, “ERK“ or “ARK“. bdf = Backward Differentiation Formulas, ADAMS = Adams-Moulton, ERK = explicite Runge-Kutta and ARK = implicite Runge-Kutta. All solvers are used in the NORMAL-Step method in a for-loop using the time-points defined in the text-file containing the states as output-points. The bdf- and ARK-Solver use the SUNLinSol_Dense as linear solver. Notably here is that for the ARK-Solver the ode system is fully implicit solved (not only part of it).
#' @example /inst/examples/optimizer_examples.r 
optimizer <- function(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, 
                      lower, upper, states, npop, ngen, error, where_to_save_output_states, where_to_save_output_parameter, solvertype) {
  start <- String_Mod(start)
  lower <- String_Mod(lower)
  upper <- String_Mod(upper)
  states <- String_Mod(states)
  where_to_save_output_states <- String_Mod(where_to_save_output_states)
  where_to_save_output_parameter <- String_Mod(where_to_save_output_parameter)
  interface_function(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, 
                     lower, upper, states, npop, ngen, error, where_to_save_output_states, where_to_save_output_parameter, solvertype)
}


#' Solves ode-system and compare result to measured states
#' @export
#' @useDynLib paropt, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @param integration_times a vector containing the time course to solve the ode-system (see Details for more Information)
#' @param ode_system the ode-system which will be integrated by the solver (see Details for more Information). 
#' @param relative_tolerance a number defining the relative tolerance used by the ode-solver.
#' @param absolute_tolerances a vector containing the absolute tolerance(s) for each state used by the ode-solver.
#' @param start a string-path to a tab seperated text-file containing values for the parameters (see Details for more Information). 
#' @param states a string-path to a tab seperated text-file containing the measured states (see Details for more Information). 
#' @param where_to_save_output_states a string-path defining a name for a textfile where the result of the integration of the states is saved. Using the previously optimized parameter in this integration.
#' @param solvertype a string defines the type of solver which should be used (bdf, ADAMS, ERK or ARK. see Details for more Information).
#' @description ode_solving solves ode system and calculates error between solved and measured values. 
#' @details The vector containing the time course to solve the ode-system should contain the same entries as the time vector in the text file containing the states (of course it can be also be a different variable instead of time). It is possible that the vector is shorter than the time vector defined in the state-file in order to optimize only a part of the problem.  
#' @details The ode system should be a Rcpp-function with a specific signature. The name of the function is free to choose. The following parameters have to be passed: a double t, a std::vector<double> params, and a Rcpp::NumericVector y.
#' @details The first entry defines the time point when the function is called.
#' @details The second argument defines the parameter which should be optimized. There exist two different types of parameters. Parameters can be either constant or variabel. In order to calculate a variable parameter at a specific timepoint the Catmull-Rom-Spline is used. This vector contains the already splined parameters, in the same order as defined in the text-files containing the start-values and the lower- and upper-boundaries. 
#' @details The last argument is a vector containing the states in the same order as defined in the text-file containing the state-information. Thus, it is obligatory that the state-derivates in the ode-system are in the same order defined as in the text-file. 
#' @details Furthermore, it is mandatory that the function return a Rcpp::NumericVector with the same dimension as the input vector containing the states. Naturally, the vector should contain the right hand side of the ode-system. 
#' @details The file containing the start values for the parameter must have the following layout. In the first column the time is defined. In the following columns the parameters are defined. Consider that the parameter order is the same as used in the ode-system. 
#' @details For constant parameters use only the first row (below the headers) if other parameters are variable use “NA“ in the following rows for the constant parameters. 
#' @details For variable parameters at least four points are needed. If a variable parameter is not available at every time point use “NA“ instead. 
#' @details Furthermore, it is notably that the time of the parameter should be within the time vector defined in the text-file containing the state information. 
#' @details The file containing the state information should contain in the first column the time. The header-name time is compulsory. The following columns contain the states. Take care that the state order is the same as defined in the ode system. If a state is not available use “NA“. This is possible for every time points except the first one. The ode solver need a start value for each state which is extracted from the first row of this file (below the headers).
#' @details The error between the solver output and the measured states is the sum of the absolute differences divided by the number of time points. It is crucial that the states are in the same order in the text file cointaining the state-information and in the ode-system to compare the states correctly!
#' @details For solving the ode system the SUNDIALS Software is used (https://computing.llnl.gov/projects/sundials). The last argument defines the solver-type which is used during optimization: “bdf“,  “ADAMS“, “ERK“ or “ARK“. bdf = Backward Differentiation Formulas, ADAMS = Adams-Moulton, ERK = explicite Runge-Kutta and ARK = implicite Runge-Kutta. All solvers are used in the NORMAL-Step method in a for-loop using the time-points defined in the text-file containing the states as output-points. The bdf- and ARK-Solver use the SUNLinSol_Dense as linear solver. Notably here is that for the ARK-Solver the ode system is fully implicit solved (not only part of it).
#' @example /inst/examples/ode_solving_examples.r 
ode_solving <- function(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, states,
                      where_to_save_output_states, solvertype) {
  start <- String_Mod(start)
  states <- String_Mod(states)
  where_to_save_output_states <- String_Mod(where_to_save_output_states)
  solve_ode_system(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, 
                     states, where_to_save_output_states, solvertype)
}



