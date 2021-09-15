library(testthat)
library(paropt)

context("test interface function")
path <- system.file("tests/testthat/files", package = "paropt")

states <- read.table(paste(path, "/states_LV.txt", sep = "/") , header = TRUE)
names(states) <- c("time", "prey", "predator")
Rcpp::cppFunction('Rcpp::NumericVector ode_system(double t, std::vector<double> params,
                               Rcpp::NumericVector states) {
                               
  Rcpp::NumericVector states_deriv(2); //two states in the example above (prey and predator)
  
  // store parameters in variables
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  
  // store states in variables
  double predator = states[0]; 
  double prey = states[1]; 
  
  // the actual ode-system            
  double ddtpredator = states_deriv[0] = predator*c*prey - predator*d;
  double ddtprey = states_deriv[1] = prey*a - prey*b*predator;
  
  return states_deriv;
}')

lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)

test_that("Test Interface Function abs tolerances", {
  expect_error(paropt:::optimizer(integration_times = states$time,
                                           ode_system = ode_system,
                                           relative_tolerance = 1e-6,
                                           absolute_tolerances = 1e-6,
                                           lb = lb, ub = ub, states = states, npop = 40, ngen = 1000,
                                           error = 0.001, solvertype = "bdf"), 
               "ERROR: absolute tolerances not defined for each state")
})


test_that("Test Interface Function correct run", {
  
  path <- system.file("tests/testthat/files", package = "paropt")
  
  states <- read.table(paste(path, "/states_LV.txt", sep = "/") , header = TRUE)
  names(states) <- c("time", "prey", "predator")
  Rcpp::cppFunction('Rcpp::NumericVector ode_system(double t, std::vector<double> params,
                               Rcpp::NumericVector states) {
                               
  Rcpp::NumericVector states_deriv(2); //two states in the example above (prey and predator)
  
  // store parameters in variables
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  
  // store states in variables
  double predator = states[0]; 
  double prey = states[1]; 
  
  // the actual ode-system            
  double ddtpredator = states_deriv[0] = predator*c*prey - predator*d;
  double ddtprey = states_deriv[1] = prey*a - prey*b*predator;
  
  return states_deriv;
}')
  
  lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
  ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)
  
  res <- paropt:::optimizer(integration_times = states$time,
                     ode_system = ode_system,
                     relative_tolerance = 1e-6,
                     absolute_tolerances = c(1e-8,1e-8), 
                     lb = lb, ub = ub, states = states, npop = 40, ngen = 1000,
                     error = 0.001, solvertype = "bdf")
  true_res <- round(c(0, 1.100091, 0.4000309, 0.100004, 0.3999953), digits = 1)
  sim_res <- round(c(res[[8]]), digits = 1)
  
  testthat::expect_equal(true_res, sim_res)
})



test_that("Test Interface Function correct run", {
  path <- system.file("tests/testthat/files", package = "paropt")

  Rcpp::sourceCpp(paste(path, "/ode.cpp", sep = ""))

  states <- read.table(paste(path, "/states_LV.txt", sep = "/") , header = TRUE)
  names(states) <- c("time", "prey", "predator")

  lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
  ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)
  set.seed(1)
  res <- paropt:::optimizer_pointer(integration_times = states$time,
                                    ode_sys = test_optimization(),
                                    relative_tolerance = 1e-6,
                                    absolute_tolerances = rep(1e-8, 2),
                                    lower = lb, upper = ub, states = states, npop = 40, ngen = 2000,
                                    error = 0.001, solvertype = "ADAMS")
  true_res <- round(c(0, 1.100091, 0.4000309, 0.100004, 0.3999953), digits = 2)
  sim_res <- round(c(res[[8]]), digits = 2)

  testthat::expect_equal(true_res, sim_res)
})
