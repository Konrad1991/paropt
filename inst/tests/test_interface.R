library(testthat)
library(ParOpt2)

context("test interface function")
path <- system.file("tests/testthat/files", package = "ParOpt2")

Rcpp::cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')
int_time <- 1:5

test_that("Test Interface Function abs tolerances", {
  expect_error(ParOpt2:::test_interface_fct(int_time, add, 1e-6, c(1e-6), paste(path, "/par.txt", sep = ""),
                                     paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                     paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: absolute tolerances not defined for each state")
})

test_that("Test Interface Function abs tolerances", {
  expect_error(ParOpt2:::test_interface_fct(int_time, add, 1e-6, c(1e-6, 1, 1), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: dimension error for absolute tolerances")
})

test_that("Test Interface Function time vector comparison states parameter", {
  expect_error(ParOpt2:::test_interface_fct(int_time, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states_wrong_time.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: Maximum of timevector of parameter smaller then corresponding timepoint of state vector")
})

test_that("Test Interface Function time vector comparison states parameter", {
  expect_error(ParOpt2:::test_interface_fct(int_time, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states_wrong_time2.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: Maximum of timevector of parameter larger then corresponding timepoint of state vector")
})

test_that("Test Interface Function integration time", {
  expect_error(ParOpt2:::test_interface_fct(1:6, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: integration_times must not be larger than time of state input")
})

test_that("Test Interface Function integration time", {
  expect_error(ParOpt2:::test_interface_fct(6:10, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: integration_times has not the same entries as the time vector of state input")
})

test_that("Test Interface Function integration time", {
  expect_error(ParOpt2:::test_interface_fct(6:10, add, 1e-6, c(1e-6), paste(path, "/test_integration_time.txt", sep = ""),
                                            paste(path, "/test_integration_time_lb.txt", sep = ""),paste(path, "/test_integration_time_ub.txt", sep = ""),
                                            paste(path, "/test_integration_time_states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: integration_times has not the same entries as the time vector of state input")
})


test_that("Test Interface Function integration time", {
  expect_error(ParOpt2:::test_interface_fct(1:6, add, 1e-6, c(1e-6, 1e-6), paste(path, "/test_strange_params.txt", sep = ""),
                                            paste(path, "/test_strange_params_lb.txt", sep = ""),paste(path, "/test_strange_params_ub.txt", sep = ""),
                                            paste(path, "/test_strange_params_states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: neither constant nor variable parameter. Variable parameters need at least four datapoints!")
})

test_that("Test Interface Function ode system", {
  expect_error(ParOpt2:::test_interface_fct(1:5, 1, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: type of odesystem should be closure")
})

Rcpp::cppFunction('int add(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector z, int v) {
  int sum = x[0]+  y[0] + z[0];
  return sum;
}')

test_that("Test Interface Function ode system number of arguments", {
  expect_error(ParOpt2:::test_interface_fct(1:5, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: odesystem should only accept three arguments")
})

Rcpp::cppFunction('int add(Rcpp::NumericVector x, Rcpp::NumericVector y, int g) {
  int sum = x[0]+  y[0] ;
  return sum;
}')

test_that("Test Interface Function ode system callable?", {
  expect_error(ParOpt2:::test_interface_fct(1:5, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: odesystem cannot be called. May be wrong types of arguments (double, Rcpp::NumericVector, Rcpp::NumericVector)?")
})

Rcpp::cppFunction('Rcpp::NumericVector add(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector z) {
  Rcpp::NumericVector ret(1);
  ret[0] = x[0]+  y[0] ;
  return ret;
}')

test_that("Test Interface Function ode system output size", {
  expect_error(ParOpt2:::test_interface_fct(1:5, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: output of odesystem is wrong! Has to be same size as number of states")
})

Rcpp::cppFunction('std::string add(double x, Rcpp::NumericVector y, Rcpp::NumericVector z) {
  Rcpp:NumericVector sum(y.length());
  sum[0] = x +  y[0];
  sum[1] = x + y[1]*z[1];
  return "test";
}')

test_that("Test Interface Function ode system output", {
  expect_error(ParOpt2:::test_interface_fct(1:5, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: output of odesystem is wrong! Has to be Rcpp::NumericVector")
})

Rcpp::cppFunction('Rcpp::NumericVector add(double t, Rcpp::NumericVector p, Rcpp::NumericVector y) {
  Rcpp::NumericVector ydot(3);
  ydot[0] = y[0];
  ydot[1] = y[1];
  ydot[2] = y[2];

  return ydot;
}')

test_that("Test Interface Function ode system output size", {
  expect_error(ParOpt2:::test_interface_fct(1:5, add, 1e-6, c(1e-6, 1e-6), paste(path, "/par.txt", sep = ""),
                                  paste(path, "/par_lb.txt", sep = ""),paste(path, "/par_ub.txt", sep = ""),
                                  paste(path, "/states.txt", sep = ""), 40, 1000, 0.1, "out.txt", "parout.txt"),
               "ERROR: output of odesystem is wrong! Has to be same size as number of states")
})


