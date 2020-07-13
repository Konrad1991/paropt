library(testthat)
library(ParOpt2)

context("test Solver")
path <- system.file("tests/testthat/files", package = "ParOpt2")
# # Test Solver =======================================================
Rcpp::cppFunction('Rcpp::NumericVector ode(double t, std::vector<double> params, Rcpp::NumericVector y) {
  Rcpp::NumericVector ydot(y.length());

  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];

  double n1 = y[0];
  double n2 = y[1];

  ydot[0] = n1*c*n2 - n1*d;
  ydot[1] = n2*a - n2*b*n1;

  return ydot;
}')

out <- read.table(paste(path, "/outputdeSolve.txt", sep = ""), header = T)
correct <- list(out[,2], out[,3])

test_that("check ode solving", {
  expect_equal(ParOpt2:::test_solver(seq(0, 100, pi/2), ode, 1e-6, c(1e-8, 1e-8), paste(path, "/paramstartLV.txt", sep = ""),
                                     paste(path, "/lbLV.txt", sep = ""),paste(path, "/ubLV.txt", sep = ""),
                                     paste(path, "/states_LV.txt", sep = ""), solvertype = "bdf"), correct,
               tolerance = 0.001)
})

test_solving_with_start_values <- function(solvertype) {
  out <- ParOpt2:::test_solve_ode_system(seq(0, 100, pi/2), ode, 1e-6, c(1e-8, 1e-8), paste(path, "/paramstartLV.txt", sep = ""),
                                  paste(path, "/states_LV.txt", sep = ""), solvertype)
  return(out)
}
test_that("check ode solving with start parameter values as text file", {
  expect_equal(test_solving_with_start_values("bdf"), correct, tolerance = 0.001)
  expect_equal(test_solving_with_start_values("ADAMS"), correct, tolerance = 0.001)
  expect_equal(test_solving_with_start_values("ERK"), correct, tolerance = 0.001)
  expect_equal(test_solving_with_start_values("ARK"), correct, tolerance = 0.001)
})

test_that("check error calculation", {
  expect_equal(ParOpt2:::test_error_calculation(paste(path,"/error.txt", sep = ""),paste(path,"/measured.txt", sep ="")), 12)
})

