library(testthat)
library(paropt)

context("test import")
path <- system.file("tests/testthat/files", package = "paropt")
df <- read.table(paste(path, "/test.txt", sep = ""), header = T)
df_lb <- read.table(paste(path, "/test_lb.txt", sep = ""), header = T)
df_ub <- read.table(paste(path, "/test_ub.txt", sep = ""), header = T)
final_list <- list(rep(150,4), rep(df[,1],4), c(df[,2], df[,3], df[,4], df[,5]),
                   c(df_lb[,2], df_lb[,3], df_lb[,4], df_lb[,5]),
                   c(df_ub[,2], df_ub[,3], df_ub[,4], df_ub[,5]))
names(final_list) <- c("cut", "time", "start", "lower", "upper")
df_NA <- read.table(paste(path, "/test_NA.txt", sep = ""), header = T)
df_lb_NA <- read.table(paste(path, "/test_lb_NA.txt", sep = ""), header = T)
df_ub_NA <- read.table(paste(path, "/test_ub_NA.txt", sep = ""), header = T)
final_list_with_NA <- list(c(149, rep(150,3)), c(df_NA[-3,1],rep(df_NA[,1],3)), c(df_NA[-3,2], df_NA[,3], df_NA[,4], df_NA[,5]),
                   c(df_lb_NA[-3,2], df_lb_NA[,3], df_lb_NA[,4], df_lb_NA[,5]),
                   c(df_ub_NA[-3,2], df_ub_NA[,3], df_ub_NA[,4], df_ub_NA[,5]))
names(final_list_with_NA) <- c("cut", "time", "start", "lower", "upper")

test_that("check Import_Parameter", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test2.txt", sep = ""),
                                     paste(path, "/test_lb.txt", sep = ""),
                                     paste(path, "/test_ub.txt", sep = "")), "ERROR: No File with start values")
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test.txt", sep = ""),
                                     paste(path, "/test_lb2.txt", sep = ""),
                                     paste(path, "/test_ub.txt", sep = "")), "ERROR: No File with lowerbound values")
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test.txt", sep = ""),
                                     paste(path, "/test_lb.txt", sep = ""),
                                     paste(path, "/test_ub2.txt", sep = "")), "ERROR: No File with upperbound values")
  expect_equal(paropt:::test_Import_Parameter(paste(path, "/test.txt", sep = ""),
                              paste(path, "/test_lb.txt", sep = ""),
                              paste(path, "/test_ub.txt", sep = "")), final_list)
  expect_equal(paropt:::test_Import_Parameter(paste(path, "/test_NA.txt", sep = ""),
                                     paste(path, "/test_lb_NA.txt", sep = ""),
                                     paste(path, "/test_ub_NA.txt", sep = "")), final_list_with_NA)
})

test_that("check Import_Parameter", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test_names.txt", sep = ""),
                                     paste(path, "/test_lb_names.txt", sep = ""),
                                     paste(path, "/test_ub_names.txt", sep = "")), "ERROR: Different of names for parameters between startvalues, lower bounds and upper bounds")
})

test_that("check Import_Parameter", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test_num_names.txt", sep = ""),
                                     paste(path, "/test_num_lb_names.txt", sep = ""),
                                     paste(path, "/test_num_ub_names.txt", sep = "")), "ERROR: Different number of names of parameters between startvalues, lower bounds and upper bounds")
})

test_that("check Import_Parameter", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test_nrow.txt", sep = ""),
                                     paste(path, "/test_lb_nrow.txt", sep = ""),
                                     paste(path, "/test_ub_nrow.txt", sep = "")), "ERROR: Different number of rows between startvalues, lower bounds and upper bounds")
})

test_that("check Import_Parameter", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test_time.txt", sep = ""),
                                     paste(path, "/test_lb_time.txt", sep = ""),
                                     paste(path, "/test_ub_time.txt", sep = "")), "ERROR:Different time vector for startvalues, lower bounds and upper bounds")
})

test_that("check Import_Parameter", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test_diffNA.txt", sep = ""),
                                     paste(path, "/test_lb_diffNA.txt", sep = ""),
                                     paste(path, "/test_ub_diffNA.txt", sep = "")), "ERROR: Different number of time_points between startvalues, lower bounds and upper bounds")
})


test_that("check if is numeric", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/test_non_numeric.txt", sep = ""),
                                               paste(path, "/test_non_numeric_lb.txt", sep = ""),
                                               paste(path, "/test_non_numeric_ub.txt", sep = "")), "\nERROR: Not a number in input. Only digits, NA, '.', '-', '+', 'e', 'E' are allowed.")
})


test_that("check if boundaries are violated", {
  expect_error(paropt:::test_Import_Parameter(paste(path, "/start_for_test_boundaries_violated.txt", sep = ""),
                                               paste(path, "/lb_over_ub.txt", sep = ""),
                                               paste(path, "/ub_lower_lb.txt", sep = "")), "ERROR: ParamClass_init boundary value error.")
})

list_states <- list(rep(150, 4), rep(df[,1], 4), c(df[,2], df[,3], df[,4], df[,5]))
names(list_states) <- c("cut", "time", "start")
list_states_NA <- list(rep(150, 4), rep(df[,1], 4), c(df_NA[,2], df_NA[,3], df_NA[,4], df_NA[,5]))
names(list_states_NA) <- c("cut", "time", "start")

test_that("check Import_States", {
  expect_equal(paropt:::test_Import_States(paste(path, "/test.txt", sep = "")), list_states)
  expect_equal(paropt:::test_Import_States(paste(path, "/test_NA.txt", sep = "")), list_states_NA)
})

lb <- data.frame(time = seq(0, 24, 2), a = 1:13, b = 1:13)
ub <- data.frame(time = seq(0, 24, 2), a = 100:112, b = 100:112)
res <- list()
res[[1]] <- c(13,13)
res[[2]] <- rep(seq(0,24,2),2)
res[[3]] <- rep(1:13, 2)
res[[4]] <- rep(100:112, 2)
names(res) <- c("cut", "time", "lower", "upper")
test_that("check import parameters", {
  expect_equal(paropt:::test_Import_Parameter_DF(lb, ub), res)
})

res2 <- list()
res2[[1]] <- c(13,13)
res2[[2]] <- rep(seq(0,24,2),2)
res2[[3]] <- rep(1:13, 2)
res2[[4]] <- c("time", "a", "b")
names(res2) <- c("cut", "time", "start", "header")
test_that("check import parameters", {
  expect_equal(paropt:::test_Import_Start_Parameter_DF(lb), res2)
})

test_that("check import states", {
  expect_equal(paropt:::test_Import_States_DF(lb), res2)
})

