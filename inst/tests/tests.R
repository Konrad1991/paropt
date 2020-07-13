library(testthat)
library(ParOpt2)

context("test basic functions")
path <- system.file("tests/testthat/files", package = "ParOpt2")

test_that("check dimensions of input file", {
  expect_equal(ParOpt2:::test_count_cols_rows(paste(path, "/test.txt", sep = "")), c(151,5))
})

test_that("check what happens when no file exist", {
  expect_error(ParOpt2:::test_no_file_exist(paste(path, "/test2.txt", sep = "")),
               "ERROR: No File")
})

test_that(("check count cols per row"), {
  expect_equal(ParOpt2:::test_check_ncols_per_row(paste(path, "/test.txt", sep = "")), rep(5, 151))
})

df <- read.table(paste(path, "/test.txt", sep = ""), header = T)
df_list <- list(df[,1], df[,2], df[,3], df[,4], df[,5])
test_that("check output of get_content", {
  expect_equal(ParOpt2:::test_get_content(paste(path, "/test.txt", sep = "")), df_list)
})

heads <- c("time", "Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
test_that("check output of get_header", {
  expect_equal(ParOpt2:::test_get_header(paste(path, "/test.txt", sep = "")), heads)
})
 
df_with_NA <- df
df_with_NA[1,1] <- NA
df_with_NA[10, 2] <- NA
df_list_without_NA <- list(df_with_NA[2:150,1],df_with_NA[-10,2],df_with_NA[,3],df_with_NA[,4],df_with_NA[,5])
test_that("check output of remove_NA", {
  expect_equal(ParOpt2:::test_remove_NA(paste(path, "/test.txt", sep = "")), df_list)
  expect_equal(ParOpt2:::test_remove_NA(paste(path, "/test_with_NA.txt", sep = "")), df_list_without_NA)
})



