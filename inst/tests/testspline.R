library(testthat)
library(ParOpt2)

context("test spline")
path <- system.file("tests/testthat/files", package = "ParOpt2")
 
df <- read.table(paste(path, "/par_spline.txt", sep = ""), header = T)
out <- ParOpt2:::test_paramsort_and_spline(seq(0, 10, 1), paste(path, "/par_spline.txt", sep = ""),
                          paste(path, "/lb_spline.txt", sep = ""), paste(path, "/ub_spline.txt", sep = ""))
out <- do.call(rbind, out)
df2 <- cbind(df[,2], df[,3])
test_that("check ode solving", {
  expect_equal(out, df2)
})


out <- ParOpt2:::test_paramsort_and_spline(seq(0, 10, 0.25), paste(path, "/par_spline.txt", sep = ""),
                                 paste(path, "/lb_spline.txt", sep = ""), paste(path, "/ub_spline.txt", sep = ""))
out <- do.call(rbind, out)
df3a <- spline(seq(0,10,1), df$a, xout = seq(0,10, 0.25))
df3b <- spline(seq(0,10,1), df$b, xout = seq(0,10, 0.25))
df3 <- cbind(df3a$y, df3b$y)
test_that("check ode solving", {
  expect_equal(out, df3, tolerance = 0.1)
})


