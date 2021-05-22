## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- States, echo = F--------------------------------------------------------

df <- data.frame(time = seq(0,10,2), predator = c(10, 2.13, 0.23, 4.3, 0.37, 8.6), prey = c(10,0.0220,3.023, 0.028, 5, 0.33))
knitr::kable(
  head(df), caption = 'Table1: Possible input for states.'
)

## ---- Example for parameter input, echo = F-----------------------------------
set.seed(1)
df <- data.frame(time = seq(0,10,2), p1 = runif(6, 0.1, 3), p2 = runif(6,0,3), p3 = runif(6, 0.1, 1),
                 p4 = c(runif(1), rep(NA, 5)))
lb <- data.frame(time = seq(0,10,2), p1 = rep(0.1,6), p2 = rep(0, 6), p3 = rep(0.1, 6), p4 = c(0,rep(NA,5)))
ub <- data.frame(time = seq(0,10,2), p1 = rep(3,6), p2 = rep(3, 6), p3 = rep(1, 6), p4 = c(1,rep(NA,5)))
knitr::kable(list(df, lb, ub),
            caption =  'Table2: Possible input for parameters: start values, lower- and upper-bounds')

## ---- Parameters for predator-prey.model, echo = F----------------------------
set.seed(1)
df <- data.frame(time = 0, a = 1.0, b = 0.5, c = 0.2, d = 0.25)
lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)
knitr::kable(list(df, lb, ub),
            caption =  'Table3: parameter-input for predator-prey-model: start values, lower- and upper-bounds')

## ---- running optimization, eval = F------------------------------------------
#  path <- system.file("examples", package = "paropt")
#  library(paropt)
#  #Rcpp::sourceCpp(paste(path,"/ode.cpp", sep = ""))
#  #if you want compile ode-system on your system (already precompiled in package)
#  df <- read.table(paste(path,"/states_LV.txt", sep = ""), header = TRUE)
#  time <- df$time
#  param_start <- paste(path, "/start.txt", sep = "")
#  param_lb <- paste(path, "/lb.txt", sep = "")
#  param_ub <- paste(path, "/ub.txt", sep = "")
#  states <- paste(path, "/states_LV.txt", sep = "")
#  state_output <- paste(tempdir(), "/final_stateoutput.txt", sep = "")
#  par_output <- paste(tempdir(), "/optimized_params.txt", sep = "")
#  set.seed(1)
#  optimizer(time, paropt:::ode_example, 1e-6, c(1e-8, 1e-8),
#            param_start, param_lb, param_ub,
#            states, npop = 40, ngen = 200, error = 1,
#            state_output, par_output, "bdf")
#  df_in_silico <- read.table(paste(tempdir(), "/final_stateoutput.txt", sep = ""), header = TRUE)

## ---- results of optimization, eval = T, echo = F-----------------------------
path <- system.file("examples", package = "paropt")
df_in_silico <- read.table(paste(path, "/final_stateoutput.txt", sep = ""), header = T)
df <- read.table(paste(path,"/states_LV.txt", sep = ""), header = T)

plot(df$time, df$n1, pch = 19, main = "predator", ylab = "predator", xlab = "time")
points(df_in_silico$time, df_in_silico$n1, pch = 19, col = "darkred")
plot(df$time, df$n2, pch = 19, main = "prey",ylab = "prey", xlab = "time")
points(df_in_silico$time, df_in_silico$n2, pch = 19, col = "darkred")

## ---- optimizer, eval = F, echo = T-------------------------------------------
#  optimizer(time, ode, 1e-6, c(1e-8, 1e-8),
#            param_start, param_lb, param_ub,
#            states, npop = 40, ngen = 200, error = 1,
#            state_output, par_output, "bdf")

