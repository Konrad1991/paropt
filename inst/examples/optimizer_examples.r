path <- system.file("examples", package = "paropt")
library(paropt)
#Rcpp::sourceCpp(paste(path,"/ode.cpp", sep = ""))
#if you want compile ode-system on your system (already precompiled in package)
df <- read.table(paste(path,"/states_LV.txt", sep = ""), header = TRUE)
time <- df$time
param_start <- paste(path, "/start.txt", sep = "")
param_lb <- paste(path, "/lb.txt", sep = "")
param_ub <- paste(path, "/ub.txt", sep = "")
states <- paste(path, "/states_LV.txt", sep = "")
state_output <- paste(tempdir(), "/final_stateoutput.txt", sep = "")
par_output <- paste(tempdir(), "/optimized_params.txt", sep = "")
set.seed(1)
optimizer(time, paropt:::ode_example, 1e-6, c(1e-8, 1e-8),
          param_start, param_lb, param_ub,
          states, npop = 10, ngen = 200, error = 3,
          state_output, par_output, "bdf")
df_in_silico <- read.table(paste(tempdir(), "/final_stateoutput.txt", sep = ""), header = TRUE)
plot(df$time, df$n1, pch = 19, main = "predator", ylab = "predator", xlab = "time")
points(df_in_silico$time, df_in_silico$n1, pch = 19, col = "darkred")
legend(1, 26, legend = c("measured", "in silico"),
       col = c("black", "darkred"),lty = 1:2, cex = 0.8)
plot(df$time, df$n2, pch = 19, main = "prey",ylab = "prey", xlab = "time")
points(df_in_silico$time, df_in_silico$n2, pch = 19, col = "darkred")
legend(1, 26, legend = c("measured", "in silico"),
       col = c("black", "darkred"),lty = 1:2, cex = 0.8)
