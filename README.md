<!-- badges: start -->
[![Travis build status](https://travis-ci.com/Konrad1991/paropt.svg?branch=Rcpp-Interface)](https://travis-ci.com/Konrad1991/paropt)
![License](https://img.shields.io/cran/l/paropt)
![Version](https://img.shields.io/cran/v/paropt)
[![](http://cranlogs.r-pkg.org/badges/last-month/paropt?color=green)](https://cran.r-project.org/package=paropt)
<!-- badges: end -->





# paropt

Parameter Optimizing of ODEs

The package *paropt* is build in order to optimize parameters of ode-systems. Thus, the aim is that the output of the integration matches previously measured states. The user has to supply an ode-system in the form of an R or Rcpp function. In case it is an R function the code is translated into C++ code (example can be seen below). 

The information about states and parameters are passed as data.frames. Additional information such as e.g. the relative tolerance are passed in R.

# Overview

The package *paropt* uses a modified particle swarm optimizer ('https://github.com/kthohr/optim') in order to find a global best solution. Furthermore, in order to evaluate each particle during optimzation four different solvers can be used all derived from SUNDIALS ('https://computing.llnl.gov/projects/sundials'). For more details see vignette. 

# Installation

*paropt* can be installed in R using install.packages("paropt") in order to use the version which is on CRAN. If you want to use the github-version (which can be used in parallel) use
remotes::install_github("Konrad1991/paropt") within R (below you can see an example showing a parallized version of the parameter optimization of the predator-prey model).

# Example

```R
remotes::install_github("Konrad1991/paropt", force = TRUE)

library(paropt)

ode <- function(t, parameter, y, ydot) {
  
  a_db = at(parameter, 1)
  b_db = at(parameter, 2)
  c_db = at(parameter, 3)
  d_db = at(parameter, 4)
  
  predator_db = at(y,1)
  prey_db = at(y, 2)
  
  ydot[1] = predator_db*prey_db*c_db - predator_db*d_db
  ydot[2] = prey_db*a_db - prey_db*predator_db*b_db
}

# compile
r <- paropt::convert(ode, verbose = TRUE)

path <- system.file("examples", package = "paropt")
states <- read.table(paste(path,"/states_LV.txt", sep = ""), header = TRUE)

# parameter
lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)

# Optimizing
set.seed(1)

start_time <- Sys.time()
df <- paropt::po(integration_times = states$time, ode_sys = r(),
                     relative_tolerance = 1e-6, absolute_tolerances = c(1e-8, 1e-8),
                     lower = lb, upper = ub, states = states, 
                     npop = 40, ngen = 1000, error = 0.0001, solvertype = "bdf")
end_time <- Sys.time()
end_time - start_time


start <- data.frame(df[[8]])
names(start) <- names(lb)
df2 <- paropt::so(integration_times = states$time, fctptr = r(),
                              relative_tolerance = 1e-6, absolute_tolerances = c(1e-8, 1e-8),
                              start = start, states = states, solvertype = "bdf")

par(mfrow = c(2,1))
plot(states$time, df$States[,1], pch = 19, type = 'l', ylab = "predator", xlab = "time", ylim = c(0, 30))
points(states$time, states$n1, pch = 19, col = "darkred", type = 'p')
points(states$time, df2$`in silico states`[,1], pch = 12, col = "darkgreen", type = 'p')
legend(80, 30, legend=c("in silico", "measured"),
       col=c("black", "darkred"), lty=1, cex=0.8)
plot(states$time, df$States[,2], pch = 19, type = 'l', ylab = "prey", xlab = "time", ylim = c(0, 65))
points(states$time, states$n2, pch = 19, col = "darkred", type = 'p')
points(states$time, df2$`in silico states`[,2], pch = 12, col = "darkgreen", type = 'p')
legend(80, 60, legend=c("in silico", "measured"),
       col=c("black", "darkred"), lty=1, cex=0.8)
```
# Further plans

- Add different error calculations
- Give users the possibility to use different interpolation functions. 

# Acknowledgment

I would like to thank all the people who supported me writing this R-Package.
First of all I'm very greatful to my supervisor Arnd Heyer who supported me in writing the package.
Also I would like to thank Satyaprakash Nayak which gave me inspiration to solve coding problems.
Moreover, I would like to thank Tobias Bartsch for elaborately testing of the package.
Furthermore, gratitude is owned to Johannes Krämer for many discussion about simulations.
Special thanks to my girlfriend Jana for supporting me all the way.
