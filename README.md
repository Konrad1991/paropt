<!-- badges: start -->
[![Travis build status](https://travis-ci.com/Konrad1991/paropt.svg?branch=Rcpp-Interface)](https://travis-ci.com/Konrad1991/paropt)
[![R-CMD-check](https://github.com/Konrad1991/paropt/workflows/R-CMD-check/badge.svg)](https://github.com/Konrad1991/paropt/actions)
https://img.shields.io/cran/l/paropt
https://img.shields.io/cran/l/paropt
<!-- badges: end -->

# paropt

Parameter Optimizing of ODEs

The package *paropt* is build in order to optimize parameters of ode-systems. Thus, the aim is that the output of the integration matches previously measured states. The user has to supply an ode-system in the form of a Rcpp-function. The information about states and parameters are passed via text-files. Additional information such as e.g. the relative tolerance are passed in R.

# Overview

The package *paropt* uses a modified particle swarm optimizer ('https://github.com/kthohr/optim') in order to find a global best solution. Furthermore, in order to evaluate each particle during optimzation four different solvers can be used all derived from SUNDIALS ('https://computing.llnl.gov/projects/sundials'). For more details see vignette. 

# Installation

*paropt* can be installed in R using install.packages("paropt") in order to use the version which is on CRAN. If you want to use the github-version (which can be used in parallel) use 
remotes::install_github("Konrad1991/paropt", ref = "Rcpp-Interface") within R (below you can see an example showing a parallized version of the parameter optimization of the predator-prey model). 

# Example

```Rcpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(paropt)]]
// [[Rcpp::plugins(cpp11)]]

typedef int (*OS)(double &t, std::vector<double> &params, std::vector<double> &states);

int ode_system(double &t, std::vector<double> &params, std::vector<double> & states) {
  
  // do not use any R-Code or R-Objects if the optimzation should run in parallel.
  // Users have to guarantee that the function can be called by several threads in parallel
  // define parameters (vector params contain the parameter in the order as defined in the corresponding textfiles)
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  
  // states have to be in the same order as specified in the textfile "states_LV.txt" (Otherwise the error-calculation does not work) 
  double n1 = states[0];
  double n2 = states[1];
  
  states[0] = n1*c*n2 - n1*d;
  states[1] = n2*a - n2*b*n1;
  
  return 0;
}


// [[Rcpp::export]]
Rcpp::XPtr<OS> test_optimization() {
  Rcpp::XPtr<OS> xpfun = Rcpp::XPtr<OS>(new OS(&ode_system));
  
  return xpfun;
}



/*** R
lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)

path <- system.file("examples", package = "paropt")
states <- read.table(paste(path,"/states_LV.txt", sep = ""), header = T)

library(paropt)
set.seed(1)
df <- optimizer_pointer(integration_times = states$time, ode_sys = test_optimization(),
                  relative_tolerance = 1e-6, absolute_tolerances = c(1e-8, 1e-8),
                  lower = lb, upper = ub, states = states, 
                  npop = 40, ngen = 30000, error = 0.0001, solvertype = "ADAMS")

par(mfrow = c(2,1))
plot(states$time, df$States[,1], pch = 19)
points(states$time, states$n1, pch = 19, col = "darkred")
plot(states$time, df$States[,2], pch = 19)
points(states$time, states$n2, pch = 19, col = "darkred")

*/
```
# Further plans

- Add different error calculations
- Give users the possibility to use their own spline (at least for the Rcpp-Interface)
- Update Documentation

# Acknowledgment

I would like to thank all the people who supported me writing this R-Package.
First of all I'm very greatful to my supervisor Arnd Heyer who supported me in writing the package.
Also I would like to thank Satyaprakash Nayak which gave me inspiration to solve coding problems.
Moreover, I would like to thank Tobias Bartsch for elaborately testing of the package.
Furthermore, gratitude is owned to Johannes Krämer for many discussion about simulations. 
Special thanks to my girlfriend Jana for supporting me all the way.

