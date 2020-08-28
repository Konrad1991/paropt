# paropt

Parameter Optimizing of ODEs

The package *paropt* is build in order to optimize parameters of ode-systems. Thus, the aim is that the output of the integration matches previously measured states. The user has to supply an ode-system in the form of a Rcpp-function. The information about states and parameters are passed via text-files. Additional information such as e.g. the relative tolerance are passed in R.

# Overview

The package *paropt* uses a particle swarm optimizer ('https://github.com/kthohr/optim') in order to find a global best solution. Furthermore, in order to evaluate each particle during optimzation four different solvers can be used all derived from SUNDIALS ('https://computing.llnl.gov/projects/sundials'). For more details see vignette. 

# Example

```Rcpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(paropt)]]
// [[Rcpp::plugins(cpp11)]]

typedef int (*OS)(double &t, std::vector<double> &params, std::vector<double> &states);

int ode_system(double &t, std::vector<double> &params, std::vector<double> & states) {
  // define parameters
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  
  double n1 = states[0];
  double n2 = states[1];
  
  states[0] = n1*c*n2 - n1*d;
  states[1] = n2*a - n2*b*n1;
  
  return 0;
}

// [[Rcpp::export]]
Rcpp::List test_integration(std::vector<double> integration_times,
                            std::string start,
                            std::string lb,
                            std::string ub,
                            std::string states,
                            std::string output,
                            std::string output_parameter,
                            std::string solvertyp) {
  std::vector<double> abs_tol = {1e-8, 1e-8};
  
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("paropt");
  Rcpp::Function f = pkg["optimizer_access_in_Rcpp"];
  
  Rcpp::XPtr<OS> xpfun = Rcpp::XPtr<OS>(new OS(&ode_system));
  
  return Rcpp::List::create(f(integration_times,
                              xpfun,
                              1e-6,
                              abs_tol,
                              start,
                              lb,
                              ub,
                              states,
                              40,
                              1500,
                              0.001,
                              output,
                              output_parameter,
                              solvertyp) );
}



/*** R
# ===========================================================
# Simple example optimizing parameter of predator-prey model:


# Get state values of predator and prey and save it in working directory as textfile:
path <- system.file("examples", package = "paropt")
df <- read.table(paste(path,"/states_LV.txt", sep = ""), header = T)

setwd("~/Parallel")
write.table(df, "states_LV.txt", quote = F, row.names = F)

# Define startvalues, lower-bounds and upper-bounds for parameters:
st <- data.frame(time = 0., a = 1.1, b = 0.4, c = 0.1, d = 0.4)
write.table(st, "start.txt", quote = F, row.names = F)

lb <- data.frame(time = 0., a = 0.8, b = 0.3, c = 0.09, d = 0.09)
write.table(lb, "lb.txt", quote = F, row.names = F)

ub <- data.frame(time = 0., a = 1.3, b = 0.7, c = 0.4, d = 0.7)
write.table(ub, "ub.txt", quote = F, row.names = F)

# Define variables containing names of textfiles:
start <-  "start.txt"
lb <-  "lb.txt"
ub <-  "ub.txt"
states <- "states_LV.txt"
output <- "output.txt"
output_par <- "output_par.txt"
solvertyp <- "bdf"

# Run simulation (parallel):
set.seed(1)
df <- read.table("states_LV.txt", header = T)
test_integration(df$time, start, lb, ub, states, output, output_par, solvertyp)

# Plot results of simulation:
df <- read.table("states_LV.txt", header = T)
is <- read.table("output.txt", header = T)

plot(df$time, df$n1)
points(is$time, is$n1, pch = 19, col = "darkred")
plot(df$time, df$n2)
points(is$time, is$n2, pch = 19, col = "darkred")
*/
```
# Further plans

- Add different error calculation
- Give users the possibility to use their own spline (at least in the Rcpp-Interface)

# Acknowledgment

I would like to thank all the people who supported me writing this R-Package.
First of all I'm very greatful to my supervisor Arnd Heyer who supported me in writing the package.
Also I would like to thank Satyaprakash Nayak which gave me inspiration to solve coding problems.
Moreover, I would like to thank Tobias Bartsch for elaborately testing of the package.
Furthermore, gratitude is owned to Johannes Krämer for many discussion about simulations. 
Special thanks to my girlfriend Jana for supporting me all the way.

