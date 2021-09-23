#include "header.hpp"
#include "paropt_types.h"
#include <thread>
#include <mutex>


Rcpp::List optimizer_pointer_new(
std::vector<double> t_lb, std::vector<double> t_ub,
settingsPSO_Rcpp_interface t_set_pso,
time_state_information_Rcpp_interface t_model,
double (*fctptr2)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc),
OS t_odes);
