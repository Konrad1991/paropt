#ifndef COMPARISONFCT
#define COMPARISONFCT

#include "header.hpp"

double comparison_fct(std::vector<double> &parameter, SEXP odesystem, time_state_information_two_stage &info_additional);
void wrapper_ode_system_two_stage(SEXP ode_system, std::vector<double> &parameter,
                                 std::vector<double> &parameter_time, std::vector<double> &parameter_cut_idx,
                                 std::vector<double> &states_measured_at_time, std::vector<double> &new_states_deriv_in_silico, double &time);

#endif // COMPARISONFCT
