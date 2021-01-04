
/*
BSD 3-Clause License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

The names of its contributors may not be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef SOLVER
#define SOLVER

#include "header.hpp"

//static int check_retval(void *returnvalue, const char *funcname, int opt);
double solver_bdf(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information &param_model);

double solver_bdf_save(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, Rcpp::NumericMatrix &DF);

double solver_adams(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information &param_model);

double solver_adams_save(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, Rcpp::NumericMatrix &DF);

double solver_erk(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information &solv_param_struc);

double solver_erk_save(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, Rcpp::NumericMatrix &DF);

double solver_ark(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information &solv_param_struc);

double solver_ark_save(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, Rcpp::NumericMatrix &DF);

int wrapper_ode_system(realtype t, N_Vector y, N_Vector ydot, void *user_data);

bool double_diff(double x, double y);

#endif // SOLVER
