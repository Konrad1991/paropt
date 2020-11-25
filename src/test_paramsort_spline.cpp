/* !!revision!!
remove compiler warnings
*/

#include "solver.hpp"

struct test_parameter{
  std::vector<double> parameter;
  std::vector<double> parameter_time;
  std::vector<int> parameter_cut_idx;
};

int test_call_paramsort(realtype t, void *user_data, std::vector<double> &splined_parameters) {

  // cast pointer to structure and store elements
  struct test_parameter *my_ode_system = (struct test_parameter*)user_data;
  std::vector<double> params = (*my_ode_system).parameter;
  std::vector<double> params_time = (*my_ode_system).parameter_time;
  std::vector<int> params_cut_idx = (*my_ode_system).parameter_cut_idx;

  //Rcpp::Rcerr << params_cut_idx.size() << std::endl;
  //Rcpp::Rcerr << std::endl;
  //for(size_t i = 0; i < params_cut_idx.size(); i++) {
    //Rcpp::Rcerr << params_cut_idx[i] << std::endl;
  //}
  // interpolate Parameter if necessary
  std::vector<double> parameter_input;
  params_sort(t, parameter_input, params, params_time, params_cut_idx);

  splined_parameters.resize(parameter_input.size());
  for(size_t i = 0; i < splined_parameters.size(); i++) {
    //Rcpp::Rcerr << parameter_input[i] << std::endl;
    splined_parameters[i] = parameter_input[i];
  }

  return 0;
}

/*
// [[Rcpp::export]]
std::vector<std::vector<double> > test_paramsort_and_spline(Rcpp::NumericVector timepoints, std::string start,std::string lower, std::string upper) {

std::vector<int> params_cut_idx_vec;
std::vector<double> params_time_combi_vec;
std::vector<double> param_combi_start;
std::vector<double> param_combi_lb;
std::vector<double> param_combi_ub;
std::vector<std::string> header_parameter;

Import_Parameter(start, lower, upper, params_cut_idx_vec, params_time_combi_vec,
param_combi_start, param_combi_lb, param_combi_ub, header_parameter);

struct test_parameter test = {param_combi_start, params_time_combi_vec, params_cut_idx_vec};
void* ptr_to_test = &test;

//Rcpp::Rcerr << params_cut_idx_vec.size() << std::endl;

std::vector<std::vector<double> > output(timepoints.length());
for(size_t i = 0; i < output.size(); i++) {
  output[i].resize(params_cut_idx_vec.size());
}

for(size_t i = 0; i < timepoints.length(); i++) {
    realtype temp = timepoints[i];
    test_call_paramsort(temp, ptr_to_test, output[i]);
}

std::vector<std::vector<double> > testimport(3);
testimport[0].resize(param_combi_start.size());
testimport[1].resize(params_time_combi_vec.size());
testimport[2].resize(params_cut_idx_vec.size());

for(size_t i = 0; i < testimport.size(); i++) {
  for(size_t j = 0; j < testimport[i].size(); j++) {
    if(i == 0) {
    testimport[i][j] = param_combi_start[j]; } else if(i == 1) {
      testimport[i][j] = params_time_combi_vec[j];
    } else if(i == 2) {
      testimport[i][j] = params_cut_idx_vec[j];
    } else {Rcpp::Rcerr << "Problem" << std::endl;}
  }
}

return output; //testimport; //output;
}
*/
