#include "header.hpp"

typedef double (*fctptr)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc);


struct SwarmStruct {
  using namespace Rcpp;
  NumericMatrix S;
  NumericMatrix velocities;
  NumericVector best_errors;
  NumericVector current_errors;
  int best_particle;
  NumericVector parameter_of_best_particle;
  NumericMatrix personal_best_parameters;
  double global_best_error;
  int n_swarm;
  int n_params;
  fctptr lossfct;
  time_state_information_Rcpp_interface model;
  OS userode;
};

void calculate_errors_init(SwarmStruct& inp) {
  using namespace Rcpp;
  for(int i = 0; i < inp.n_swarm; i++) {
    inp.best_errors[i] = inp.lossfct(inp.S[i,_], inp.userode, inp.model);
    inp.current_errors[i] = inp.best_errors[i];
  }

}
