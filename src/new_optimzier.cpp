//#include "header.hpp"

#include <Rcpp.h>

struct time_state_information_Rcpp_interface {
  std::vector<double> init_state;
  std::vector<double> par_times;
  std::vector<int> param_idx_cuts;
  std::vector<double> state_measured;
  std::vector<double> state_times;
  std::vector<int> state_idx_cut;
  std::vector<double> integration_times;
  double reltol;
  std::vector<double> absolute_tolerances;
};

typedef int (*OS)(double &t, std::vector<double> &params, std::vector<double> &states);
typedef double (*fctptr)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc);


struct SwarmStruct {
  Rcpp::NumericMatrix S;
  Rcpp::NumericMatrix velocities;
  Rcpp::NumericVector best_errors;
  Rcpp::NumericVector current_errors;
  int best_particle;
  Rcpp::NumericVector parameter_of_best_particle;
  Rcpp::NumericMatrix personal_best_parameters;
  double global_best_error;
  int n_swarm;
  int n_params;
  fctptr lossfct;
  time_state_information_Rcpp_interface model;
  OS userode;
};

void calculate_errors_init(SwarmStruct& inp) {
  using namespace Rcpp;
  std::vector<double> temp(inp.n_params);
  for(int i = 0; i < inp.n_swarm; i++) {
    for(int j = 0; j < inp.n_params; j++) {
      temp[j] = inp.S(i, j);
    }
    inp.best_errors[i] = inp.lossfct(temp, inp.userode, inp.model);
    inp.current_errors[i] = inp.best_errors[i];
  }

}


void random_num_vec(Rcpp::NumericVector& inp) {
  for(int i = 0; i < inp.size(); i++) {
    inp[i] = R::runif(0, 1);
  }

}


void init_fct(SwarmStruct sw, int n_swarm,  Rcpp::NumericVector lb, Rcpp::NumericVector ub, OS ode_system, fctptr solver) {

  Rcpp::RNGScope scope;
  sw.n_swarm = n_swarm;
  sw.n_params = lb.size();
  sw.S = Rcpp::NumericMatrix(n_swarm, lb.size());
  sw.velocities = Rcpp::NumericMatrix(n_swarm, lb.size());
  sw.personal_best_parameters = Rcpp::NumericMatrix(n_swarm, lb.size());
  sw.best_errors  = Rcpp::NumericVector(n_swarm);
  sw.current_errors = Rcpp::NumericVector(n_swarm);
  sw.parameter_of_best_particle = Rcpp::NumericVector(lb.size());
  sw.userode = ode_system;
  sw.lossfct = solver;

  Rcpp::NumericVector temp(lb.size());

  for(int i = 0; i < n_swarm; i++) {
    random_num_vec(temp);
    sw.S(i, Rcpp::_) = lb + (ub - lb)*temp;
    sw.personal_best_parameters(i, Rcpp::_) = sw.S(i, Rcpp::_);
  }

  calculate_errors_init(sw);

  NumericVector::iterator it = std::min_element(sw.best_errors.begin(), sw.best_errors.end());
  int which_min = it - sw.best_errors.begin();
  sw.best_particle = which_min;
  sw.global_best_error = sw.best_errors(sw.best_particle);

}
