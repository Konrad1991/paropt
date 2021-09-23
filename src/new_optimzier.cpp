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

struct settingsPSO_Rcpp_interface {
  double err_tol;
  int pso_n_pop;
  int pso_n_gen;
  double pso_par_initial_w;
  double pso_par_w_max;
  double pso_par_w_min;
  double pso_par_w_damp;
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

  Rcpp::NumericVector::iterator it = std::min_element(sw.best_errors.begin(), sw.best_errors.end());
  int which_min = it - sw.best_errors.begin();
  sw.best_particle = which_min;
  sw.global_best_error = sw.best_errors(sw.best_particle);

}

struct neighbours {
  std::vector<int> neigh;
};

struct neighbour {
  int num_particle;
  std::vector<neighbours> N;
  std::vector<int> K;
};


void generate_random_int(Rcpp::IntegerVector& inp, std::vector<int>& res) {

  Rcpp::Vector<INTSXP> temp(1);
  for(int i = 0; i < res.size(); i++) {
    temp =  Rcpp::sample(inp, 1);
    res[i] = temp[0];
  }
}

int generate_random_int(Rcpp::IntegerVector& inp) {
  Rcpp::Vector<INTSXP> temp = Rcpp::sample(inp, 1);
   return temp[0];
}


void calc_neighbours(neighbour NE, int num_particle) {
  Rcpp::IntegerVector K = {0, 1, 2, 3};

  Rcpp::IntegerVector range(num_particle);
  for(int i = 0; i < range.size(); i++) {
    range[i] = i;
  }

  if(NE.N.size() != num_particle) {
      NE.num_particle = num_particle;
      NE.N.resize(NE.num_particle);
      NE.K.resize(NE.num_particle);
  }

  neighbours K0; K0.neigh = {0};
  neighbours K1; K1.neigh = {0, 1};
  neighbours K2; K2.neigh = {0, 1, 2};
  neighbours K3; K3.neigh = {0, 1, 2, 3};


  for(int i = 0; i < num_particle; i++) {
    NE.K[i] = generate_random_int(K);

    switch (NE.K[i]) {
      case 0:
        K0.neigh[0] = i;
        NE.N[i].neigh = K0.neigh;
        break;
      case 1:
        generate_random_int(range, K1.neigh);
        NE.N[i].neigh = K1.neigh;
        break;
      case 2:
        generate_random_int(range, K2.neigh);
        NE.N[i].neigh = K2.neigh;
        break;
      case 3:
        generate_random_int(range, K3.neigh);
        NE.N[i].neigh = K3.neigh;
        break;
    }
  }

}

void calculate_errors(SwarmStruct inp) {
  std::vector<double> temp(inp.n_params);
  for(int i = 0; i < inp.n_swarm; i++) {
    for(int j = 0; j < inp.n_params; j++) {
      temp[j] = inp.S(i, j);
    }
    inp.current_errors[i] = inp.lossfct(temp, inp.userode, inp.model);
  }
}

void check_boundaries(SwarmStruct inp, int n_params,
                      Rcpp::NumericVector lb, Rcpp::NumericVector ub, int index) {

    for(int i = 0; i < n_params; i++) {
      if(inp.S(index, i) > ub[i]) {
        inp.S(index, i) = ub[i];
      } else if(inp.S(index, i) < ub[i]) {
        inp.S(index, i) = lb[i];
      }
    }

}

Rcpp::List optimizer_pointer_new(
std::vector<double> t_lb, std::vector<double> t_ub,
settingsPSO_Rcpp_interface t_set_pso,
time_state_information_Rcpp_interface t_model,
double (*fctptr2)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc),
OS t_odes)
{

  int n_params = t_lb.size();
  int n_swarm = t_set_pso.pso_n_pop;
  int n_generations = t_set_pso.pso_n_gen;
  double desired_error = t_set_pso.err_tol;

  Rcpp::NumericVector lb(n_params);
  Rcpp::NumericVector ub(n_params);
  for(int i = 0; i < n_params; i++) {
    lb[i] = t_lb[i];
    ub[i] = t_ub[i];
  }
  fctptr solver = fctptr2;
  OS ode_system = t_odes;
  time_state_information_Rcpp_interface model = t_model;

  bool convergence = false;

  Rcpp::NumericVector local_best_errors(n_swarm);
  neighbour hood;
  std::vector<double> current_erros_of_hood;
  SwarmStruct SW;
  init_fct(SW, n_swarm, lb, ub, ode_system, solver);

  int position_local_best;
  int position_global_best;
  double global_best;
  Rcpp::NumericVector local_best_parameters(n_params);
  Rcpp::NumericVector rand1(n_params);
  Rcpp::NumericVector rand2(n_params);

  double cog;
  double soc;
  double initial_cog = 2.5;
  double initial_soc = 0.5;
  double final_cog = 0.5;
  double final_soc = 2.5;
  double par_w = 0.5;
  double par_w_max = 0.9;
  double par_w_min = 0.4;

  int i = 0;
  while(i < n_generations) {

    calc_neighbours(hood, n_swarm);

    par_w = par_w_max - i*(par_w_max - par_w_min)/n_generations;
    cog = initial_cog - (initial_cog - final_cog) * (i + 1) / n_generations;
    soc = initial_soc - (initial_soc - final_soc) * (i + 1) / n_generations;

    for(int j = 0; j < n_swarm; j++) {
      current_erros_of_hood.resize(hood.K[j]);
      for(int k = 0; k < hood.K[j]; k++) {
        current_erros_of_hood[k] = SW.best_errors[hood.N[j].neigh[k]];
      }
       position_local_best = std::min_element(current_erros_of_hood.begin(),
                                    current_erros_of_hood.end()) - current_erros_of_hood.begin();

      local_best_parameters = SW.S(position_local_best, Rcpp::_);
      rand1 = Rcpp::runif(n_params)*cog;
      rand2 = Rcpp::runif(n_params)*soc;

      SW.velocities(j, Rcpp::_) = par_w*SW.velocities(j, Rcpp::_) +
                    rand1*(SW.personal_best_parameters(j, Rcpp::_) - SW.S(j, Rcpp::_)) +
                    rand2*(SW.S(position_local_best, Rcpp::_) - SW.S(j, Rcpp::_));

      SW.S(j, Rcpp::_) = SW.S(j, Rcpp::_) + SW.velocities(j, Rcpp::_);

      check_boundaries(SW, n_params, lb, ub, j);

      calculate_errors(SW);

      if(SW.current_errors[j] < SW.best_errors[j]) {
        SW.best_errors[j] = SW.current_errors[j];
        SW.personal_best_parameters(j, Rcpp::_) = SW.S(j, Rcpp::_);
      }

    }

    position_global_best = std::min_element(SW.best_errors.begin(),
                                 SW.best_errors.end()) - SW.best_errors.begin();
    global_best = SW.best_errors(position_global_best);

    if(global_best < SW.global_best_error) {
      SW.global_best_error = global_best;
      SW.best_particle = position_global_best;
      SW.parameter_of_best_particle = SW.S(position_global_best, Rcpp::_);
      if(SW.global_best_error <= desired_error) {
        break;
      }
    }

    i++;
  }

  return Rcpp::List::create(Rcpp::Named("global best error") = SW.global_best_error,
                            Rcpp::Named("best parameters") = SW.parameter_of_best_particle);;
}
