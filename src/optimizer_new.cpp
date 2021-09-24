#include "optimizer_new.hpp"

typedef int (*OS)(double &t, std::vector<double> &params, std::vector<double> &states);
typedef double (*fctptr)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc);


void generate_random_int(Rcpp::IntegerVector& inp, std::vector<int>& res) {

  double temp;
  for(int i = 0; i < res.size(); i++) {
    GetRNGstate();
    temp =  R::runif(0, 1);
    PutRNGstate();
    res[i] = floor(temp*inp.size());
  }
}

int generate_random_int(Rcpp::IntegerVector& inp) {
  GetRNGstate();
  double temp = R::runif(0, 1);
  PutRNGstate();
   return floor(temp*inp.size());
}

void generate_random_double(Rcpp::NumericVector& inp) {
  for(int i = 0; i < inp.size(); i++) {
    GetRNGstate();
    inp[i] = R::runif(0, 1);
    PutRNGstate();
  }
}

struct SwarmStruct {
  Rcpp::NumericMatrix S; // shuffle!
  Rcpp::NumericMatrix velocities; // shuffle!
  Rcpp::NumericVector best_errors; // shuffle!
  Rcpp::NumericVector current_errors; // shuffle!
  int best_particle;
  Rcpp::NumericVector parameter_of_best_particle;
  Rcpp::NumericMatrix personal_best_parameters; // shuffle!
  double global_best_error;
  int n_swarm;
  int n_params;
  fctptr lossfct;
  time_state_information_Rcpp_interface model;
  OS userode;
};

void shuffle(SwarmStruct& inp) {
  int indx1;
  int indx2;

  double swap1;
  Rcpp::NumericVector swap2(inp.n_params);

  Rcpp::IntegerVector range(inp.n_swarm);
  for(int i = 0; i < range.size(); i++) {
    range[i] = i;
  }

  for(int i = 0; i < 1000; i++) { // 10000 better
    indx1 = arma::randi<int>(arma::distr_param(0, inp.n_swarm-1)); //generate_random_int(range);
    indx2 = arma::randi<int>(arma::distr_param(0, inp.n_swarm-1)); //generate_random_int(range);
    swap2 = inp.S(indx1, Rcpp::_);
    inp.S(indx1, Rcpp::_) = inp.S(indx2, Rcpp::_);
    inp.S(indx2, Rcpp::_) = swap2;
    swap2 = inp.velocities(indx1, Rcpp::_);
    inp.velocities(indx1, Rcpp::_) = inp.velocities(indx2, Rcpp::_);
    inp.velocities(indx2, Rcpp::_) = swap2;
    swap2 = inp.personal_best_parameters(indx1, Rcpp::_);
    inp.personal_best_parameters(indx1, Rcpp::_) = inp.personal_best_parameters(indx2, Rcpp::_);
    inp.personal_best_parameters(indx2, Rcpp::_) = swap2;

    swap1 = inp.best_errors[indx1];
    inp.best_errors[indx1] = inp.best_errors[indx2];
    inp.best_errors[indx2] = swap1;
    swap1 = inp.current_errors[indx1];
    inp.current_errors[indx1] = inp.current_errors[indx2];
    inp.current_errors[indx2] = swap1;
  }
}

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
    GetRNGstate();
    inp[i] = R::runif(0, 1);
    PutRNGstate();
  }

}


void init_fct(SwarmStruct& sw, int n_swarm,  Rcpp::NumericVector& lb, Rcpp::NumericVector& ub, OS ode_system, fctptr solver,
              time_state_information_Rcpp_interface& solv_param_struc) {

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
  sw.model = solv_param_struc;

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

void calc_neighbours(neighbour& NE, int num_particle) {
  Rcpp::IntegerVector K = {0, 1, 2, 3};

  Rcpp::IntegerVector range(num_particle);
  for(int i = 0; i < range.size(); i++) {
    range[i] = i;
  }

  NE.num_particle = num_particle;
  NE.N.resize(NE.num_particle);
  NE.K.resize(NE.num_particle);

  neighbours K0;
  K0.neigh.resize(1);
  neighbours K1;
  K1.neigh.resize(2);
  neighbours K2;
  K2.neigh.resize(3);
  neighbours K3;
  K3.neigh.resize(4);


  for(int i = 0; i < num_particle; i++) {
    NE.K[i] = arma::randi<int>(arma::distr_param(0,3)); //generate_random_int(K);

    if(NE.K[i] == 0) {
      K0.neigh[0] = i;
      NE.N[i].neigh.resize(1);
      NE.N[i].neigh = K0.neigh;
    } else if(NE.K[i] == 1) {
      generate_random_int(range, K1.neigh);
      K1.neigh = informants;
      NE.N[i].neigh.resize(2);
      NE.N[i].neigh = K1.neigh;
    } else if(NE.K[i] == 2) {
      generate_random_int(range, K2.neigh);
      NE.N[i].neigh.resize(3);
      NE.N[i].neigh = K2.neigh;
    } else if(NE.K[i] == 3) {
      generate_random_int(range, K3.neigh);
      NE.N[i].neigh.resize(4);
      NE.N[i].neigh = K3.neigh;
    }

  }

}


void check_boundaries(SwarmStruct& inp, int n_params,
                      Rcpp::NumericVector& lb, Rcpp::NumericVector& ub, int index) {

    for(int i = 0; i < n_params; i++) {
      if(inp.S(index, i) > ub[i]) {
        inp.S(index, i) = ub[i];
      } else if(inp.S(index, i) < lb[i]) {
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

  int convergence_check = 0;

  Rcpp::NumericVector local_best_errors(n_swarm);
  neighbour hood;
  std::vector<double> current_erros_of_hood;
  SwarmStruct SW;

  init_fct(SW, n_swarm, lb, ub, ode_system, solver, model);

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
  RcppThread::ThreadPool pool;
  std::vector<std::future<double> > futures(n_swarm);
  std::vector<double> temp(n_params);


  while(i < n_generations) {

    if(i == 0 || convergence_check >= 1) {
      shuffle(SW);
      calc_neighbours(hood, n_swarm);
      convergence_check = 0;
    }

    par_w = par_w_max - static_cast<double>(i)*(par_w_max - par_w_min)/static_cast<double>(n_generations);
    cog = initial_cog - (initial_cog - final_cog) * (i + 1) / n_generations;
    soc = initial_soc - (initial_soc - final_soc) * (i + 1) / n_generations;
    double v = 0.9;
    double chi;
    chi = 2*v;
    GetRNGstate();
    double omega = cog*R::runif(0, 1) + soc*R::runif(0,1);
    if(omega < 4.) {
      omega = 4.;
    }
    PutRNGstate();
    chi = chi/std::abs(2. - omega - std::sqrt(omega*(omega -4.) ));

    for(int j = 0; j < n_swarm; j++) {

      current_erros_of_hood.resize(hood.K[j]);

      for(int k = 0; k < current_erros_of_hood.size(); k++) {
        current_erros_of_hood[k] = SW.best_errors[hood.N[j].neigh[k]];
      }

      position_local_best = std::min_element(current_erros_of_hood.begin(),
                                    current_erros_of_hood.end()) - current_erros_of_hood.begin();

      local_best_parameters = SW.personal_best_parameters(hood.N[j].neigh[position_local_best], Rcpp::_);

      generate_random_double(rand1);
      generate_random_double(rand2);

      SW.velocities(j, Rcpp::_) = par_w*SW.velocities(j, Rcpp::_) +
                    cog*rand1*(SW.personal_best_parameters(j, Rcpp::_) - SW.S(j, Rcpp::_) ) +
                    soc*rand2*(local_best_parameters - SW.S(j, Rcpp::_) );

      SW.velocities(j, Rcpp::_) = SW.velocities(j, Rcpp::_)*chi;
      SW.S(j, Rcpp::_) = SW.S(j, Rcpp::_) + SW.velocities(j, Rcpp::_);

      check_boundaries(SW, n_params, lb, ub, j);
  }

  for(int j = 0; j < n_swarm; j++) {
      for(int k = 0; k < n_params; k++) {
        temp[k] = SW.S(j, k);
      }

      futures[j] = pool.pushReturn(SW.lossfct, std::ref(temp), SW.userode, std::ref(SW.model) );

      SW.current_errors[j] = futures[j].get();
  }

  for(int j = 0; j < n_swarm; j++) {
      if(SW.current_errors[j] < SW.best_errors[j]) {
        SW.best_errors[j] = SW.current_errors[j];
        SW.personal_best_parameters(j, Rcpp::_) = SW.S(j, Rcpp::_);
      }
  }


    position_global_best = std::min_element(SW.best_errors.begin(),
                                 SW.best_errors.end()) - SW.best_errors.begin();
    global_best = SW.best_errors[position_global_best];

    if(global_best < SW.global_best_error) {
      SW.global_best_error = global_best;
      SW.best_particle = position_global_best;
      SW.parameter_of_best_particle = SW.S(position_global_best, Rcpp::_);
      if(SW.global_best_error <= desired_error) {
        break;
      }
    } else {
      convergence_check++;
    }

    i++;

    if(i % 50 == 0) {
        Rcpp::Rcout << "generation is: " << i << "\t" << "global best error is:  " << SW.global_best_error << std::endl;
    }


    Rcpp::checkUserInterrupt();
  }

  Rcpp::Rcout << SW.parameter_of_best_particle << std::endl;

  return Rcpp::List::create(Rcpp::Named("global best error") = SW.global_best_error,
                            Rcpp::Named("best parameters") = SW.parameter_of_best_particle);;
}
