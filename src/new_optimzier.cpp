#include "header.hpp"



typedef VD std::vector<double>;
typedef MD std::vector<std::vector< double> >;
typedef VI std::vector<int>;
typedef MI std::vector<std::vector< int> >;
typedef double (*fctptr)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc);

void pop_loop() {

}

void calc_neighberhood() {
  int idx1, idx2;
  double best_idx1, best_idx2;
  double values_idx1, values_idx2;

  VD swap(NUM_PAR);

  for(int i = 0; i < 1000; i++) { // z = number of switches
    GetRNGstate();
    idx1 = round(R::runif(0, SIZE_SWARM - 1));
    PutRNGstate();
    GetRNGstate();
    idx2 = round(R::runif(0, SIZE_SWARM - 1));
    PutRNGstate();

    for(int j = 0; j < NUM_PAR; j++) {
      swap[j] = SWARM[idx1][j];
      SWARM[idx1][j] = SWARM[idx2][j];
      SWARM[idx2][j] = swap[j];

      swap[j] = VELOCITIES[idx1][j];
      VELOCITIES[idx1][j] = VELOCITIES[idx2][j];
      VELOCITIES[idx2][j] = swap[j];
    }


    best_vals_idx1 = best_vals(idx1);
    best_vals_idx2 = best_vals(idx2);
    best_vals(idx1) = best_vals_idx2;
    best_vals(idx2) = best_vals_idx1;
    objfn_vals_idx1 = objfn_vals(idx1);
    objfn_vals_idx2 = objfn_vals(idx2);
    objfn_vals(idx1) = objfn_vals_idx2;
    objfn_vals(idx2) = objfn_vals_idx1;
  }
}




  // =============================
  // Calculate neighberhood
  // =============================
  int n_informants;
  int K = 3; // Neighberhoodsize K = 3
  for(int i = 0; i < n_pop; i++) { // 0. Loop ueber Reihen
    GetRNGstate();
      n_informants = arma::randi<int>(arma::distr_param(0,K));
    PutRNGstate();
    GetRNGstate();
      arma::ivec informants = arma::randi(n_informants, arma::distr_param(0, n_pop-1));
    PutRNGstate();
    GetRNGstate();
      arma::ivec temp_row = arma::zeros<arma::ivec>(n_pop);
    PutRNGstate();
      temp_row(i) = 1;
      for(int m = 0; m < n_informants; m++) {
        temp_row(informants(m)) = 1;
      }
      arma::irowvec test = temp_row.t();
      nbhood_precursor.row(i) = temp_row.t();
  }
      std::vector<int> size_counter(n_pop);
      for(int i = 0; i < n_pop; i++) {
          int temp_counter = 0;
          for(int j = 0; j < n_pop; j++) {
            if(nbhood_precursor(j,i) == 1) {
              temp_counter = temp_counter + 1;
            }
          }
          size_counter[i] = temp_counter;
          neighberhood[i].resize(size_counter[i]);
          int m_counter = 0;
          for(int m = 0; m < n_pop; m++) {
            if(nbhood_precursor(m,i) == 1) {
              neighberhood[i][m_counter] = m;
              m_counter = m_counter + 1;
            }
          }
      }
  convergence_check = 0;
}


result_pso particle_swarm(VD lb, VD ub, settingsPSO_Rcpp_interface set_pso, time_state_information_Rcpp_interface info, fctptr FP, OS ode) {

  /*
  variables needed for pso
  */
  int NUM_PAR = lb.size();

  int SIZE_SWARM = set_pso.pso_n_pop;
  int NUM_GEN = set_pso.pso_n_gen;
  double TOL = set_pso.err_tol;
  double W = set_pso.pso_par_initial_w;
  double WMAX = set_pso.pso_par_w_max;
  double WMIN = set_pso.pso_par_w_min;
  double DAMP = set_pso.pso_par_w_damp;

  double COG = 2.05;
  double SOC = 2.05;
  double COG_INIT = 2.5;
  double COG_FINAL = 0.5;
  double SOC_INIT = 0.5;
  double SOC_FINAL = 2.5;

  double current_error;
  VD PARTICLE_VALUES(SIZE_SWARM);
  MD SWARM(SIZE_SWARM);
  MD VELOCITIES(SIZE_SWARM);
  for(int i = 0; i < SIZE_SWARM; i++) {
    SWARM.resize(NUM_PAR);
    VELOCITIES.resize(NUM_PAR);
    for(int j = 0; j < NUM_PAR; j++) {
      SWARM[i][j] = 0.;
      VELOCITIES[i][j] = 0.;
    }
  }
  VD PARAM_TEMP(NUM_PAR);
  MI nbhood_precursor(SIZE_SWARM);
  for(int i = 0; i < SIZE_SWARM; i++) {
    nbhood_precursor.resize(SIZE_SWARM);
    for(int j = 0; j < SIZE_SWARM; j++) {
      if(i == j) {
        nbhood_precursor[i][j] = 1;
      }
    }
  }
  MI neighberhood(SIZE_SWARM);

  RcppThread::ThreadPool pool;
  std::vector<std::future<double> > futures(SIZE_SWARM);


  // initialize
  // =============================
  for(int i = 0; i < SIZE_SWARM; i++) {
    for(int j = 0; j < NUM_PAR; j++) {
      GetRNGstate();
      SWARM[i][j] = lb[i] + (ub[i] - lb[i])*R::runif(0, 1);
      PutRNGstate();

      if(SWARM[i][j] > ub[j]) {
        SWARM[i][j] = ub[j];
      } else if(SWARM[i][j] < lb[j]) {
        SWARM[i][j] = lb[j];
      }
      PARAM_TEMP[j] = SWARM[i][j];
    }
    current_error = FP(PARAM_TEMP, ode, info);
    PARTICLE_VALUES[i] = current_error;
  }

  VD GLOBAL_BEST_ERRORS(SIZE_SWARM);
  for(int i = 0; i < SIZE_SWARM; i++) {
    GLOBAL_BEST_ERRORS[i] = PARTICLE_VALUES[i];
  }

  MD GLOBAL_BEST_PARAMETERS(SIZE_SWARM);
  for(int i = 0; i < SIZE_SWARM; i++) {
    GLOBAL_BEST_PARAMETERS.resize(NUM_PAR);
    for(int j = 0; j < NUM_PAR; j++) {
      GLOBAL_BEST_PARAMETERS[i][j] = SWARM[i][j];
    }
  }

  int index_global_best = std::min_element(GLOBAL_BEST_ERRORS.begin(), GLOBAL_BEST_ERRORS.end());
  double BEST_ERROR = GLOBAL_BEST_ERRORS[index_global_best];
  VD BEST_PARAMETERS(NUM_PAR);

  for(int i = 0; i < NUM_PAR; i++) {
    BEST_PARAMETERS[i] = SWARM[index_global_best][i];
  }

  double err = err_tol + 1.0;
  int convergence_check = 0;
  int convergence_check_stop = 0;



}
