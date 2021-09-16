#include "header.hpp"



typedef VD std::vector<double>;
typedef MD std::vector<std::vector< double> >;
typedef VI std::vector<int>;
typedef MI std::vector<std::vector< int> >;
typedef double (*fctptr)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc);

/*
calculate one generation
*/
void pop_loop() {

}

/*
shuffles rows of SWARM, VELOCITIES, GLOBAL_BEST_ERRORS, GLOBAL_BEST_PARAMETERS, BEST_PARAMETERS
*/
void shuffling(MD& sw, MD& velo, VD& errors, VD& params) {

  Rcpp::IntegerVector pop_num( (sw.size() - 1) );
  for(int i = 0; i < pop_num.size(); i++) {
    pop_num[i] = i;
  }

  int idx1;
  int idx2;
  double swap;

  for(size_t z = 0; z < 1000; z++) { // z = number of switches
    GetRNGstate();
    idx1 = Rcpp::sample(pop_num, 1);
    PutRNGstate();
    GetRNGstate();
    idx2 = Rcpp::sample(pop_num, 1);
    PutRNGstate();

    for(size_t i = 0; i < sw[0].size(); i++) {
      swap = sw[idx1][i];
      sw[idx1][i] = sw[idx2][i];
      sw[idx2][i] = swap;

      swap = velo[idx1][i];
      velo[idx1][i] = velo[idx2][i];
      velo[idx2][i] = swap;


    }

  }
}


/*
calculate neighberhood
*/
void calc_neighberhood() {

  // 1. shuffle

  // 2. fill nbhood_precursor

  // 3. fill neighberhood
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


  for(int i = 0; i < NUM_GEN; i++) {

    // 1. calc_neighberhood


    // 2. update parameter

    // 3. call pop_loop

    // 4. update and check if current_err < err....


  }

}
