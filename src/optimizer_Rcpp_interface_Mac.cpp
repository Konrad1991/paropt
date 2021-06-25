/* !!revision!!
completly set up a new PSO.

1. classical pso
2. pso as clerc propose it
3. Initialisation in own fct
4. calculation of neighberhood in own fct

5. implement TPO

6. preparation of multithreading have to be improved.
   Maybe it is sometimes better not to use all threads. Due to large overhead.
   Especially a problem if solving is fast.
   Thus, it is necessary to test how fast solving is and based on this choose number of threads.
   Give user possibility to define number of threads which should be used
Maybe it is possible to have only one PSO; not well arranged though


*/



/*################################################################################
  ##
  ##   Copyright (C) 2016-2018 Keith O'Hara
  ##
  ##   This file is part of the OptimLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * Particle Swarm Optimization (PSO)
*/

/*
OptimLib: A Lightweight C++ Library of Numerical Optimization Methods for Nonlinear Functions
Copyright 2016-2018 Keith O'Hara

This product includes software developed by Keith O'Hara (http://www.kthohr.com)
*/

/*
// ===============================================================
The particle swarm algorithm in this file is a modified version of the code from Keith O'Hara.
For details see: https://github.com/kthohr/optim
// ===============================================================
*/

#include "optimizer_Rcpp_interface_Mac.hpp"


typedef double (*fctptr2)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc);
//std::mutex mtx;

void helper_fct (arma::mat inp, arma::vec & errors, int index, int number_of_rows, const int nvals, OS odes, time_state_information_Rcpp_interface model,
                 fctptr2 fctptr){
/*
  assert(nvals > 0);
  std::vector<double> param_temp(nvals);
  double prop_objfn_val = 0.;
  for(int j = 0; j < number_of_rows; j++) {
  mtx.lock();
    for(int i = 0; i < nvals; i++) {
      param_temp[i] = inp(j, i);
    }

  OS ODE = odes;
  time_state_information_Rcpp_interface MODEL = model;
  mtx.unlock();

  prop_objfn_val = fctptr(param_temp, ODE, MODEL);

  mtx.lock();
  errors(index + j) = prop_objfn_val;
  mtx.unlock();
  }
  */
}


Optimizer_Rcpp_interface::Optimizer_Rcpp_interface(
        std::vector<double> t_lb,
        std::vector<double> t_ub,
        settingsPSO_Rcpp_interface t_set_pso,
        time_state_information_Rcpp_interface t_model,
        double (*fctptr2)(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc) ,
        OS t_odes ) { //Rcpp::XPtr<OS> t_odes

    m_lb.resize(t_lb.size());
    m_ub.resize(t_ub.size());
    for(size_t i = 0; i < t_lb.size(); i++) {
        m_lb[i] = t_lb[i];
        m_ub[i] = t_ub[i];
    }
    setpso = t_set_pso;
    model = t_model;
    fctptr = fctptr2;
    odes = t_odes;
}

// ========================================================================
// ========================================================================
// ========================================================================

void Optimizer_Rcpp_interface::get_best_particle_param_values(std::vector<double> & temp) {
    temp.resize(best_particle_param_values.size());
    for(size_t   i = 0; i < temp.size(); i++) {
        temp[i] = best_particle_param_values[i];
    }
}

// ========================================================================
// ========================================================================
// ========================================================================

double Optimizer_Rcpp_interface::pso() { // (labled with ! need check)

  // define parameters
  // ============================
  const int n_vals = m_lb.size();
  const double err_tol = setpso.err_tol;
  const int n_pop = setpso.pso_n_pop;
  const int n_gen = setpso.pso_n_gen;
  double par_w = setpso.pso_par_initial_w;
  const double par_w_max = setpso.pso_par_w_max;
  const double par_w_min = setpso.pso_par_w_min;
  //const double par_damp = setpso.pso_par_w_damp;
  double par_c_cog = 2.05; //2.0
  double par_c_soc = 2.05; //2.0
  const double par_initial_c_cog = 2.5; // 0.5
  const double par_final_c_cog = 0.5; // 2.5
  const double par_initial_c_soc = 0.5; // 2.5
  const double par_final_c_soc = 2.5; // 0.5
  //const double c = 1.193; // value from Akman2018: 1.193
  double prop_objfn_val;
  arma::vec objfn_vals(n_pop); // personal best fitness vector
  GetRNGstate();
  arma::mat P(n_pop, n_vals, arma::fill::randu);
  PutRNGstate();
  arma::mat V = arma::zeros(n_pop,n_vals);
  std::vector<double> param_temp;
  param_temp.resize(n_vals);
  arma::imat nbhood_precursor(n_pop, n_pop, arma::fill::eye); // Matrix with elements of main diagonal set to 1; remeaining entries = 0
  std::vector<std::vector<int> > neighberhood(n_pop); // neighberhood containts for each particle informants
  arma::rowvec real_values(n_vals);
  arma::rowvec gravity(n_vals);
  arma::vec lower_start_bounds(m_lb.size());
  arma::vec upper_start_bounds(m_ub.size());

  // preparation for multithreading
  // =============================
  int hardware_con = std::thread::hardware_concurrency(); // does not work on server?
  int supported_threads = hardware_con == 0 ? 2 : hardware_con;
  int max_amount_of_threads = (supported_threads > n_pop) ? n_pop : hardware_con;
  int quotient = n_pop/max_amount_of_threads;
  int remainder = n_pop % max_amount_of_threads;
  int threads = max_amount_of_threads;
  std::vector<int> indices(threads);
  std::vector<int> sizes(threads);

  // fill sizes
  if(remainder == 0) {
    for(int i = 0; i < threads; i++) {
      sizes[i] = quotient;
    }
  } else {
    for(int i = 0; i < threads; i++) {
      if(i < threads-1) {
        sizes[i] = quotient;
      } else {
        sizes[i] = quotient + remainder;
      }
    }
  }
  indices[0] = 0;
  for(unsigned int i = 1; i < indices.size(); i++) {
    indices[i] = indices[i-1] + sizes[i-1];
  }
  std::vector<arma::Mat<double> > sub_mats(max_amount_of_threads);
  // =============================

  // define borders for start values
  // =============================
  for(size_t i = 0; i < m_lb.size(); i++) {
      lower_start_bounds(i) = m_lb[i];
      upper_start_bounds(i) = m_ub[i];
  }
  // =============================

  // initialize
  // =============================
    for (int i = 0; i < n_pop; i++){
      if(i == 0) {}
      GetRNGstate();
      P.row(i) = lower_start_bounds.t() + (upper_start_bounds.t() - lower_start_bounds.t())%arma::randu(1,n_vals);
      PutRNGstate();
      for(size_t k = 0; k < param_temp.size(); k++) {
            if(P(i,k) > upper_start_bounds(k)) {
            P(i,k) = upper_start_bounds(k);
          } else if(P(i,k) < lower_start_bounds(k)) {
            P(i,k) = lower_start_bounds(k);
          }
          param_temp[k] = P(i,k);
      }
      prop_objfn_val = fctptr(param_temp, odes, model);;
      objfn_vals(i) = prop_objfn_val;
    }

  arma::vec best_vals = objfn_vals; // best_vals = global best solutions
  arma::mat best_vecs = P;
  double global_best_val = objfn_vals.min(); // global best
  arma::rowvec global_best_vec = P.row( objfn_vals.index_min() ); // global best solution
  // =============================

  size_t iter = 0;
  double err = err_tol + 1.0;
  int convergence_check = 0;
  int convergence_check_stop = 0;
  // begin loop
  // =============================
  while (err > err_tol && iter < static_cast<size_t>(n_gen)) {

    // Calculate neighberhood
    // =============================
    if(iter == 0 || convergence_check >= 1) { // 10
      // Randomize indizes of
      // =============================
      int idx1, idx2;
      double best_vals_idx1, best_vals_idx2;
      double objfn_vals_idx1, objfn_vals_idx2;

      for(size_t z = 0; z < 1000; z++) { // z = number of switches
        GetRNGstate();
        idx1 = arma::randi<int>(arma::distr_param(0, n_pop-1));
        PutRNGstate();
        GetRNGstate();
        idx2 = arma::randi<int>(arma::distr_param(0, n_pop-1));
        PutRNGstate();
        P.swap_rows(idx1, idx2);
        best_vecs.swap_rows(idx1, idx2);
        V.swap_rows(idx1, idx2);
        best_vals_idx1 = best_vals(idx1);
        best_vals_idx2 = best_vals(idx2);
        best_vals(idx1) = best_vals_idx2;
        best_vals(idx2) = best_vals_idx1;
        objfn_vals_idx1 = objfn_vals(idx1);
        objfn_vals_idx2 = objfn_vals(idx2);
        objfn_vals(idx1) = objfn_vals_idx2;
        objfn_vals(idx2) = objfn_vals_idx1;
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

    iter++;
    // updating par_w_min, par_w_max, c_sog and par_c_cog
    // =============================
    par_w = par_w_max - iter*(par_w_max - par_w_min)/n_gen;
    par_c_cog = par_initial_c_cog - (par_initial_c_cog - par_final_c_cog) * (iter + 1) / n_gen;
    par_c_soc = par_initial_c_soc - (par_initial_c_soc - par_final_c_soc) * (iter + 1) / n_gen;

    // =============================
    // population loop Nr.1
    // =============================
    for (int i=0; i < n_pop; i++)
    {
      int best_neighberhood_particel;
      arma::vec current_nbhood(neighberhood[i].size());
      for(size_t j = 0; j < neighberhood[i].size(); j++) {
        current_nbhood(j) = best_vals(neighberhood[i][j]);
      }
      int temp_fittness_index = current_nbhood.index_min();
      best_neighberhood_particel = neighberhood[i][temp_fittness_index];
      arma::vec local_best_vec = best_vecs.row(best_neighberhood_particel).t();

      double v = 0.9; 
      double chi;
      chi = 2*v;
      GetRNGstate();
      double omega = par_c_cog*arma::randu() + par_c_soc*arma::randu();
      if(omega < 4.) {
        omega = 4.;
      }
      PutRNGstate();
      chi = chi/std::abs(2. - omega - std::sqrt(omega*(omega -4.) ));


      GetRNGstate();
      V.row(i) = par_w*V.row(i) + par_c_cog*arma::randu(1,n_vals)%(best_vecs.row(i) - P.row(i)) + par_c_soc*arma::randu(1,n_vals)%(local_best_vec.t() - P.row(i));
      PutRNGstate();
      
      V.row(i) = V.row(i)*chi;
      P.row(i) += V.row(i);

      // check if boundaries are violated
      // =============================
        for(int k = 0; k < n_vals; k++){
          if(P(i,k) > upper_start_bounds(k)) {
            P(i,k) = upper_start_bounds(k);
          } else if(P(i,k) < lower_start_bounds(k)) {
            P(i,k) = lower_start_bounds(k);
          }
        }


    }

    std::vector<std::vector<double> > parameter(n_pop);

    for(int i = 0; i < n_pop; i++) {
      parameter[i].resize(n_vals);
      for(int l = 0; l < n_vals; l++) {
        parameter[i][l] = P(i,l);
      }
    }

    
    //objfn_vals.zeros(); // delete!

    //#pragma omp parallel for shared(objfn_vals)
      for(int o = 0; o < n_pop; o++) {
        double current_val = fctptr(parameter[o], odes, model);
        objfn_vals(o) = current_val;
      }


  /*
      for(int i = 0; i < n_pop; i++) {
        double current_val = fctptr(parameter[i], odes, model);
        objfn_vals(i) = current_val;
      }
  */

      /*
      // =============================
      // population loop Nr.2
      for (int i=0; i < n_pop; i++)
      {
      // evaluate updated particle
      // =============================
        for(size_t l = 0; l < param_temp.size(); l++) {
          param_temp[l] = P(i,l);
        }
        prop_objfn_val = fctptr(param_temp, odes, model);
        objfn_vals(i) = prop_objfn_val;
      }
      // =============================
      */

      // population loop Nr.3
      for (int i=0; i < n_pop; i++)
      {
      // Update personal best particle
      // =============================
      if (objfn_vals(i) < best_vals(i))
      {
        best_vals(i) = objfn_vals(i);
        best_vecs.row(i) = P.row(i);
      }
      // =============================
      }
    // =============================

    // Update global best particle
    // =============================
    int min_objfn_val_index = best_vals.index_min();
    double min_objfn_val = best_vals(min_objfn_val_index);

    if (min_objfn_val < global_best_val)
    {
      global_best_val = min_objfn_val;
      global_best_vec = best_vecs.row( min_objfn_val_index );
      convergence_check_stop = 0;
    } else {
      convergence_check = convergence_check + 1;
      convergence_check_stop = convergence_check_stop + 1;
    }
    if(convergence_check_stop > ( (n_gen/n_pop)*10 )) { // 200
      int convergence_threshold = (n_gen/n_pop)*10;
      Rcpp::Rcerr << "No convergence for: " << convergence_threshold << " generations. Optimizing stoped" << std::endl;
      break;
    }
    if(iter % 50 == 0) {
    Rcpp::Rcerr << "global best val" << "\t" << global_best_val << std::endl;
    Rcpp::Rcerr << "Generation number:" << "\t" << iter << std::endl;
    }
    // =============================
    err = global_best_val;
    Rcpp::checkUserInterrupt();
  }
  // =============================

  Rcpp::Rcerr<<"Global best fitness value" <<"\t" << global_best_val <<"\n";
  best_particle_param_values.resize(global_best_vec.size());
  for(size_t i = 0; i < best_particle_param_values.size(); i++) {
      best_particle_param_values[i] = global_best_vec(i);
  }
  return global_best_val;
}
