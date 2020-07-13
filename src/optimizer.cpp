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

#include "optimizer.hpp"

Optimizer::Optimizer(
        std::vector<double> t_start,
        std::vector<double> t_lb,
        std::vector<double> t_ub,
        settingsPSO t_set_pso,
        time_state_information t_model,
        double (*fctptr2)(std::vector<double> &param_vec, SEXP ode_system, time_state_information &model),
        SEXP t_odes) {

    m_start.resize(t_start.size());
    m_lb.resize(t_lb.size());
    m_ub.resize(t_ub.size());
    for(size_t i = 0; i < t_lb.size(); i++) {
        m_start[i] = t_start[i];
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

void Optimizer::get_best_particle_param_values(std::vector<double> & temp) {
    temp.resize(best_particle_param_values.size());
    for(size_t   i = 0; i < temp.size(); i++) {
        temp[i] = best_particle_param_values[i];
    }
}

// ========================================================================
// ========================================================================
// ========================================================================

double Optimizer::pso() { // (labled with ! need check)

  // define parameters
  // =============================
  bool classical_pso = true;
  const int n_vals = m_lb.size();
  const double err_tol = setpso.err_tol;
  const int n_pop = setpso.pso_n_pop;
  const int n_gen = setpso.pso_n_gen;
  const int inertia_method = 1; //1 or 2
  double par_w = setpso.pso_par_initial_w;
  const double par_w_max = setpso.pso_par_w_max;
  const double par_w_min = setpso.pso_par_w_min;
  const double par_damp = setpso.pso_par_w_damp;
  const int velocity_method = 2; //1 or 2
  double par_c_cog = 2.05; //2.0
  double par_c_soc = 2.05; //2.0
  const double par_initial_c_cog = 2.5; // 0.5
  const double par_final_c_cog = 0.5; // 2.5
  const double par_initial_c_soc = 0.5; // 2.5
  const double par_final_c_soc = 2.5; // 0.5
  const double c = 1.193; // value from Akman2018: 1.193
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
  arma::vec start_values(m_start.size());
  arma::vec lower_start_bounds(m_lb.size());
  arma::vec upper_start_bounds(m_ub.size());
  // =============================

  // define borders for start values
  // =============================
  for(size_t i = 0; i < m_lb.size(); i++) {
      start_values(i) = m_start[i];
      lower_start_bounds(i) = m_lb[i];
      upper_start_bounds(i) = m_ub[i];
  }
  // =============================

  // initialize
  // =============================

  if(classical_pso == true) {
    for (int i = 0; i < n_pop; i++){
      if(i == 0) {
      P.row(i) = start_values.t();}
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
  } else {
      for (int i = 0; i < n_pop; i++){
      real_values = lower_start_bounds.t() + (upper_start_bounds.t() - lower_start_bounds.t())%P.row(i);
        for(size_t k = 0; k < param_temp.size(); k++) {
            param_temp[k] = real_values(k);
        }
      prop_objfn_val = fctptr(param_temp, odes, model);;
      objfn_vals(i) = prop_objfn_val;
      }
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
    if (inertia_method == 1) {
        //par_w = par_w_min + (par_w_max - par_w_min) * (iter + 1) / n_gen; //!
        // par_w = par_w_max - (par_w_max - par_w_min) * (iter + 1) / n_gen; //!
        par_w = par_w_max - iter*(par_w_max - par_w_min)/n_gen;
        //par_w = 0.721; // Akman2018
        //par_w = arma::randu(); //!
    } else {
        par_w *= par_damp;
    }
    if (velocity_method == 2)
    {
        par_c_cog = par_initial_c_cog - (par_initial_c_cog - par_final_c_cog) * (iter + 1) / n_gen;
        par_c_soc = par_initial_c_soc - (par_initial_c_soc - par_final_c_soc) * (iter + 1) / n_gen;
    }
    // =============================
    // population loop
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

      if(classical_pso == true) {
            GetRNGstate();
            V.row(i) = par_w*V.row(i) + par_c_cog*arma::randu(1,n_vals)%(best_vecs.row(i) - P.row(i)) + par_c_soc*arma::randu(1,n_vals)%(local_best_vec.t() - P.row(i));
            PutRNGstate();
            P.row(i) += V.row(i);
      } else {
        if(best_neighberhood_particel != i) {
          //gravity = P.row(i) + c*(best_vecs.row(i) + local_best_vec.t() - 2.*P.row(i))/3.;
          gravity = P.row(i) + c*(best_vecs.row(i) + local_best_vec.t() - 2.*P.row(i))/3.;
        } else {
          //gravity = P.row(i) + c*(best_vecs.row(i) - P.row(i))/2.;
          gravity = P.row(i) + c*(best_vecs.row(i) + local_best_vec.t() - 2.*P.row(i))/3.;
        }
        GetRNGstate();
        double radius = arma::norm(P.row(i) - gravity);
        PutRNGstate();
        GetRNGstate();
        arma::rowvec temp_p = arma::randu(1,n_vals);
        PutRNGstate();
        GetRNGstate();
        temp_p = arma::randu()*radius*temp_p/arma::norm(temp_p) + gravity;
        PutRNGstate();
        // Update velocity
        V.row(i) = par_w*V.row(i) + (temp_p - P.row(i));
        // Update position
        //P.row(i) = P.row(i) + V.row(i);
        P.row(i) = par_w*V.row(i) + temp_p;
      }

      // check if boundaries are violated
      // =============================
      if(classical_pso == true) {
        for(int k = 0; k < n_vals; k++){
          if(P(i,k) > upper_start_bounds(k)) {
            P(i,k) = upper_start_bounds(k);
          } else if(P(i,k) < lower_start_bounds(k)) {
            P(i,k) = lower_start_bounds(k);
          }
        }
      } else {
        for(int k = 0; k < n_vals; k++) {
          if(P(i,k) > 1.) {
            P(i,k) = 1;
          } else if(P(i,k) < 0.) {
            P(i,k) = 0.;
          }
        }
      }
      // =============================

      // evaluate updated particle
      // =============================
      if(classical_pso == true) {
        for(size_t l = 0; l < param_temp.size(); l++) {
          param_temp[l] = P(i,l);
        }
        prop_objfn_val = fctptr(param_temp, odes, model);
        objfn_vals(i) = prop_objfn_val;
      } else {
        real_values = lower_start_bounds.t() + (upper_start_bounds.t() - lower_start_bounds.t())%P.row(i);

        for(size_t l = 0; l < param_temp.size(); l++) {
          if(real_values(l) > upper_start_bounds(l)) {
            real_values(l) = upper_start_bounds(l);
          } else if(real_values(l) < lower_start_bounds(l)) {
            real_values(l) = lower_start_bounds(l);
          }
            param_temp[l] = real_values(l);
        }
      prop_objfn_val = fctptr(param_temp, odes, model);
      objfn_vals(i) = prop_objfn_val;
      }
      // =============================

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
    if(convergence_check_stop > 200) {
      Rcpp::Rcerr << "No convergence for 200 generations. Optimizing stoped" << std::endl;
      break;
    }
    Rcpp::Rcerr << "global best val" << "\t" << global_best_val << std::endl;
    Rcpp::Rcerr << "Generation number:" << "\t" << iter << std::endl;
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
