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

#ifndef OPTIMIZER
#define OPTIMIZER

#include "header.hpp"

class Optimizer {
private:
    std::vector<double> m_lb;
    std::vector<double> m_ub;
    settingsPSO setpso;
    time_state_information model;
    double best_particle_objective_value;
    std::vector<double> best_particle_param_values;
    double (*fctptr)(std::vector<double> &param_vec, SEXP ode_system, time_state_information &model);
    SEXP odes;
public:
    Optimizer (
        std::vector<double> t_lb,
        std::vector<double> t_ub,
        settingsPSO t_set_pso,
        time_state_information t_model,
        double (*fctptr)(std::vector<double> &param_vec, SEXP ode_system, time_state_information &model),
        SEXP t_odes);
    ~ Optimizer () {}

    //double pso(std::vector<double> m_lb,std::vector<double> m_ub, settingsPSO setpso, time_state_information model, double (*fctptr)(time_state_information model, //std::vector<double> param_vec), std::vector<double> best_particle_param_values);

    double pso();

    //double get_best_particle_objective_value () const {
      //  pso(m_lb, m_ub, setpso, model, fctptr, &best_particle_param_values);
        //return best_particle_objective_value;};
     void get_best_particle_param_values (std::vector<double> & temp);
};

#endif // OPTIMIZER
