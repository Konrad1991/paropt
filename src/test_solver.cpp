#include "solver.hpp"

static int check_retval(void *returnvalue, const char *funcname, int opt);

struct usr_data{
  Rcpp::Function ode_systen;
  std::vector<double> parameter;
  std::vector<double> parameter_time;
  std::vector<int> parameter_cut_idx;
};

double solver_return_data(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, std::vector<std::vector<double> > &nc) {

    std::vector<double> init_state = solv_param_struc.init_state;
    std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
    std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
    std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
    std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
    std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
    Rcpp::NumericVector integration_times = solv_param_struc.integration_times;

    nc.resize(hs_cut_idx_vec.size());
    for(size_t i = 0; i < nc.size(); i++) {
      nc[i].resize(integration_times.size());
    }

    for(size_t i = 0; i < nc.size(); i++) {
      nc[i][0] = solv_param_struc.init_state[i];
    }

    // Begin Solver
    int NEQ = hs_cut_idx_vec.size();
    realtype reltol, t, tout;
    N_Vector y, abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    int retval, retvalr, iout;

    y = abstol = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;

    y = N_VNew_Serial(NEQ);
    if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
    abstol = N_VNew_Serial(NEQ);
    if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

    for (int i = 0; i < NEQ; ++i) {
      NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
      NV_Ith_S(y, i) = solv_param_struc.init_state[i];
    }

    reltol = solv_param_struc.reltol;

    cvode_mem = CVodeCreate(CV_BDF);
    if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    double sum_of_least_squares = 0.;

    struct usr_data my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
    void*ptr_to_my_ode_system = &my_ode_system;
    retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
    if (check_retval((void *)cvode_mem, "CVodeSetUserData", 0)) return(1);

    retval = CVodeInit(cvode_mem, wrapper_ode_system, integration_times[0], y);
    if (check_retval(&retval, "CVodeInit", 1)) return(1);

    retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

    A = SUNDenseMatrix(NEQ, NEQ);
    if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

    LS = SUNLinSol_Dense(y, A);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

    int CVodetmpcount=0;
    double return_time;
    float return_steps=1.;

    std::vector<double> temp_measured(init_state.size());
    int num_states = init_state.size();

    for ( int ti = 1; ti < integration_times.size(); ti++) {
       for (int j = 1; j <= return_steps; j++) {
             return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
             retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
           for (size_t n = 0; n < NV_LENGTH_S(y); n++) {
             nc[n][ti] = NV_Ith_S(y,n);
           }
       }
    }

    N_VDestroy(y);
    N_VDestroy(abstol);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

   return 0.;
}

double ADAMS_return_data(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, std::vector<std::vector<double> > &nc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  Rcpp::NumericVector integration_times = solv_param_struc.integration_times;

  nc.resize(hs_cut_idx_vec.size());
  for(size_t i = 0; i < nc.size(); i++) {
    nc[i].resize(integration_times.size());
  }

  for(size_t i = 0; i < nc.size(); i++) {
    nc[i][0] = solv_param_struc.init_state[i];
  }

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t, tout;
   N_Vector y, abstol;
   void *cvode_mem;
   int retval, retvalr, iout;

   y = abstol = NULL;
   cvode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   cvode_mem = CVodeCreate(CV_ADAMS);
   if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

   double sum_of_least_squares = 0.;

     struct usr_data my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
     if (check_retval((void *)cvode_mem, "CVodeSetUserData", 0)) return(1);

     retval = CVodeInit(cvode_mem, wrapper_ode_system, integration_times[0], y);
     if (check_retval(&retval, "CVodeInit", 1)) return(1);

     retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
     if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

     retval = CVDiag(cvode_mem);
     if(check_retval(&retval, "CVDiag", 1)) return(1);

     int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     int num_states = init_state.size();

       for ( int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
                   for (size_t n = 0; n < NV_LENGTH_S(y); n++) {
                     nc[n][ti] = NV_Ith_S(y,n);
                   }
           }
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       CVodeFree(&cvode_mem);

   return 0.;
}

double ERK_return_data(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, std::vector<std::vector<double> > &nc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  Rcpp::NumericVector integration_times = solv_param_struc.integration_times;

  nc.resize(hs_cut_idx_vec.size());
  for(size_t i = 0; i < nc.size(); i++) {
    nc[i].resize(integration_times.size());
  }

  for(size_t i = 0; i < nc.size(); i++) {
    nc[i][0] = solv_param_struc.init_state[i];
  }

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t, tout;
   N_Vector y, abstol;
   void *arkode_mem;
   int retval, retvalr, iout;

   y = abstol = NULL;
   arkode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   arkode_mem = ERKStepCreate(wrapper_ode_system, integration_times[0], y);
   if (check_retval((void *)arkode_mem, "ERKCreate", 0)) return(1);

   double sum_of_least_squares = 0.;

     struct usr_data my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = ERKStepSetUserData(arkode_mem, ptr_to_my_ode_system);
     if (check_retval(&retval, "ERKStepSetUserData", 1)) return(1);

     retval = ERKStepSVtolerances(arkode_mem, reltol, abstol);
     if (check_retval(&retval, "ERKStepSVtolerances", 1)) return(1);

     int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     int num_states = init_state.size();

       for ( int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = ERKStepEvolve(arkode_mem, return_time, y, &t, ARK_NORMAL);
                   for (size_t n = 0; n < NV_LENGTH_S(y); n++) {
                     nc[n][ti] = NV_Ith_S(y,n);
                   }
           }
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       ERKStepFree(&arkode_mem);

   return 0.;
}

double ARK_return_data(std::vector<double> &param_combi_start, SEXP ode_system, time_state_information solv_param_struc, std::vector<std::vector<double> > &nc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  Rcpp::NumericVector integration_times = solv_param_struc.integration_times;

  nc.resize(hs_cut_idx_vec.size());
  for(size_t i = 0; i < nc.size(); i++) {
    nc[i].resize(integration_times.size());
  }

  for(size_t i = 0; i < nc.size(); i++) {
    nc[i][0] = solv_param_struc.init_state[i];
  }

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t, tout;
   N_Vector y, abstol;
   SUNMatrix A = NULL;
   SUNLinearSolver LS = NULL;
   void *arkode_mem;
   int retval, retvalr, iout;

   y = abstol = NULL;
   arkode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   A = SUNDenseMatrix(NEQ, NEQ);
   if (check_retval((void *)A, "SUNDenseMatrix", 0)) return 1;
   LS = SUNLinSol_Dense(y, A);
   if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return 1;

   // right hand side => y' = f(t,y); f(t,y) = f_E + f_I
   // f_E = explizit Part; f_I = implicit part
   // f_E is zero => fully implicit
   arkode_mem = ARKStepCreate(NULL, wrapper_ode_system, integration_times[0], y);
   if (check_retval((void *)arkode_mem, "ARKCreate", 0)) return(1);

   double sum_of_least_squares = 0.;

     struct usr_data my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = ARKStepSetUserData(arkode_mem, ptr_to_my_ode_system);
     if (check_retval(&retval, "ARKStepSetUserData", 1)) return(1);

     retval = ARKStepSVtolerances(arkode_mem, reltol, abstol);
     if (check_retval(&retval, "ARKStepSStolerances", 1)) return(1);

     retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
     if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return(1);

     int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     int num_states = init_state.size();

       for ( int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = ARKStepEvolve(arkode_mem, return_time, y, &t, ARK_NORMAL);
                   for (size_t n = 0; n < NV_LENGTH_S(y); n++) {
                     nc[n][ti] = NV_Ith_S(y,n);
                   }
           }
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       SUNMatDestroy(A);
       SUNLinSolFree(LS);
       ARKStepFree(&arkode_mem);

   return 0.;
}

// [[Rcpp::export]]
std::vector<std::vector<double> > test_solver(Rcpp::NumericVector integration_times,
                   SEXP ode_system, double relative_tolerance,
                   Rcpp::NumericVector absolute_tolerances, std::string start,
                   std::string lower, std::string upper, std::string states,
                  std::string solvertype) {

    std::vector<int> params_cut_idx_vec;
    std::vector<double> params_time_combi_vec;
    std::vector<double> param_combi_start;
    std::vector<double> param_combi_lb;
    std::vector<double> param_combi_ub;
    std::vector<std::string> header_parameter;

    Import_Parameter(start, lower, upper, params_cut_idx_vec, params_time_combi_vec,
    param_combi_start, param_combi_lb, param_combi_ub, header_parameter);

    std::vector<int> hs_cut_idx_vec;
    std::vector<double> hs_time_combi_vec;
    std::vector<double> hs_harvest_state_combi_vec;
    std::vector<std::string> header_states;
    Import_states(states, hs_cut_idx_vec, hs_time_combi_vec, hs_harvest_state_combi_vec, header_states);

    int tmpcount=0;

    std::vector<double> init_state ( hs_cut_idx_vec.size() );
    for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
     init_state[i] = hs_harvest_state_combi_vec[tmpcount];
     tmpcount += hs_cut_idx_vec[i];
    }

    time_state_information param_model;

    param_model.init_state = init_state;
    param_model.par_times = params_time_combi_vec;
    param_model.param_idx_cuts = params_cut_idx_vec;
    param_model.state_measured = hs_harvest_state_combi_vec;
    param_model.state_times = hs_time_combi_vec;
    param_model.state_idx_cut = hs_cut_idx_vec;
    param_model.integration_times = integration_times;
    param_model.reltol = relative_tolerance;
    param_model.absolute_tolerances = absolute_tolerances;

    std::vector<std::vector<double> > nc;

    if(solvertype == "bdf") {
    solver_return_data(param_combi_start, ode_system, param_model, nc);}
    else if(solvertype == "ADAMS") {
    ADAMS_return_data(param_combi_start, ode_system, param_model, nc);}
    else if(solvertype == "ERK") {
    ERK_return_data(param_combi_start, ode_system, param_model, nc);}
    else if(solvertype == "ARK") {
    ARK_return_data(param_combi_start, ode_system, param_model, nc);}
    else {
      Rcpp::Rcerr << "Unknown Solvertyp" << std::endl;
    }
    return nc;
}

// [[Rcpp::export]]
std::vector<std::vector<double> > test_solve_ode_system(Rcpp::NumericVector integration_times,
                   SEXP ode_system, double relative_tolerance,
                   Rcpp::NumericVector absolute_tolerances,
                    std::string start, std::string states, std::string solvertype) {

  std::vector<std::vector<double> > nc;
  std::vector<int> params_cut_idx_vec;
  std::vector<double> params_time_combi_vec;
  std::vector<double> param_combi_start;
  std::vector<std::string> header_parameter;

  Import_start_parameter(start, params_cut_idx_vec, params_time_combi_vec,  param_combi_start, header_parameter);

  std::vector<int> hs_cut_idx_vec;
  std::vector<double> hs_time_combi_vec;
  std::vector<double> hs_harvest_state_combi_vec;
  std::vector<std::string> header_states;

  Import_states(states, hs_cut_idx_vec, hs_time_combi_vec, hs_harvest_state_combi_vec, header_states);

  int tmpcount=0;

  std::vector<double> init_state ( hs_cut_idx_vec.size() );
  for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
   init_state[i] = hs_harvest_state_combi_vec[tmpcount];
   tmpcount += hs_cut_idx_vec[i];
  }

  if(init_state.size() > absolute_tolerances.length()) {
   Rcpp::stop("\nERROR: absolute tolerances not defined for each state");
   exit (EXIT_FAILURE);
  }

  if(init_state.size() < absolute_tolerances.length()) {
   Rcpp::stop("\nERROR: dimension error for absolute tolerances");
     exit (EXIT_FAILURE);
   }

   std::vector<double>::iterator max_time_param_vector = std::max_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
   std::vector<double>::iterator min_time_param_vector = std::min_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
   std::vector<double>::iterator max_time_harvest_vector = std::max_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());
   std::vector<double>::iterator min_time_harvest_vector = std::min_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());

   if(*max_time_param_vector != 0) {
     if(*max_time_param_vector > *max_time_harvest_vector) {
       Rcpp::stop("\nERROR: Maximum of timevector of parameter larger then corresponding timepoint of state vector");
       exit (EXIT_FAILURE);
     } else if(*max_time_param_vector < *max_time_harvest_vector) {
       Rcpp::stop("\nERROR: Maximum of timevector of parameter smaller then corresponding timepoint of state vector");
       exit (EXIT_FAILURE);
     }
     if(*min_time_param_vector > *min_time_harvest_vector) {
       Rcpp::stop("\nERROR: Minimum of timevector of parameter larger then corresponding timepoint of state vector");
       exit (EXIT_FAILURE);
     } else if(*min_time_param_vector < *min_time_harvest_vector) {
       Rcpp::stop("\nERROR: Minimum of timevector of parameter smaller then corresponding timepoint of state vector");
       exit (EXIT_FAILURE);
     }
   } else {
   // all parameter are scalare
   }
   std::vector<double> integration_time_assumption(hs_cut_idx_vec[0]);
   for(size_t i = 0; i < hs_cut_idx_vec[0];i++) {
     integration_time_assumption[i] = hs_time_combi_vec[i];
   }
   if(integration_times.length() > integration_time_assumption.size()) {
     Rcpp::stop("\nERROR: integration_times must not be larger than time of state input");
     exit (EXIT_FAILURE);
   }

   if(TYPEOF(ode_system) != CLOSXP) {
     Rcpp::stop("\nERROR: type of odesystem should be closure");
     exit (EXIT_FAILURE);
   }

   time_state_information param_model;

   param_model.init_state = init_state;
   param_model.par_times = params_time_combi_vec;
   param_model.param_idx_cuts = params_cut_idx_vec;
   param_model.state_measured = hs_harvest_state_combi_vec;
   param_model.state_times = hs_time_combi_vec;
   param_model.state_idx_cut = hs_cut_idx_vec;
   param_model.integration_times = integration_times;
   param_model.reltol = relative_tolerance;
   param_model.absolute_tolerances = absolute_tolerances;

   if(solvertype == "bdf") {
   solver_return_data(param_combi_start, ode_system, param_model, nc);}
   else if(solvertype == "ADAMS") {
   ADAMS_return_data(param_combi_start, ode_system, param_model, nc);}
   else if(solvertype == "ERK") {
   ERK_return_data(param_combi_start, ode_system, param_model, nc);}
   else if(solvertype == "ARK") {
   ARK_return_data(param_combi_start, ode_system, param_model, nc);}
   else {
     Rcpp::Rcerr << "Unknown Solvertyp" << std::endl;
   }

return nc;
}


static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    Rcpp::Rcerr << "SUNDIALS_ERROR:" << " " << funcname << " " << "failed - returned NULL pointer" << std::endl;
    return(1); }
  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      Rcpp::Rcerr << "SUNDIALS_ERROR:" << " " << funcname << " " << "failed with retval = " << " " << *retval << std::endl;
      return(1); }}
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    Rcpp::Rcerr << "MEMORY_ERROR:" << " " << funcname << " " << "failed - returned NULL pointer" << std::endl;
    return(1); }

  return(0);
}
