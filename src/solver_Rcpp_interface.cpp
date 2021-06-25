/* !!revision!!
remove RKs?

error check if error is inf, NA etc.
*/

/*
BSD 3-Clause License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

The names of its contributors may not be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "solver_Rcpp_interface.hpp"

typedef std::vector<double> VD;
typedef std::vector<int> VI;
typedef std::vector<std::vector<double> > MD;
typedef std::vector<std::vector<int> > MI;
typedef std::vector<std::string> VS;
typedef Rcpp::DataFrame DF;

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

static int check_retval_Rcpp_interface(void *returnvalue, const char *funcname, int opt);

bool double_diff_Rcpp_interface(double x, double y) {
  bool equal = false;
  double tolerance = 0.001;
  if(std::abs(x - y) < tolerance) {
    equal = true;
  }
  return equal;
}

// Rcpp::Rcerr is not thread safe!!!
void own_error_handler_Rcpp_interface(int error_code, const char *module, const char *function, char *msg, void *usr_data) {
  if(error_code < 0) {
    //Rcpp::Rcerr << "Error:" << " " << error_code << " " << msg << " " << "In module:" << " " << module << " " << "and in fct:" << " " << function << std::endl;
  }
  else if(error_code == 0) {
    //Rcpp::warning("\nSundials: function encounters error: %s", function);
  } else {
    //Rcpp::warning("\nSundials-warning during integration with: error_code: %i", error_code);
  }
}


struct usr_data_Rcpp_interface{
  OS ode_system; //Rcpp::XPtr<OS> ode_system;
  std::vector<double> parameter;
  std::vector<double> parameter_time;
  std::vector<int> parameter_cut_idx;
};

int wrapper_ode_system_Rcpp_interface(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  // cast pointer to structure and store elements
  struct usr_data_Rcpp_interface *my_ode_system = (struct usr_data_Rcpp_interface*)user_data;
  OS odes;
  odes = (*my_ode_system).ode_system; //*(*my_ode_system).ode_system;
  std::vector<double> params = (*my_ode_system).parameter;
  std::vector<double> params_time = (*my_ode_system).parameter_time;
  std::vector<int> params_cut_idx = (*my_ode_system).parameter_cut_idx;

  // current state values ( = y)
  std::vector<double> current_states(NV_LENGTH_S(y));
  for(unsigned int i = 0; i < current_states.size(); i++) {
    current_states[i] = NV_Ith_S(y,i);
  }

  // interpolate Parameter if necessary
  std::vector<double> parameter_input;
  params_sort(t, parameter_input, params, params_time, params_cut_idx);

  // extract time
  double time = t;

  odes(time, parameter_input, current_states);

  for(unsigned int i = 0; i < current_states.size(); i++) {
    NV_Ith_S(ydot, i) = current_states[i];
  }

  return 0;
}

double solver_bdf_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t;// tout;
   N_Vector y, abstol;
   SUNMatrix A;
   SUNLinearSolver LS;
   void *cvode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   A = NULL;
   LS = NULL;
   cvode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   cvode_mem = CVodeCreate(CV_BDF);
   if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_cvode_mem = &cvode_mem;
   void* ptr_to_nothing = &cvode_mem;
   retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

   struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
   void*ptr_to_my_ode_system = &my_ode_system;
   retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
   if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeSetUserData", 0)) return(1);

   retval = CVodeInit(cvode_mem, wrapper_ode_system_Rcpp_interface, integration_times[0], y);
   if (check_retval_Rcpp_interface(&retval, "CVodeInit", 1)) return(1);

   retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
   if (check_retval_Rcpp_interface(&retval, "CVodeSVtolerances", 1)) return(1);

   A = SUNDenseMatrix(NEQ, NEQ);
   if(check_retval_Rcpp_interface((void *)A, "SUNDenseMatrix", 0)) return(1);

   LS = SUNLinSol_Dense(y, A);
   if(check_retval_Rcpp_interface((void *)LS, "SUNLinSol_Dense", 0)) return(1);

   retval = CVodeSetLinearSolver(cvode_mem, LS, A);
   if(check_retval_Rcpp_interface(&retval, "CVodeSetLinearSolver", 1)) return(1);

   //int CVodetmpcount=0;
   double return_time;
   float return_steps=1.;
   std::vector<double> temp_measured(init_state.size());
   //int num_states = init_state.size();

       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       //sum_of_least_squares += std::abs((1./temp_measured[n])*std::abs(NV_Ith_S(y,n) - temp_measured[n]) );
                       }
                   }
                   if(retval < 0) {
                     sum_of_least_squares = 1.79769e+308; //1000000.; //1.79769e+308
                     break;
                   }
           }
       }

    N_VDestroy(y);
    N_VDestroy(abstol);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_bdf_save_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface solv_param_struc, Rcpp::NumericMatrix &DF) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t; // tout;
   N_Vector y, abstol;
   SUNMatrix A;
   SUNLinearSolver LS;
   void *cvode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   A = NULL;
   LS = NULL;
   cvode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   cvode_mem = CVodeCreate(CV_BDF);
   if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_cvode_mem = &cvode_mem;
   void* ptr_to_nothing = &cvode_mem;
   retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeSetUserData", 0)) return(1);

     retval = CVodeInit(cvode_mem, wrapper_ode_system_Rcpp_interface, integration_times[0], y);
     if (check_retval_Rcpp_interface(&retval, "CVodeInit", 1)) return(1);

     retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "CVodeSVtolerances", 1)) return(1);

     A = SUNDenseMatrix(NEQ, NEQ);
     if(check_retval_Rcpp_interface((void *)A, "SUNDenseMatrix", 0)) return(1);

     LS = SUNLinSol_Dense(y, A);
     if(check_retval_Rcpp_interface((void *)LS, "SUNLinSol_Dense", 0)) return(1);

     retval = CVodeSetLinearSolver(cvode_mem, LS, A);
     if(check_retval_Rcpp_interface(&retval, "CVodeSetLinearSolver", 1)) return(1);

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();


          for(int i = 0; i < NV_LENGTH_S(y); i++) {
              DF(0 , i) = NV_Ith_S(y,i);
          }

          int counter = 1;
            for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
                for (int j = 1; j <= return_steps; j++) {
                    return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
                    retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL); //
                        for (int n = 0; n < NV_LENGTH_S(y); n++) {
                          DF(counter, n) = NV_Ith_S(y,n);
                            temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                            if(std::isnan(temp_measured[n])) { }
                            else {
                            sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                            //sum_of_least_squares += std::abs((1./temp_measured[n])*std::abs(NV_Ith_S(y,n) - temp_measured[n]) );
                            //Rcpp::Rcerr << NV_Ith_S(y,n) << "\t" << temp_measured[n] << "\t" << return_time << std::endl;
                            }
                        }
                        if (check_retval_Rcpp_interface(&retval, "CVode", 1)) {
                            break;}
                }
                counter++;
            }

       N_VDestroy(y);
       N_VDestroy(abstol);
       CVodeFree(&cvode_mem);
       SUNLinSolFree(LS);
       SUNMatDestroy(A);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_adams_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t;// tout;
   N_Vector y, abstol;
   void *cvode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   cvode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   cvode_mem = CVodeCreate(CV_ADAMS);
   if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_cvode_mem = &cvode_mem;
   void* ptr_to_nothing = &cvode_mem;
   retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeSetUserData", 0)) return(1);

     retval = CVodeInit(cvode_mem, wrapper_ode_system_Rcpp_interface, integration_times[0], y);
     if (check_retval_Rcpp_interface(&retval, "CVodeInit", 1)) return(1);

     retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "CVodeSVtolerances", 1)) return(1);

     retval = CVDiag(cvode_mem);
     if(check_retval_Rcpp_interface(&retval, "CVDiag", 1)) return(1);

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();

       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       }
                   }
                   if(retval < 0) {
                     sum_of_least_squares = 1.79769e+308;//1000000.;
                     break;
                   }
           }
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       CVodeFree(&cvode_mem);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_adams_save_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface solv_param_struc,  Rcpp::NumericMatrix &DF) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t; // tout;
   N_Vector y, abstol;
   void *cvode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   cvode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   cvode_mem = CVodeCreate(CV_ADAMS);
   if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_cvode_mem = &cvode_mem;
   void* ptr_to_nothing = &cvode_mem;
   retval = CVodeSetErrHandlerFn(cvode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = CVodeSetUserData(cvode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)cvode_mem, "CVodeSetUserData", 0)) return(1);

     retval = CVodeInit(cvode_mem, wrapper_ode_system_Rcpp_interface, integration_times[0], y);
     if (check_retval_Rcpp_interface(&retval, "CVodeInit", 1)) return(1);

     retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "CVodeSVtolerances", 1)) return(1);

     retval = CVDiag(cvode_mem);
     if(check_retval_Rcpp_interface(&retval, "CVDiag", 1)) return(1);

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();

     for(int i = 0; i < NV_LENGTH_S(y); i++) {
         DF(0 , i) = NV_Ith_S(y,i);
     }

     int counter = 1;
       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = CVode(cvode_mem, return_time, y, &t, CV_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                     DF(counter, n) = NV_Ith_S(y,n);
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       //Rcpp::Rcerr << NV_Ith_S(y,n) << "\t" << temp_measured[n] << "\t" << return_time << std::endl;
                       }
                   }
                   if (check_retval_Rcpp_interface(&retval, "CVode", 1)) {
                       break;}

           }
           counter++;
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       CVodeFree(&cvode_mem);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_erk_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t; // tout;
   N_Vector y, abstol;
   void *arkode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   arkode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   arkode_mem = ERKStepCreate(wrapper_ode_system_Rcpp_interface, integration_times[0], y);
   if (check_retval_Rcpp_interface((void *)arkode_mem, "ERKCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_arkode_mem = &arkode_mem;
   void* ptr_to_nothing = &arkode_mem;
   retval = ARKStepSetErrHandlerFn(arkode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = ERKStepSetUserData(arkode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)arkode_mem, "ERKStepSetUserData", 0)) return(1);

     retval = ERKStepSVtolerances(arkode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "ERKStepSVtolerances", 1)) return(1);

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();

       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = ERKStepEvolve(arkode_mem, return_time, y, &t, ARK_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       }
                   }
                   if(retval < 0) {
                     sum_of_least_squares = 1.79769e+308; //1000000.;
                     break;
                   }
           }
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       ERKStepFree(&arkode_mem);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_erk_save_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface solv_param_struc,  Rcpp::NumericMatrix &DF) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t; // tout;
   N_Vector y, abstol;
   void *arkode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   arkode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   arkode_mem = ERKStepCreate(wrapper_ode_system_Rcpp_interface, integration_times[0], y);
   if (check_retval_Rcpp_interface((void *)arkode_mem, "ERKCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_arkode_mem = &arkode_mem;
   void* ptr_to_nothing = &arkode_mem;
   retval = ARKStepSetErrHandlerFn(arkode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = ERKStepSetUserData(arkode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)arkode_mem, "ERKStepSetUserData", 0)) return(1);

     retval = ERKStepSVtolerances(arkode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "ERKStepSVtolerances", 1)) return(1);

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();

     std::ofstream myfile;
     for(int i = 0; i < NV_LENGTH_S(y); i++) {
         DF(0 , i) = NV_Ith_S(y,i);
     }

     int counter = 1;

       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = ERKStepEvolve(arkode_mem, return_time, y, &t, ARK_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                      DF(counter, n) = NV_Ith_S(y,n);
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       //Rcpp::Rcerr << NV_Ith_S(y,n) << "\t" << temp_measured[n] << "\t" << return_time << std::endl;
                       }
                   }
                   if (check_retval_Rcpp_interface(&retval, "CVode", 1)) {
                       break;}
           }
           counter++;
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       ERKStepFree(&arkode_mem);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_ark_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface &solv_param_struc) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t; // tout;
   N_Vector y, abstol;
   SUNMatrix A = NULL;
   SUNLinearSolver LS = NULL;
   void *arkode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   arkode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   A = SUNDenseMatrix(NEQ, NEQ);
   if (check_retval_Rcpp_interface((void *)A, "SUNDenseMatrix", 0)) return 1;
   LS = SUNLinSol_Dense(y, A);
   if (check_retval_Rcpp_interface((void*)LS, "SUNLinSol_Dense", 0)) return 1;

   // right hand side => y' = f(t,y); f(t,y) = f_E + f_I
   // f_E = explizit Part; f_I = implicit part
   // f_E is zero => fully implicit
   arkode_mem = ARKStepCreate(NULL, wrapper_ode_system_Rcpp_interface, integration_times[0], y);
   if (check_retval_Rcpp_interface((void *)arkode_mem, "ARKStepCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_arkode_mem = &arkode_mem;
   void* ptr_to_nothing = &arkode_mem;
   retval = ARKStepSetErrHandlerFn(arkode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = ARKStepSetUserData(arkode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)arkode_mem, "ARKStepSetUserData", 0)) return(1);

     retval = ARKStepSVtolerances(arkode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "ARKStepSVtolerances", 1)) return(1);

     retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
     if (check_retval_Rcpp_interface(&retval, "ARKStepSetLinearSolver", 1)) return(1);

     //retval = ARKStepSetLinear(arkode_mem, 1);
     //if (check_retval_Rcpp_interface(&retval, "ARKStepSetLinear", 1)) return 1;

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();

       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = ARKStepEvolve(arkode_mem, return_time, y, &t, ARK_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       }
                   }
                   if(retval < 0) {
                     sum_of_least_squares = 1.79769e+308; //1000000.;
                     break;
                   }
           }
       }

       N_VDestroy(y);
       N_VDestroy(abstol);
       SUNMatDestroy(A);
       SUNLinSolFree(LS);
       ARKStepFree(&arkode_mem);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

double solver_ark_save_Rcpp_interface(std::vector<double> &param_combi_start, OS ode_system, time_state_information_Rcpp_interface solv_param_struc,  Rcpp::NumericMatrix &DF) {

  std::vector<double> init_state = solv_param_struc.init_state;
  std::vector<double> params_time_combi_vec = solv_param_struc.par_times;
  std::vector<int> params_cut_idx_vec = solv_param_struc.param_idx_cuts;
  std::vector<double> hs_harvest_state_combi_vec = solv_param_struc.state_measured;
  std::vector<double> hs_time_combi_vec = solv_param_struc.state_times;
  std::vector<int> hs_cut_idx_vec = solv_param_struc.state_idx_cut;
  std::vector<double> integration_times = solv_param_struc.integration_times;

    // Begin Solver
   int NEQ = hs_cut_idx_vec.size();
   realtype reltol, t; // tout;
   N_Vector y, abstol;
   SUNMatrix A = NULL;
   SUNLinearSolver LS = NULL;
   void *arkode_mem;
   int retval; // retvalr, iout;

   y = abstol = NULL;
   arkode_mem = NULL;

   y = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)y, "N_VNew_Serial", 0)) return(1);
   abstol = N_VNew_Serial(NEQ);
   if (check_retval_Rcpp_interface((void *)abstol, "N_VNew_Serial", 0)) return(1);

   for (int i = 0; i < NEQ; ++i) {
     NV_Ith_S(abstol, i) = solv_param_struc.absolute_tolerances[i];
     NV_Ith_S(y, i) = solv_param_struc.init_state[i];
   }

   reltol = solv_param_struc.reltol;

   A = SUNDenseMatrix(NEQ, NEQ);
   if (check_retval_Rcpp_interface((void *)A, "SUNDenseMatrix", 0)) return 1;
   LS = SUNLinSol_Dense(y, A);
   if (check_retval_Rcpp_interface((void*)LS, "SUNLinSol_Dense", 0)) return 1;

   // right hand side => y' = f(t,y); f(t,y) = f_E + f_I
   // f_E = explizit Part; f_I = implicit part
   // f_E is zero => fully implicit
   arkode_mem = ARKStepCreate(NULL, wrapper_ode_system_Rcpp_interface, integration_times[0], y);
   if (check_retval_Rcpp_interface((void *)arkode_mem, "ARKCreate", 0)) return(1);

   // set error handler to Rcpp::Rcerr
   //void* ptr_to_arkode_mem = &arkode_mem;
   void* ptr_to_nothing = &arkode_mem;
   retval = ARKStepSetErrHandlerFn(arkode_mem, own_error_handler_Rcpp_interface, ptr_to_nothing);

   double sum_of_least_squares = 0.;

     struct usr_data_Rcpp_interface my_ode_system = {ode_system, param_combi_start, params_time_combi_vec, params_cut_idx_vec};
     void*ptr_to_my_ode_system = &my_ode_system;
     retval = ARKStepSetUserData(arkode_mem, ptr_to_my_ode_system);
     if (check_retval_Rcpp_interface((void *)arkode_mem, "aRKStepSetUserData", 0)) return(1);

     retval = ARKStepSVtolerances(arkode_mem, reltol, abstol);
     if (check_retval_Rcpp_interface(&retval, "ARKStepSVtolerances", 1)) return(1);

     retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
     if (check_retval_Rcpp_interface(&retval, "ARKStepSetLinearSolver", 1)) return(1);

     //int CVodetmpcount=0;
     double return_time;
     float return_steps=1.;

     std::vector<double> temp_measured(init_state.size());
     //int num_states = init_state.size();

     for(int i = 0; i < NV_LENGTH_S(y); i++) {
         DF(0 , i) = NV_Ith_S(y,i);
     }

     int counter = 1;

       for (unsigned int ti = 1; ti < integration_times.size(); ti++) {
           for (int j = 1; j <= return_steps; j++) {
               return_time = integration_times[ti-1] +j/return_steps*(integration_times[ti]-integration_times[ti-1]);
               retval = ARKStepEvolve(arkode_mem, return_time, y, &t, ARK_NORMAL);
                   for (int n = 0; n < NV_LENGTH_S(y); n++) {
                     DF(counter, n) = NV_Ith_S(y,n);
                       temp_measured[n] =  hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
                       if(std::isnan(temp_measured[n])) { }
                       else {
                       sum_of_least_squares += std::abs(NV_Ith_S(y,n) - temp_measured[n]);
                       //Rcpp::Rcerr << NV_Ith_S(y,n) << "\t" << temp_measured[n] << "\t" << return_time << std::endl;
                       }
                   }
                   if (check_retval_Rcpp_interface(&retval, "CVode", 1)) {
                     break;}
           }
           counter++;
       }
       N_VDestroy(y);
       N_VDestroy(abstol);
       SUNMatDestroy(A);
       SUNLinSolFree(LS);
       ARKStepFree(&arkode_mem);

 return sum_of_least_squares/static_cast<double>(integration_times.size());
}

static int check_retval_Rcpp_interface(void *returnvalue, const char *funcname, int opt)
{
  int *retval;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    //Rcpp::Rcerr << "SUNDIALS_ERROR:" << " " << funcname << " " << "failed - returned NULL pointer" << std::endl;
    return(1); }
  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      //Rcpp::Rcerr << "SUNDIALS_ERROR:" << " " << funcname << " " << "failed with retval = " << " " << *retval << std::endl;
      return(1); }}
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    //Rcpp::Rcerr << "MEMORY_ERROR:" << " " << funcname << " " << "failed - returned NULL pointer" << std::endl;
    return(1); }

  return(0);
}
