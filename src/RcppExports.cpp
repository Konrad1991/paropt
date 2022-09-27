// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/paropt_types.h"
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// wrapper_optimizer
Rcpp::List wrapper_optimizer(vd& init_state, vd& par_times, vi& param_idx_cuts, vd& lb_, vd& ub_, vd& state_measured, vi& state_idx_cuts, vd& integration_times, double reltol, vd& absolute_tolerances, Rcpp::XPtr<OS> fct, int nswarm, int ngen, double error, int solvertype);
RcppExport SEXP _paropt_wrapper_optimizer(SEXP init_stateSEXP, SEXP par_timesSEXP, SEXP param_idx_cutsSEXP, SEXP lb_SEXP, SEXP ub_SEXP, SEXP state_measuredSEXP, SEXP state_idx_cutsSEXP, SEXP integration_timesSEXP, SEXP reltolSEXP, SEXP absolute_tolerancesSEXP, SEXP fctSEXP, SEXP nswarmSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vd& >::type init_state(init_stateSEXP);
    Rcpp::traits::input_parameter< vd& >::type par_times(par_timesSEXP);
    Rcpp::traits::input_parameter< vi& >::type param_idx_cuts(param_idx_cutsSEXP);
    Rcpp::traits::input_parameter< vd& >::type lb_(lb_SEXP);
    Rcpp::traits::input_parameter< vd& >::type ub_(ub_SEXP);
    Rcpp::traits::input_parameter< vd& >::type state_measured(state_measuredSEXP);
    Rcpp::traits::input_parameter< vi& >::type state_idx_cuts(state_idx_cutsSEXP);
    Rcpp::traits::input_parameter< vd& >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< vd& >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS> >::type fct(fctSEXP);
    Rcpp::traits::input_parameter< int >::type nswarm(nswarmSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< int >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(wrapper_optimizer(init_state, par_times, param_idx_cuts, lb_, ub_, state_measured, state_idx_cuts, integration_times, reltol, absolute_tolerances, fct, nswarm, ngen, error, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// wrapper_solver
Rcpp::List wrapper_solver(vd& init_state, vd& par_times, vi& param_idx_cuts, vd& parameter_vec, vd& state_measured, vi& state_idx_cuts, vd& integration_times, double reltol, vd& absolute_tolerances, Rcpp::XPtr<OS> fct, int solvertype);
RcppExport SEXP _paropt_wrapper_solver(SEXP init_stateSEXP, SEXP par_timesSEXP, SEXP param_idx_cutsSEXP, SEXP parameter_vecSEXP, SEXP state_measuredSEXP, SEXP state_idx_cutsSEXP, SEXP integration_timesSEXP, SEXP reltolSEXP, SEXP absolute_tolerancesSEXP, SEXP fctSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vd& >::type init_state(init_stateSEXP);
    Rcpp::traits::input_parameter< vd& >::type par_times(par_timesSEXP);
    Rcpp::traits::input_parameter< vi& >::type param_idx_cuts(param_idx_cutsSEXP);
    Rcpp::traits::input_parameter< vd& >::type parameter_vec(parameter_vecSEXP);
    Rcpp::traits::input_parameter< vd& >::type state_measured(state_measuredSEXP);
    Rcpp::traits::input_parameter< vi& >::type state_idx_cuts(state_idx_cutsSEXP);
    Rcpp::traits::input_parameter< vd& >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< vd& >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS> >::type fct(fctSEXP);
    Rcpp::traits::input_parameter< int >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(wrapper_solver(init_state, par_times, param_idx_cuts, parameter_vec, state_measured, state_idx_cuts, integration_times, reltol, absolute_tolerances, fct, solvertype));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_paropt_wrapper_optimizer", (DL_FUNC) &_paropt_wrapper_optimizer, 15},
    {"_paropt_wrapper_solver", (DL_FUNC) &_paropt_wrapper_solver, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_paropt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
