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

// solve_ode_system_pointer
Rcpp::List solve_ode_system_pointer(std::vector<double> integration_times, Rcpp::XPtr<OS> fctptr, double relative_tolerance, std::vector<double> absolute_tolerances, Rcpp::DataFrame start, Rcpp::DataFrame states, std::string solvertype);
RcppExport SEXP _paropt_solve_ode_system_pointer(SEXP integration_timesSEXP, SEXP fctptrSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS> >::type fctptr(fctptrSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type start(startSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_ode_system_pointer(integration_times, fctptr, relative_tolerance, absolute_tolerances, start, states, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// optimizer
Rcpp::List optimizer(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, Rcpp::DataFrame lb, Rcpp::DataFrame ub, Rcpp::DataFrame states, int npop, int ngen, double error, std::string solvertype);
RcppExport SEXP _paropt_optimizer(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP statesSEXP, SEXP npopSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer(integration_times, ode_system, relative_tolerance, absolute_tolerances, lb, ub, states, npop, ngen, error, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// po
Rcpp::List po(std::vector<double> integration_times, Rcpp::XPtr<OS2> ode_sys, double relative_tolerance, std::vector<double> absolute_tolerances, Rcpp::DataFrame lower, Rcpp::DataFrame upper, Rcpp::DataFrame states, int npop, int ngen, double error, std::string solvertype);
RcppExport SEXP _paropt_po(SEXP integration_timesSEXP, SEXP ode_sysSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP statesSEXP, SEXP npopSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS2> >::type ode_sys(ode_sysSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(po(integration_times, ode_sys, relative_tolerance, absolute_tolerances, lower, upper, states, npop, ngen, error, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// so
Rcpp::List so(std::vector<double> integration_times, Rcpp::XPtr<OS2> fctptr, double relative_tolerance, std::vector<double> absolute_tolerances, Rcpp::DataFrame start, Rcpp::DataFrame states, std::string solvertype);
RcppExport SEXP _paropt_so(SEXP integration_timesSEXP, SEXP fctptrSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS2> >::type fctptr(fctptrSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type start(startSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(so(integration_times, fctptr, relative_tolerance, absolute_tolerances, start, states, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// ode_example
Rcpp::NumericVector ode_example(double t, std::vector<double> params, Rcpp::NumericVector y);
RcppExport SEXP _paropt_ode_example(SEXP tSEXP, SEXP paramsSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ode_example(t, params, y));
    return rcpp_result_gen;
END_RCPP
}
// optimizer_pointer
Rcpp::List optimizer_pointer(std::vector<double> integration_times, Rcpp::XPtr<OS> ode_sys, double relative_tolerance, std::vector<double> absolute_tolerances, Rcpp::DataFrame lower, Rcpp::DataFrame upper, Rcpp::DataFrame states, int npop, int ngen, double error, std::string solvertype);
RcppExport SEXP _paropt_optimizer_pointer(SEXP integration_timesSEXP, SEXP ode_sysSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP statesSEXP, SEXP npopSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS> >::type ode_sys(ode_sysSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_pointer(integration_times, ode_sys, relative_tolerance, absolute_tolerances, lower, upper, states, npop, ngen, error, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// solve_ode_system
Rcpp::List solve_ode_system(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, Rcpp::DataFrame start, Rcpp::DataFrame states, std::string solvertype);
RcppExport SEXP _paropt_solve_ode_system(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type start(startSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_ode_system(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, states, solvertype));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_paropt_solve_ode_system_pointer", (DL_FUNC) &_paropt_solve_ode_system_pointer, 7},
    {"_paropt_optimizer", (DL_FUNC) &_paropt_optimizer, 11},
    {"_paropt_po", (DL_FUNC) &_paropt_po, 11},
    {"_paropt_so", (DL_FUNC) &_paropt_so, 7},
    {"_paropt_ode_example", (DL_FUNC) &_paropt_ode_example, 3},
    {"_paropt_optimizer_pointer", (DL_FUNC) &_paropt_optimizer_pointer, 11},
    {"_paropt_solve_ode_system", (DL_FUNC) &_paropt_solve_ode_system, 7},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_paropt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
