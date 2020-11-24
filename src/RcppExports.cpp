// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/paropt.h"
#include "../inst/include/paropt_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// function_access
Rcpp::List function_access(std::vector<double> integration_times, Rcpp::XPtr<OS> fctptr, double relative_tolerance, std::vector<double> absolute_tolerances, std::string start, std::string states, std::string where_to_save_output_states, std::string solvertype);
RcppExport SEXP _paropt_function_access(SEXP integration_timesSEXP, SEXP fctptrSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP statesSEXP, SEXP where_to_save_output_statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS> >::type fctptr(fctptrSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_states(where_to_save_output_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(function_access(integration_times, fctptr, relative_tolerance, absolute_tolerances, start, states, where_to_save_output_states, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// interface_function
Rcpp::List interface_function(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, std::string start, std::string lower, std::string upper, std::string states, int npop, int ngen, double error, std::string where_to_save_output_states, std::string where_to_save_output_parameter, std::string solvertype);
RcppExport SEXP _paropt_interface_function(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP statesSEXP, SEXP npopSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP where_to_save_output_statesSEXP, SEXP where_to_save_output_parameterSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< std::string >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_states(where_to_save_output_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_parameter(where_to_save_output_parameterSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(interface_function(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, lower, upper, states, npop, ngen, error, where_to_save_output_states, where_to_save_output_parameter, solvertype));
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
// optimizer_access_in_Rcpp
Rcpp::List optimizer_access_in_Rcpp(std::vector<double> integration_times, Rcpp::XPtr<OS> ode_sys, double relative_tolerance, std::vector<double> absolute_tolerances, std::string start, std::string lower, std::string upper, std::string states, int npop, int ngen, double error, std::string where_to_save_output_states, std::string where_to_save_output_parameter, std::string solvertype);
RcppExport SEXP _paropt_optimizer_access_in_Rcpp(SEXP integration_timesSEXP, SEXP ode_sysSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP statesSEXP, SEXP npopSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP where_to_save_output_statesSEXP, SEXP where_to_save_output_parameterSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<OS> >::type ode_sys(ode_sysSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< std::string >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_states(where_to_save_output_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_parameter(where_to_save_output_parameterSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_access_in_Rcpp(integration_times, ode_sys, relative_tolerance, absolute_tolerances, start, lower, upper, states, npop, ngen, error, where_to_save_output_states, where_to_save_output_parameter, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// solve_ode_system
Rcpp::List solve_ode_system(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, std::string start, std::string states, std::string where_to_save_output_states, std::string solvertype);
RcppExport SEXP _paropt_solve_ode_system(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP statesSEXP, SEXP where_to_save_output_statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_states(where_to_save_output_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_ode_system(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, states, where_to_save_output_states, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// test_paramsort_and_spline
std::vector<std::vector<double> > test_paramsort_and_spline(Rcpp::NumericVector timepoints, std::string start, std::string lower, std::string upper);
RcppExport SEXP _paropt_test_paramsort_and_spline(SEXP timepointsSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type timepoints(timepointsSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< std::string >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(test_paramsort_and_spline(timepoints, start, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// test_solver
std::vector<std::vector<double> > test_solver(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, std::string start, std::string lower, std::string upper, std::string states, std::string solvertype);
RcppExport SEXP _paropt_test_solver(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< std::string >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(test_solver(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, lower, upper, states, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// test_solve_ode_system
std::vector<std::vector<double> > test_solve_ode_system(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, std::string start, std::string states, std::string solvertype);
RcppExport SEXP _paropt_test_solve_ode_system(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP statesSEXP, SEXP solvertypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type solvertype(solvertypeSEXP);
    rcpp_result_gen = Rcpp::wrap(test_solve_ode_system(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, states, solvertype));
    return rcpp_result_gen;
END_RCPP
}
// test_no_file_exist
int test_no_file_exist(std::string importfile);
RcppExport SEXP _paropt_test_no_file_exist(SEXP importfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type importfile(importfileSEXP);
    rcpp_result_gen = Rcpp::wrap(test_no_file_exist(importfile));
    return rcpp_result_gen;
END_RCPP
}
// test_count_cols_rows
Rcpp::NumericVector test_count_cols_rows(std::string importfile);
RcppExport SEXP _paropt_test_count_cols_rows(SEXP importfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type importfile(importfileSEXP);
    rcpp_result_gen = Rcpp::wrap(test_count_cols_rows(importfile));
    return rcpp_result_gen;
END_RCPP
}
// test_check_ncols_per_row
std::vector<int> test_check_ncols_per_row(std::string import_path);
RcppExport SEXP _paropt_test_check_ncols_per_row(SEXP import_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type import_path(import_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(test_check_ncols_per_row(import_path));
    return rcpp_result_gen;
END_RCPP
}
// test_get_content
std::vector<std::vector<double> > test_get_content(std::string importfile);
RcppExport SEXP _paropt_test_get_content(SEXP importfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type importfile(importfileSEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_content(importfile));
    return rcpp_result_gen;
END_RCPP
}
// test_get_header
std::vector<std::string> test_get_header(std::string importfile);
RcppExport SEXP _paropt_test_get_header(SEXP importfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type importfile(importfileSEXP);
    rcpp_result_gen = Rcpp::wrap(test_get_header(importfile));
    return rcpp_result_gen;
END_RCPP
}
// test_remove_NA
std::vector<std::vector<double> > test_remove_NA(std::string importfile);
RcppExport SEXP _paropt_test_remove_NA(SEXP importfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type importfile(importfileSEXP);
    rcpp_result_gen = Rcpp::wrap(test_remove_NA(importfile));
    return rcpp_result_gen;
END_RCPP
}
// test_Import_Parameter
Rcpp::List test_Import_Parameter(std::string start, std::string lower, std::string upper);
RcppExport SEXP _paropt_test_Import_Parameter(SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< std::string >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Import_Parameter(start, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// test_Import_States
Rcpp::List test_Import_States(std::string start);
RcppExport SEXP _paropt_test_Import_States(SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Import_States(start));
    return rcpp_result_gen;
END_RCPP
}
// test_error_calculation
double test_error_calculation(std::string measured, std::string solver_output);
RcppExport SEXP _paropt_test_error_calculation(SEXP measuredSEXP, SEXP solver_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type measured(measuredSEXP);
    Rcpp::traits::input_parameter< std::string >::type solver_output(solver_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(test_error_calculation(measured, solver_output));
    return rcpp_result_gen;
END_RCPP
}
// test_interface_fct
int test_interface_fct(Rcpp::NumericVector integration_times, SEXP ode_system, double relative_tolerance, Rcpp::NumericVector absolute_tolerances, std::string start, std::string lower, std::string upper, std::string states, int npop, int ngen, double error, std::string where_to_save_output_states, std::string where_to_save_output_parameter);
RcppExport SEXP _paropt_test_interface_fct(SEXP integration_timesSEXP, SEXP ode_systemSEXP, SEXP relative_toleranceSEXP, SEXP absolute_tolerancesSEXP, SEXP startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP statesSEXP, SEXP npopSEXP, SEXP ngenSEXP, SEXP errorSEXP, SEXP where_to_save_output_statesSEXP, SEXP where_to_save_output_parameterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type integration_times(integration_timesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ode_system(ode_systemSEXP);
    Rcpp::traits::input_parameter< double >::type relative_tolerance(relative_toleranceSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type absolute_tolerances(absolute_tolerancesSEXP);
    Rcpp::traits::input_parameter< std::string >::type start(startSEXP);
    Rcpp::traits::input_parameter< std::string >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< std::string >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< std::string >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_states(where_to_save_output_statesSEXP);
    Rcpp::traits::input_parameter< std::string >::type where_to_save_output_parameter(where_to_save_output_parameterSEXP);
    rcpp_result_gen = Rcpp::wrap(test_interface_fct(integration_times, ode_system, relative_tolerance, absolute_tolerances, start, lower, upper, states, npop, ngen, error, where_to_save_output_states, where_to_save_output_parameter));
    return rcpp_result_gen;
END_RCPP
}
// test_Import_Parameter_DF
Rcpp::List test_Import_Parameter_DF(Rcpp::DataFrame lb, Rcpp::DataFrame ub);
RcppExport SEXP _paropt_test_Import_Parameter_DF(SEXP lbSEXP, SEXP ubSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ub(ubSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Import_Parameter_DF(lb, ub));
    return rcpp_result_gen;
END_RCPP
}
// test_Import_Start_Parameter_DF
Rcpp::List test_Import_Start_Parameter_DF(Rcpp::DataFrame start);
RcppExport SEXP _paropt_test_Import_Start_Parameter_DF(SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Import_Start_Parameter_DF(start));
    return rcpp_result_gen;
END_RCPP
}
// test_Import_States_DF
Rcpp::List test_Import_States_DF(Rcpp::DataFrame states);
RcppExport SEXP _paropt_test_Import_States_DF(SEXP statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type states(statesSEXP);
    rcpp_result_gen = Rcpp::wrap(test_Import_States_DF(states));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _paropt_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _paropt_RcppExport_registerCCallable() { 
    R_RegisterCCallable("paropt", "_paropt_RcppExport_validate", (DL_FUNC)_paropt_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_paropt_function_access", (DL_FUNC) &_paropt_function_access, 8},
    {"_paropt_interface_function", (DL_FUNC) &_paropt_interface_function, 14},
    {"_paropt_ode_example", (DL_FUNC) &_paropt_ode_example, 3},
    {"_paropt_optimizer_access_in_Rcpp", (DL_FUNC) &_paropt_optimizer_access_in_Rcpp, 14},
    {"_paropt_solve_ode_system", (DL_FUNC) &_paropt_solve_ode_system, 8},
    {"_paropt_test_paramsort_and_spline", (DL_FUNC) &_paropt_test_paramsort_and_spline, 4},
    {"_paropt_test_solver", (DL_FUNC) &_paropt_test_solver, 9},
    {"_paropt_test_solve_ode_system", (DL_FUNC) &_paropt_test_solve_ode_system, 7},
    {"_paropt_test_no_file_exist", (DL_FUNC) &_paropt_test_no_file_exist, 1},
    {"_paropt_test_count_cols_rows", (DL_FUNC) &_paropt_test_count_cols_rows, 1},
    {"_paropt_test_check_ncols_per_row", (DL_FUNC) &_paropt_test_check_ncols_per_row, 1},
    {"_paropt_test_get_content", (DL_FUNC) &_paropt_test_get_content, 1},
    {"_paropt_test_get_header", (DL_FUNC) &_paropt_test_get_header, 1},
    {"_paropt_test_remove_NA", (DL_FUNC) &_paropt_test_remove_NA, 1},
    {"_paropt_test_Import_Parameter", (DL_FUNC) &_paropt_test_Import_Parameter, 3},
    {"_paropt_test_Import_States", (DL_FUNC) &_paropt_test_Import_States, 1},
    {"_paropt_test_error_calculation", (DL_FUNC) &_paropt_test_error_calculation, 2},
    {"_paropt_test_interface_fct", (DL_FUNC) &_paropt_test_interface_fct, 13},
    {"_paropt_test_Import_Parameter_DF", (DL_FUNC) &_paropt_test_Import_Parameter_DF, 2},
    {"_paropt_test_Import_Start_Parameter_DF", (DL_FUNC) &_paropt_test_Import_Start_Parameter_DF, 1},
    {"_paropt_test_Import_States_DF", (DL_FUNC) &_paropt_test_Import_States_DF, 1},
    {"_paropt_RcppExport_registerCCallable", (DL_FUNC) &_paropt_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_paropt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
