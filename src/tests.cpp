/* !!revision!!
remove compiler warnings
*/

//#include "optimizer.hpp"
#include "basic_functions.hpp"
#include "params.hpp"
#include "state.hpp"
#include "param_interpolation.hpp"
#include "solver.hpp"

// [[Rcpp::export]]
int test_no_file_exist(std::string importfile) {
  bool is_file = false;
  is_file = check_file_exists(importfile);

  if(is_file == false) {
    Rcpp::stop("\nERROR: No File");
    //exit (EXIT_FAILURE);
  }
  return 0;
}


// [[Rcpp::export]]
Rcpp::NumericVector test_count_cols_rows (std::string importfile) {
    int row, col;
    count_cols_rows(importfile, row, col);
    Rcpp::NumericVector out(2);
    out[0] = row;
    out[1] = col;
    return out;
}

// [[Rcpp::export]]
std::vector<int> test_check_ncols_per_row(std::string import_path) {
  int nrow, ncol;
  count_cols_rows(import_path, nrow, ncol);
  int n_col = 0;
  std::string dummyline;
  std::ifstream myfile (import_path);

  int numCols;
  int numRows = 0;

  std::vector<int> temp;

  if (myfile.is_open())
  {
    while (getline (myfile,dummyline) )
     {
     numCols = counte_cols(dummyline);
     temp.push_back(numCols);
     numRows++;
     }
   myfile.close();
  }
  else Rcpp::Rcerr << "Unable to open file";

  for(int i = 0; i < temp.size(); i++) { // temp.size() - 1 --> last line often empty
  	if(ncol != temp[i]) {
      Rcpp::stop("\nERROR: Not the same number of columns in each row");
      //exit (EXIT_FAILURE);
  	}
  }
 return temp;
}

// [[Rcpp::export]]
std::vector<std::vector<double> > test_get_content (std::string importfile) {
    int row, col;
    count_cols_rows(importfile, row, col);
    std::vector<std::vector<double> > nc;
    get_content(importfile, nc, row, col);
    return nc;
}

// [[Rcpp::export]]
std::vector<std::string> test_get_header(std::string importfile) {
  int row, col;
  count_cols_rows(importfile, row, col);
  std::vector<std::string> test;
  get_header(importfile, test, col);
  return test;
}

// [[Rcpp::export]]
std::vector<std::vector<double> > test_remove_NA(std::string importfile) {
  int row, col;
  count_cols_rows(importfile, row, col);
  std::vector<std::vector<double> > nc;
  get_content(importfile, nc, row, col);
  std::vector<std::vector<double> > nc_without_NA;
  remove_NA(nc_without_NA, nc, row, col);
  return nc_without_NA;
}

// [[Rcpp::export]]
Rcpp::List test_Import_Parameter(std::string start, std::string lower, std::string upper) {

  std::vector<int> params_cut_idx_vec;
  std::vector<double> params_time_combi_vec;
  std::vector<double> param_combi_start;
  std::vector<double> param_combi_lb;
  std::vector<double> param_combi_ub;

  // Check if file exists
  // ===================================================
  bool is_file_start = false;
  is_file_start = check_file_exists(start);
  if(is_file_start == false) {
    Rcpp::stop("\nERROR: No File with start values");
    //exit (EXIT_FAILURE);
  }

  bool is_file_lb = false;
  is_file_lb = check_file_exists(lower);
  if(is_file_lb == false) {
    Rcpp::stop("\nERROR: No File with lowerbound values");
    //exit (EXIT_FAILURE);
  }

  bool is_file_ub = false;
  is_file_ub = check_file_exists(upper);
  if(is_file_ub == false) {
    Rcpp::stop("\nERROR: No File with upperbound values");
    //exit (EXIT_FAILURE);
  }
  // ===================================================

  // Import start, lower and upper
  // ===================================================
  int ncol, nrow;
  int temp_nrow;
  count_cols_rows(start, nrow, ncol);
  check_ncols_per_row(start, temp_nrow, ncol);
  std::vector<std::vector<double> > nc;
  get_content(start, nc, nrow, ncol);
  std::vector<std::string> test;
  get_header(start, test, ncol);
  std::vector<std::vector<double> > nc_without_NA;
  remove_NA(nc_without_NA, nc, nrow, ncol);

  // lower_bounds
  int ncol_lb, nrow_lb;
  int temp_nrow_lb;
  count_cols_rows(lower, nrow_lb, ncol_lb);
  check_ncols_per_row(lower, temp_nrow_lb, ncol_lb);
  std::vector<std::vector<double> > nc_lb;
  get_content(lower, nc_lb, nrow_lb, ncol_lb);
  std::vector<std::string> test_lb;
  get_header(lower, test_lb, ncol_lb);
  std::vector<std::vector<double> > nc_without_NA_lb;
  remove_NA(nc_without_NA_lb, nc_lb, nrow_lb, ncol_lb);

  // upper_bounds
  int ncol_ub, nrow_ub;
  int temp_nrow_ub;
  count_cols_rows(upper, nrow_ub, ncol_ub);
  check_ncols_per_row(upper, temp_nrow_ub, ncol_ub);
  std::vector<std::vector<double> > nc_ub;
  get_content(upper, nc_ub, nrow_ub, ncol_ub);
  std::vector<std::string> test_ub;
  get_header(upper, test_ub, ncol_ub);
  std::vector<std::vector<double> > nc_without_NA_ub;
  remove_NA(nc_without_NA_ub, nc_ub, nrow_ub, ncol_ub);
  // ===================================================

  // Error checks
  // ===================================================
  if(test.size() != test_lb.size() || test.size() != test_ub.size()) {
    Rcpp::stop("\nERROR: Different number of names of parameters between startvalues, lower bounds and upper bounds");
    //exit (EXIT_FAILURE);
  }

  for(size_t i = 0; i < test.size(); i++) {
    if(test[i] != test_lb[i] || test[i] != test_ub[i]) {
      Rcpp::stop("\nERROR: Different of names for parameters between startvalues, lower bounds and upper bounds");
      //exit (EXIT_FAILURE);
    }
  }

  if(ncol != ncol_lb || ncol != ncol_ub) {
    Rcpp::stop("\nERROR: Different number of cols between startvalues, lower bounds and upper bounds");
    //exit (EXIT_FAILURE);
  }
  if(nrow != nrow_lb || nrow != nrow_ub) {
    Rcpp::stop("\nERROR: Different number of rows between startvalues, lower bounds and upper bounds");
    //exit (EXIT_FAILURE);
  }

  if(nc.size() != nc_lb.size() || nc.size() != nc_ub.size()) {
    Rcpp::stop("\nERROR: Different number of cols of imported numeric values between startvalues, lower bounds and upper bounds");
    //exit (EXIT_FAILURE);
  }

  for(size_t i = 0; i < ncol; i++) {
      if(nc[i].size() != nc_lb[i].size() || nc[i].size() != nc_ub[i].size()) {
        Rcpp::Rcerr << "number of rows start values:" << "\t" << nc[i].size() << std::endl;
        Rcpp::Rcerr << "number of rows lb values:" << "\t" << nc_lb[i].size() << std::endl;
        Rcpp::Rcerr << "number of rows ub values:" << "\t" << nc_ub[i].size() << std::endl;
        Rcpp::stop("\nERROR: Different number of rows of imported numeric values between startvalues, lower bounds and upper bounds");
        //exit (EXIT_FAILURE);
      }
    }

  if(nc_without_NA.size() != nc_without_NA_lb.size() || nc_without_NA.size() != nc_without_NA_ub.size()) {
    Rcpp::stop("\nERROR: Different number of cols of imported numeric values without NA between startvalues, lower bounds and upper bounds");
    //exit (EXIT_FAILURE);
  }
  // ===================================================

  // extract time vector and error check
  // ===================================================
  bool time_check = false;
  time_check = check_time_column(nc[0], "time");
  if(time_check == false) {
    Rcpp::stop("\nERROR: time column of start values not correct");
    //exit (EXIT_FAILURE);
  }
  std::vector<double> time(nc[0].size());

  bool time_check_lb = false;
  time_check_lb = check_time_column(nc_lb[0], "time");
  if(time_check_lb == false) {
    Rcpp::stop("\nERROR: time column of lower bounds not correct");
    //exit (EXIT_FAILURE);
  }

  bool time_check_ub = false;
  time_check_ub = check_time_column(nc_ub[0], "time");
  if(time_check_ub == false) {
    Rcpp::stop("\nERROR: time column of upper bounds not correct");
    //exit (EXIT_FAILURE);
  }

  for(size_t i = 0; i < time.size(); i++) {
    if(nc[0][i] != nc_lb[0][i] || nc[0][i] != nc_ub[0][i]) {
      Rcpp::stop("\nERROR:Different time vector for startvalues, lower bounds and upper bounds");
      //exit (EXIT_FAILURE);
    }
    time[i] = nc[0][i];
  }
  // ===================================================

  // extract number of rows without NA and timepoints
  // ===================================================
  std::vector<int> number_rows_without_NA(ncol);
  std::vector<int> number_rows_without_NA_lb(ncol_lb);
  std::vector<int> number_rows_without_NA_ub(ncol_ub);
  std::vector<std::vector<int> > time_points(ncol);
  std::vector<std::vector<int> > time_points_lb(ncol_lb);
  std::vector<std::vector<int> > time_points_ub(ncol_ub);
  for(size_t i = 0; i < ncol; i++) { // i = 1
    for(size_t j = 0; j < (nrow-1); j++) {
      if(std::isnan(nc[i][j])) { }
      else {
        number_rows_without_NA[i] = number_rows_without_NA[i] + 1;
        time_points[i].push_back(j);
      }

      if(std::isnan(nc_lb[i][j])) { }
      else {
        number_rows_without_NA_lb[i] = number_rows_without_NA_lb[i] + 1;
        time_points_lb[i].push_back(j);
      }

      if(std::isnan(nc_ub[i][j])) { }
      else {
        number_rows_without_NA_ub[i] = number_rows_without_NA_ub[i] + 1;
        time_points_ub[i].push_back(j);
      }
    }
  }
  // ===================================================

  // Error checks for: number of rows without NA and timepoints
  // ===================================================
  if(time_points.size() != time_points_lb.size() || time_points.size() != time_points_ub.size()) {
    Rcpp::stop("\nERROR: Different number of time_points between startvalues, lower bounds and upper bounds: check if NAs are at the same position in all files");
    //exit (EXIT_FAILURE);
  }

  for(size_t i = 0; i < time_points.size(); i++) {
    if(time_points[i].size() != time_points_lb[i].size() || time_points[i].size() != time_points_ub[i].size()) {
      Rcpp::stop("\nERROR: Different number of time_points between startvalues, lower bounds and upper bounds");
      //exit (EXIT_FAILURE);
    }
    for(size_t j = 0; j < time_points[i].size(); j++) {
      if(time_points[i][j] != time_points_lb[i][j] || time_points[i][j] != time_points_ub[i][j]) {
        Rcpp::stop("\nERROR: Different entries in time_points between startvalues, lower bounds and upper bounds");
        //exit (EXIT_FAILURE);
      }
    }
  }

  for(size_t i = 0; i < number_rows_without_NA.size(); i++) {
    if(number_rows_without_NA[i] != number_rows_without_NA_lb[i] || number_rows_without_NA[i] != number_rows_without_NA_ub[i]) {
      Rcpp::stop("\nERROR: NA values not at the same position in startvalues, lower bounds and upper bounds");
      //exit (EXIT_FAILURE);
    }
  }

  // ===================================================

  // define combination of parameter, lower_bounds, upper_bounds and time and store it in std::list
  // ===================================================
  std::list<ParamClass> paramlist;
  for(size_t i = 1; i < ncol; i++) {
    std::vector<double> temp_param(number_rows_without_NA[i]);
    std::vector<double> temp_lb(number_rows_without_NA[i]);
    std::vector<double> temp_ub(number_rows_without_NA[i]);
    std::vector<double> temp_time(number_rows_without_NA[i]);
    for(size_t j = 0; j < number_rows_without_NA[i]; j++) {
      temp_param[j] = nc_without_NA[i][j];
      temp_lb[j] = nc_without_NA_lb[i][j];
      temp_ub[j] = nc_without_NA_ub[i][j];
      temp_time[j] = nc[0][time_points[i][j]];
    }
    if(temp_time.size() == 1) {
      Rcpp::Rcerr << test[i] << ":" << "\t" << "Parameter is const" << std::endl;
      ParamClass temp(temp_param[0], temp_lb[0], temp_ub[0]);
      paramlist.push_back(temp);
    } else {
      Rcpp::Rcerr << test[i] << ":" << "\t" << "Parameter is variabel" << std::endl;
      ParamClass temp(temp_time.size(), temp_time, temp_param, temp_lb, temp_ub);
      paramlist.push_back(temp);
    }
  }

  ParamOrderClass ParamOrder(paramlist);
  ParamOrder.cut_idx(params_cut_idx_vec);
  ParamOrder.get_time_combi(params_time_combi_vec);
  ParamOrder.get_param_combi(param_combi_start);
  ParamOrder.get_lb_combi(param_combi_lb);
  ParamOrder.get_up_combi(param_combi_ub);

  return Rcpp::List::create(Rcpp::Named("cut") = params_cut_idx_vec,
                     Rcpp::Named("time") = params_time_combi_vec,
                     Rcpp::Named("start") = param_combi_start,
                     Rcpp::Named("lower") = param_combi_lb,
                     Rcpp::Named("upper") = param_combi_ub);
}


// [[Rcpp::export]]
Rcpp::List test_Import_States(std::string start) {

std::vector<int> hs_cut_idx_vec;
std::vector<double> hs_time_combi_vec;
std::vector<double> hs_harvest_state_combi_vec;

  bool is_file_start = false;
  is_file_start = check_file_exists(start);
  if(is_file_start == false) {
    Rcpp::stop("\nERROR: No File with state values");
    //exit (EXIT_FAILURE);
  }

  int ncol, nrow;
  int temp_nrow;
  count_cols_rows(start, nrow, ncol);
  check_ncols_per_row(start, temp_nrow, ncol);
  std::vector<std::vector<double> > nc;
  get_content(start, nc, nrow, ncol);
  std::vector<std::string> test;
  get_header(start, test, ncol);

  // extract time vector
  bool time_check = false;
  time_check = check_time_column(nc[0], "time");
  std::vector<double> time(nc[0].size());

  for(size_t i = 0; i < time.size(); i++) {
    time[i] = nc[0][i];
  }

  int size_of_states = time.size();
  // define combination of state and time and store it in std::list
  std::list<HarvestStateClass> paramlist;
  for(size_t i = 1; i < ncol; i++) {
    std::vector<double> temp_state(size_of_states);
    std::vector<double> temp_time(size_of_states);
    for(size_t j = 0; j < size_of_states; j++) {
      temp_state[j] = nc[i][j];
      temp_time[j] = nc[0][j];
    }
      HarvestStateClass temp(temp_time, temp_state);
      paramlist.push_back(temp);
  }

  HarvestStateOrderClass HarvestOrder(paramlist);
  HarvestOrder.cut_idx(hs_cut_idx_vec);
  HarvestOrder.get_harvest_time_combi(hs_time_combi_vec);
  HarvestOrder.get_harvest_state_combi(hs_harvest_state_combi_vec);

  return Rcpp::List::create(Rcpp::Named("cut") = hs_cut_idx_vec,
                     Rcpp::Named("time") = hs_time_combi_vec,
                     Rcpp::Named("start") = hs_harvest_state_combi_vec);
}

// [[Rcpp::export]]
double test_error_calculation(std::string measured, std::string solver_output) {

  std::vector<int> hs_cut_idx_vec;
  std::vector<double> hs_time_combi_vec;
  std::vector<double> hs_harvest_state_combi_vec;
  std::vector<std::string> header_states;

  Import_states(measured, hs_cut_idx_vec, hs_time_combi_vec, hs_harvest_state_combi_vec, header_states);

  int row, col;
  count_cols_rows(solver_output, row, col);
  std::vector<std::vector<double> > nc;
  get_content(solver_output, nc, row, col);

  double error = 0;
  std::vector<double> integration_times = {1., 2., 3., 4., 5.};
  std::vector<double> temp_measured(integration_times.size());
  for(int ti = 1; ti < integration_times.size(); ti++) {
    Rcpp::Rcerr << "time:" << "\t" << integration_times[ti] << std::endl;
    for(size_t n = 0; n < hs_cut_idx_vec.size(); n++) {
        temp_measured[n] = hs_harvest_state_combi_vec[hs_cut_idx_vec[n] * n + ti];
        error += std::abs(temp_measured[n] - nc[n+1][ti]);
        Rcpp::Rcerr << temp_measured[n] << "\t" << nc[n+1][ti] << std::endl;
    }
  }

  return error;
}

// [[Rcpp::export]]
int test_interface_fct(Rcpp::NumericVector integration_times,
                   SEXP ode_system, double relative_tolerance,
                   Rcpp::NumericVector absolute_tolerances, std::string start,
                     std::string lower, std::string upper, std::string states,
                   int npop, int ngen, double error, std::string where_to_save_output_states,
                 std::string where_to_save_output_parameter) {

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

if(init_state.size() > absolute_tolerances.length()) {
  Rcpp::stop("\nERROR: absolute tolerances not defined for each state");
  //exit (EXIT_FAILURE);
}

if(init_state.size() < absolute_tolerances.length()) {
  Rcpp::stop("\nERROR: dimension error for absolute tolerances");
  //exit (EXIT_FAILURE);
}

std::vector<double>::iterator max_time_param_vector = std::max_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
std::vector<double>::iterator min_time_param_vector = std::min_element(params_time_combi_vec.begin(), params_time_combi_vec.end());
std::vector<double>::iterator max_time_harvest_vector = std::max_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());
std::vector<double>::iterator min_time_harvest_vector = std::min_element(hs_time_combi_vec.begin(), hs_time_combi_vec.end());

/*
if(*max_time_param_vector > *max_time_harvest_vector) {
  Rcpp::stop("\nERROR: Maximum of timevector of parameter larger then corresponding timepoint of state vector");
  //exit (EXIT_FAILURE);
} else if(*max_time_param_vector < *max_time_harvest_vector) {
  Rcpp::stop("\nERROR: Maximum of timevector of parameter smaller then corresponding timepoint of state vector");
  //exit (EXIT_FAILURE);
}
if(*min_time_param_vector > *min_time_harvest_vector) {
  Rcpp::stop("\nERROR: Minimum of timevector of parameter larger then corresponding timepoint of state vector");
  //exit (EXIT_FAILURE);
} else if(*min_time_param_vector < *min_time_harvest_vector) {
  Rcpp::stop("\nERROR: Minimum of timevector of parameter smaller then corresponding timepoint of state vector");
  //exit (EXIT_FAILURE);
}
*/
bool max_time_diff_zero = double_diff(*max_time_param_vector, 0);
bool max_time_param_vs_max_time_harvest = double_diff(*max_time_param_vector, *max_time_harvest_vector);
bool min_time_param_vs_min_time_harvest = double_diff(*min_time_param_vector, *min_time_harvest_vector);
if(!max_time_diff_zero) {
  if(*max_time_param_vector > *max_time_harvest_vector) { //check || or &&
    Rcpp::stop("\nERROR: Maximum of timevector of parameter larger then corresponding timepoint of state vector");
  } else if(*max_time_param_vector < *max_time_harvest_vector) {
    Rcpp::stop("\nERROR: Maximum of timevector of parameter smaller then corresponding timepoint of state vector");
  }
  if(*min_time_param_vector > *min_time_harvest_vector) {
    Rcpp::stop("\nERROR: Minimum of timevector of parameter larger then corresponding timepoint of state vector");
  } else if(*min_time_param_vector < *min_time_harvest_vector) {
    Rcpp::stop("\nERROR: Minimum of timevector of parameter smaller then corresponding timepoint of state vector");
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
  //exit (EXIT_FAILURE);
}

bool check_entries_time;
for(size_t i = 0; i < integration_times.length(); i++) {
  check_entries_time = double_diff(integration_times[i],integration_time_assumption[i]);
  if(!check_entries_time) {
    Rcpp::stop("\nERROR: integration_times has not the same entries as the time vector of state input");
    //exit (EXIT_FAILURE);
  }
}

for(size_t i = 0; i < params_cut_idx_vec.size(); i++) {
  if(params_cut_idx_vec[i] == 1 || params_cut_idx_vec[i] >=4) {
    // everything is fine. 4 values needed for spline
  } else{
    Rcpp::stop("\nERROR: neither constant nor variable parameter. Variable parameters need at least four datapoints!");
  }
}

if(TYPEOF(ode_system) != CLOSXP) {
  Rcpp::stop("\nERROR: type of odesystem should be closure");
  //exit (EXIT_FAILURE);
}

Rcpp::NumericVector new_states(init_state.size());
realtype time_param_sort = params_time_combi_vec[0];
std::vector<double> parameter_input;
double time = params_time_combi_vec[0];

params_sort(time_param_sort, parameter_input, param_combi_start, params_time_combi_vec, params_cut_idx_vec);

Rcpp::NumericVector current_states(hs_cut_idx_vec.size());
for (size_t i = 0; i < hs_cut_idx_vec.size(); i++) {
  current_states[i] = init_state[i];
}

Rcpp::Function odes = ode_system;
Rcpp::Function fo("formals");
Rcpp::List test_num_arguments = fo(odes);
if(test_num_arguments.length() != 3) {
  Rcpp::stop("\nERROR: odesystem should only accept three arguments");
  //exit (EXIT_FAILURE);
}

try {odes(time, parameter_input, current_states); } catch (...) {
Rcpp::stop("\nERROR: odesystem cannot be called. May be wrong types of arguments (double, std::vector<double>, Rcpp::NumericVector)?");
//exit (EXIT_FAILURE);
}

bool output_numeric = false;
switch(TYPEOF(odes(time, parameter_input, current_states))) {
  case REALSXP: {
    Rcpp::Rcerr << "correct output" << std::endl;
    output_numeric = true;
  }
}

if(output_numeric == true) {
  try {new_states = odes(time, parameter_input, current_states); } catch (...) {
    Rcpp::stop("\nERROR: output of odesystem is wrong! Has to be Rcpp::NumericVector");
    //exit (EXIT_FAILURE);
  }
} else {
  Rcpp::stop("\nERROR: output of odesystem is wrong! Has to be Rcpp::NumericVector");
}

if(new_states.length() != init_state.size()) {
  Rcpp::stop("\nERROR: output of odesystem is wrong! Has to be same size as number of states");
  //exit (EXIT_FAILURE);
}

return 0;
}
