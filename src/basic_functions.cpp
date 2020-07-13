//#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <vector>
#include <cctype>
#include <cmath>
#include <stdio.h>
#include "state.hpp"
#include "params.hpp"

#include <Rcpp.h>

#define NA std::nan("l")

// check wether file exists
bool check_file_exists(std::string import_path) {
  std::ifstream myfile (import_path);
  bool rt = false;
  if(myfile.is_open()) {
    rt = true;
  }
  return rt;
}

// count columns of text file
int counte_cols(std::string line) {
std::istringstream iss(line);
int columns = 0;
do
{
    std::string sub;
    iss >> sub;
    if (sub.length())
        ++columns;
}
while(iss);

return columns;
}

// count rows and columns of text file
void count_cols_rows(std::string import_path, int &nrow, int &ncol) {

//int n_row = 0;
//int n_col = 0;
std::string dummyline;
std::ifstream myfile (import_path);

int numCols = 0; // !
int numRows = 0;

if (myfile.is_open())
{
  while (getline (myfile,dummyline) )
   {
   numCols = counte_cols(dummyline);
   numRows++;
   }
 myfile.close();
}
else Rcpp::Rcerr << "Unable to open file";
myfile.close();
nrow = numRows;
ncol = numCols;

}

//Check if each row contains the same number of columns
void check_ncols_per_row(std::string import_path, int &nrow, int &ncol) {

//int n_col = 0;
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

for(size_t i = 0; i < temp.size(); i++) { // temp.size() - 1 --> last line often empty
	if(ncol != temp[i]) {
    Rcpp::stop("\nERROR: Not the same number of columns in each row");
    //exit (EXIT_FAILURE);
	}
}
}

// used for conversion from string to double
template<class T> T fromString(const std::string& s)
{
    std::size_t dot, min_idx, e_idx, E_idx, plus_idx;
    dot = s.find('.');
    min_idx = s.find('-');
    e_idx = s.find('e');
    E_idx = s.find('E');
    plus_idx = s.find('+');
     for(std::string::size_type i = 0; i < s.size(); i++) {
       if(std::isdigit(s[i]) || i == dot || i == min_idx || i == e_idx || i == E_idx || i == plus_idx) {
         // it is a digit or "."
       } else {
         Rcpp::stop("\nERROR: Not a number in input. As non-numeric input only NA allowed");
       }
     }

     std::istringstream stream (s);
     T t;
     stream >> t;
     return t;
}


// import matrix (not the first line = headers) into std::vector<std::vector<double>>
void get_content(std::string import_path, std::vector<std::vector<double> > &numeric_content, int n_row, int n_col) {

	numeric_content.resize(n_col); //n_row
	for(size_t i = 0; i < numeric_content.size(); i++) {
		numeric_content[i].resize((n_row-1));
	}

std::vector<double> temp(n_col*(n_row-1));

std::string dummyline;
std::ifstream myfile (import_path);

size_t line_counter = 0;

if (myfile.is_open())
{
  size_t counter = 0;
  while (getline (myfile,dummyline) )
   {

   if(line_counter != 0) {
   std::istringstream iss(dummyline);

   	do
   	{
   		std::string sub;
   		iss >> sub;
   		if(sub.length()) {
   			if(sub == "NA") {
   				temp[counter] = NA;
   				counter++;
   			} else {
   			    temp[counter] = fromString<double>(sub);
   		        counter++; }
   		}
   	}
   	while(iss);
   }

   line_counter++;
   }
 myfile.close();
}
else Rcpp::Rcerr << "Unable to open file";

int size_counter = 0;

for(int i = 0; i < (n_row - 1); i++) {
	for(int j = 0; j < n_col; j++) {
		numeric_content[j][i] = temp[size_counter];
		size_counter++;
	}
}

}

// import headers into std::vector<std::string>
void get_header (std::string import_path, std::vector<std::string> &header, int n_col)
{
header.resize(n_col);

std::string dummyline;
std::ifstream myfile (import_path);

size_t line_counter = 0;

if (myfile.is_open())
{
  size_t counter = 0;
  while (getline (myfile,dummyline) )
   {

   if(line_counter == 0) {
   std::istringstream iss(dummyline);

   	do
   	{
   		std::string sub;
   		iss >> sub;
   		if(sub.length()) {
			header[counter] = sub;
			counter++;
   		}
   	}
   	while(iss);
   }

   line_counter++;
   }
 myfile.close();
}
else Rcpp::Rcerr << "Unable to open file";

}

// identify the time column (name time is obligatory and also the position)
bool check_time_column(std::vector<double> time_vector, std::string time) {
  bool name = false;
  bool values_all_numeric = false;
  bool ret = false;
  int counter = 0;
	if(time == "time") {
    name = true;
  } else {
			Rcpp::Rcerr << "Error: First column has to be the Name time" << std::endl;}
	for(size_t i = 0; i < time_vector.size(); i++) {
		if(std::isnan(time_vector[i])) {
		Rcpp::Rcerr << "Error: Time vector is not allowed to contain NAs" << std::endl;
		break;
    }
    counter = counter + 1;
	}

  if(static_cast<unsigned int>(counter) == time_vector.size()) {
      values_all_numeric = true;
  }

  if(name == true && values_all_numeric == true) {
    ret = true;
  }

  return ret;
}

// generate nc without NA
void remove_NA(std::vector<std::vector<double> > &nc_without_NA, std::vector<std::vector<double> > &nc, int nrow, int ncol) {
nc_without_NA.resize(ncol);
std::vector<int> number_rows_without_NA(ncol);

for(int i = 0; i < ncol; i++) {
  for(int j = 0; j < (nrow-1); j++) {
    if(std::isnan(nc[i][j])) {
    }
    else {
      number_rows_without_NA[i] = number_rows_without_NA[i] + 1;
    }
  }
}
for(int i = 0; i < ncol; i++) {
  nc_without_NA[i].resize(number_rows_without_NA[i]);
}
int p = 0;
for(int i = 0; i < ncol; i++) {
  for(int j = 0; j < (nrow-1); j++) {
    if(std::isnan(nc[i][j])) { }
    else {
      nc_without_NA[i][p] = nc[i][j];
      p = p + 1;
    }
  }
  p = 0;
}
}

void Import_Parameter(std::string start, std::string lower, std::string upper,
std::vector<int> &params_cut_idx_vec,
std::vector<double> &params_time_combi_vec,
std::vector<double> &param_combi_start,
std::vector<double> &param_combi_lb,
std::vector<double> &param_combi_ub,
std::vector<std::string> &test) {

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

  for(int i = 0; i < ncol; i++) {
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
    Rcpp::stop("\nERROR: time column of startvalues not correct");
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
  for(int i = 0; i < ncol; i++) { // i = 1
    for(int j = 0; j < (nrow-1); j++) {
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
    Rcpp::stop("\nERROR: Different number of time_points between startvalues, lower bounds and upper bounds");
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
  for(int i = 1; i < ncol; i++) {
    std::vector<double> temp_param(number_rows_without_NA[i]);
    std::vector<double> temp_lb(number_rows_without_NA[i]);
    std::vector<double> temp_ub(number_rows_without_NA[i]);
    std::vector<double> temp_time(number_rows_without_NA[i]);
    for(int j = 0; j < number_rows_without_NA[i]; j++) {
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
  // ===================================================

}

void Import_start_parameter(std::string start, std::vector<int> &params_cut_idx_vec,
std::vector<double> &params_time_combi_vec, std::vector<double> &param_combi_start,
std::vector<std::string> &test) {

  // Check if file exists
  // ===================================================
  bool is_file_start = false;
  is_file_start = check_file_exists(start);
  if(is_file_start == false) {
    Rcpp::stop("\nERROR: No File with start values");
    //exit (EXIT_FAILURE);
  }
  // ===================================================

  // Import start
  // ===================================================
  int ncol, nrow;
  int temp_nrow;
  count_cols_rows(start, nrow, ncol);
  check_ncols_per_row(start, temp_nrow, ncol);
  std::vector<std::vector<double> > nc;
  get_content(start, nc, nrow, ncol);
  get_header(start, test, ncol);
  std::vector<std::vector<double> > nc_without_NA;
  remove_NA(nc_without_NA, nc, nrow, ncol);

  // ===================================================

  // extract time vector
  // ===================================================
  bool time_check = false;
  time_check = check_time_column(nc[0], "time");
  if(time_check == false) {
    Rcpp::stop("\nERROR: time column of states not correct");
    //exit (EXIT_FAILURE);
  }
  std::vector<double> time(nc[0].size());

  for(size_t i = 0; i < time.size(); i++) {
    time[i] = nc[0][i];
  }
  // ===================================================

  // extract number of rows without NA and timepoints
  // ===================================================
  std::vector<int> number_rows_without_NA(ncol);
  std::vector<std::vector<int> > time_points(ncol);
  for(int i = 0; i < ncol; i++) { // i = 1
    for(int j = 0; j < (nrow-1); j++) {
      if(std::isnan(nc[i][j])) { }
      else {
        number_rows_without_NA[i] = number_rows_without_NA[i] + 1;
        time_points[i].push_back(j);
      }
    }
  }
  // ===================================================

  // define combination of parameter and time and store it in std::list
  // ===================================================
  std::list<ParamClass> paramlist;
  for(int i = 1; i < ncol; i++) {
    std::vector<double> temp_param(number_rows_without_NA[i]);
    std::vector<double> temp_time(number_rows_without_NA[i]);
    for(int j = 0; j < number_rows_without_NA[i]; j++) {
      temp_param[j] = nc_without_NA[i][j];;
      temp_time[j] = nc[0][time_points[i][j]];
    }
    if(temp_time.size() == 1) {
      Rcpp::Rcerr << test[i] << ":" << "\t" << "Parameter is const" << std::endl;
      ParamClass temp(temp_param[0]);
      paramlist.push_back(temp);
    } else {
      Rcpp::Rcerr << test[i] << ":" << "\t" << "Parameter is variabel" << std::endl;
      ParamClass temp(temp_time.size(), temp_time, temp_param);
      paramlist.push_back(temp);
    }
  }

  ParamOrderClass ParamOrder(paramlist);
  ParamOrder.cut_idx(params_cut_idx_vec);
  ParamOrder.get_time_combi(params_time_combi_vec);
  ParamOrder.get_param_combi(param_combi_start);
  // ===================================================

}

void Import_states(std::string start,
std::vector<int> &hs_cut_idx_vec,
std::vector<double> &hs_time_combi_vec,
std::vector<double> &hs_harvest_state_combi_vec,
  std::vector<std::string> &test) {

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

  for(int i = 1; i < (ncol); i++) {
    if(std::isnan(nc[i][0])) {
      Rcpp::stop("\nError: first line of states contains NA");
    }
  }

  get_header(start, test, ncol);

  // extract time vector
  bool time_check = false;
  time_check = check_time_column(nc[0], "time");
  if(time_check == false) {
    Rcpp::warning("\nERROR: time column of startvalues not correct");
  }
  std::vector<double> time(nc[0].size());

  for(size_t i = 0; i < time.size(); i++) {
    time[i] = nc[0][i];
  }

  int size_of_states = time.size();
  // define combination of state and time and store it in std::list
  std::list<HarvestStateClass> paramlist;
  for(int i = 1; i < ncol; i++) {
    std::vector<double> temp_state(size_of_states);
    std::vector<double> temp_time(size_of_states);
    for(int j = 0; j < size_of_states; j++) {
      temp_state[j] = nc[i][j];
      temp_time[j] = nc[0][j];
    }
      Rcpp::Rcerr << test[i] << std::endl;
      HarvestStateClass temp(temp_time, temp_state);
      paramlist.push_back(temp);
  }

  HarvestStateOrderClass HarvestOrder(paramlist);
  HarvestOrder.cut_idx(hs_cut_idx_vec);
  HarvestOrder.get_harvest_time_combi(hs_time_combi_vec);
  HarvestOrder.get_harvest_state_combi(hs_harvest_state_combi_vec);
}
