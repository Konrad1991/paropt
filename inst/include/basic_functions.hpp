#ifndef BASICFUNCTIONS
#define BASICFUNCTIONS

#include "header.hpp"

// check wether file exists
bool check_file_exists(std::string import_path);

// count columns of text file
int counte_cols(std::string line);

// count rows and columns of text file
void count_cols_rows(std::string import_path, int &nrow, int &ncol);

//Check if each row contains the same number of columns
void check_ncols_per_row(std::string import_path, int &nrow, int &ncol);

// import matrix (not the first line = headers) into std::vector<std::vector<double>>
void get_content(std::string import_path, std::vector<std::vector<double> > &numeric_content, int n_row, int n_col);

// import headers into std::vector<std::string>
void get_header (std::string import_path, std::vector<std::string> &header, int n_col);


// identify the time column (name time is obligatory and also the position)
bool check_time_column(std::vector<double> time_vector, std::string time);

// generate nc without NA
void remove_NA(std::vector<std::vector<double> > &nc_without_NA, std::vector<std::vector<double> > &nc, int nrow, int ncol);

void Import_Parameter(std::string start, std::string lower, std::string upper,
std::vector<int> &params_cut_idx_vec,
std::vector<double> &params_time_combi_vec,
std::vector<double> &param_combi_start,
std::vector<double> &param_combi_lb,
std::vector<double> &param_combi_ub,
std::vector<std::string> &test);

void Import_start_parameter(std::string start, std::vector<int> &params_cut_idx_vec,
std::vector<double> &params_time_combi_vec, std::vector<double> &param_combi_start,
std::vector<std::string> &test);

void Import_states(std::string start,
std::vector<int> &hs_cut_idx_vec,
std::vector<double> &hs_time_combi_vec,
std::vector<double> &hs_harvest_state_combi_vec,
std::vector<std::string> &test);

#endif // BASICFUNCTIONS
