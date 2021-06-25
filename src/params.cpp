/********************************************************************
//
//  file: params.cpp
//
//  contents:
//    - Definition of ParamClass and ParamOrderClass
//
*********************************************************************/
#include "params.hpp"

// ParamClass
// ========================================================================
// ========================================================================
// ========================================================================
ParamClass::ParamClass(
  double t_par_val
){
  // Initialize Constant parameter, which requires one upper and
  // lower bound for the whole time, and no additional time vec
  m_no_spl_pts = 1;

  m_time_vec.resize(1);
  m_time_vec[0] = 0.;
  //
  m_spl_pts_arr.resize(1);
  m_spl_pts_arr[0] = t_par_val;
  //
}
//
ParamClass::ParamClass(
  int t_no_spl_pts,
  std::vector<double> t_time_vec,
  std::vector<double> t_par_vec
){
  if (
      static_cast<unsigned int>(t_no_spl_pts)!=t_time_vec.size()
    || static_cast<unsigned int>(t_no_spl_pts)!=t_par_vec.size()
  ){
    //Rcpp::Rcerr << "[params.h] ParamClass_init() wrong size of arguments." << std::endl;
    Rcpp::stop("\nERROR: ParamClass_init wrong size of arguments.");
    //exit (EXIT_FAILURE);
  }
  //
  m_no_spl_pts=t_no_spl_pts;
  //
  m_time_vec.resize(t_no_spl_pts);
  m_time_vec = t_time_vec;
  //
  m_spl_pts_arr.resize(t_no_spl_pts);
  for ( int j = 0; j < t_no_spl_pts; ++j ) {
    m_spl_pts_arr[j] = t_par_vec[j];
  }
}
//
ParamClass::ParamClass (
    double t_low_bound,
    double t_up_bound
){
  // Initialize Constant parameter, which requires one upper and
  // lower bound for the whole time, and no additional time vec
  m_no_spl_pts = 1;
  //
  m_time_vec.resize(1);
  m_time_vec[0] = 0.;
  //
  if ( t_low_bound>t_up_bound ) {
    //Rcpp::Rcerr << "[params.h] ParamClass_init() boundary value error." << std::endl;
    Rcpp::stop("\nERROR: ParamClass_init boundary value error.");
    //exit (EXIT_FAILURE);
  }
  m_lb_arr.resize(1);
  m_ub_arr.resize(1);
  //
  m_lb_arr[0] = t_low_bound;
  m_ub_arr[0] = t_up_bound;
}
//
ParamClass::ParamClass (
    int t_no_spl_pts,
    std::vector<double> t_time_vec,
    std::vector<double> t_low_bound,
    std::vector<double> t_up_bound
){
  if (
      static_cast<unsigned int>(t_no_spl_pts)!=t_time_vec.size()
    || static_cast<unsigned int>(t_no_spl_pts)!=t_low_bound.size()
    || static_cast<unsigned int>(t_no_spl_pts)!=t_up_bound.size()
  ){
    //Rcpp::Rcerr << "[params.h] ParamClass_init() wrong size of arguments." << std::endl;
    Rcpp::stop("\nERROR: ParamClass_init wrong size of arguments.");
    //exit (EXIT_FAILURE);
  }
  //
  m_no_spl_pts=t_no_spl_pts;
  //
  m_time_vec.resize(t_no_spl_pts);
  m_time_vec = t_time_vec;
  //
  for ( int i = 0; i < t_no_spl_pts; ++i ) {
    if ( t_low_bound[i]>t_up_bound[i] ) {
      //Rcpp::Rcerr << "[params.h] ParamClass_init() boundary value error." << std::endl;
      Rcpp::stop("\nERROR: ParamClass_init boundary value error.");
      //exit (EXIT_FAILURE);
    }
  }

  m_lb_arr.resize(m_no_spl_pts);
  m_ub_arr.resize(m_no_spl_pts);

  for(int i = 0; i < m_no_spl_pts; i++) {
    m_lb_arr[i] = t_low_bound[i];
    m_ub_arr[i] = t_up_bound[i];
  }
}
//
ParamClass::ParamClass (
    double t_par_val,
    double t_low_bound,
    double t_up_bound
){
  m_no_spl_pts = 1;
  //
  m_time_vec.resize(1);
  m_time_vec[0] = 0.;
  //
  if ( t_low_bound>t_up_bound ) {
    //Rcpp::Rcerr << "[params.h] ParamClass_init() boundary value error." << std::endl;
    Rcpp::stop("\nERROR: ParamClass_init boundary value error.");
    //exit (EXIT_FAILURE);
  }
  m_spl_pts_arr.resize(1);
  m_lb_arr.resize(1);
  m_ub_arr.resize(1);
  //
  m_spl_pts_arr[0] = t_par_val;
  m_lb_arr[0] = t_low_bound;
  m_ub_arr[0] = t_up_bound;
}
//
ParamClass::ParamClass (
    int t_no_spl_pts,
    std::vector<double> t_time_vec,
    std::vector<double> t_par_vec,
    std::vector<double> t_low_bound,
    std::vector<double> t_up_bound
){
     if (
      static_cast<unsigned int>(t_no_spl_pts)!=t_time_vec.size()
    || static_cast<unsigned int>(t_no_spl_pts)!=t_low_bound.size()
    || static_cast<unsigned int>(t_no_spl_pts)!=t_up_bound.size()
    || static_cast<unsigned int>(t_no_spl_pts)!=t_par_vec.size()
  ){
    //Rcpp::Rcerr << "[params.h] ParamClass_init() wrong size of arguments." << std::endl;
    Rcpp::stop("\nERROR: ParamClass_init wrong size of arguments.");
    //exit (EXIT_FAILURE);
  }
  //
  m_no_spl_pts=t_no_spl_pts;
  //
  m_time_vec.resize(t_no_spl_pts);
  m_time_vec = t_time_vec;
  //
  for ( int i = 0; i < t_no_spl_pts; ++i ) {
    if ( t_low_bound[i]>t_up_bound[i] ) {
      //Rcpp::Rcerr << "[params.h] ParamClass_init() boundary value error." << std::endl;
      Rcpp::stop("\nERROR: ParamClass_init boundary value error.");
      //exit (EXIT_FAILURE);
    }
  }

  m_lb_arr.resize(m_no_spl_pts);
  m_ub_arr.resize(m_no_spl_pts);
  m_spl_pts_arr.resize(m_no_spl_pts);
  for(int i = 0; i < m_no_spl_pts; i++) {
    m_lb_arr[i] = t_low_bound[i];
    m_ub_arr[i] = t_up_bound[i];
    m_spl_pts_arr[i] = t_par_vec[i];
  }
}

// ParamOrderClass
// ========================================================================
// ========================================================================
// ========================================================================

ParamOrderClass::ParamOrderClass (std::list<ParamClass> t_ParamList){
  m_ParamList = t_ParamList;
}

void ParamOrderClass::c_count(int &t_col_no) {
  //int cerrPOScrc=0;

  std::list<ParamClass>::iterator it;
  t_col_no = 0;

  for (it = m_ParamList.begin(); it != m_ParamList.end(); ++it){
    t_col_no += it->get_no_spl_pts();
  }
}

// ================================================================= //

void ParamOrderClass::cut_idx(std::vector<int> &t_cut_idx_vec){
  std::list<ParamClass>::iterator it;
  int Cmax;
  c_count(Cmax);
  int Cc=0;
  for (it = m_ParamList.begin(); it != m_ParamList.end(); ++it){
    Cc+=1;
  }
  t_cut_idx_vec.resize(Cc);
  Cc=0;
  for (it = m_ParamList.begin(); it != m_ParamList.end(); ++it){
    t_cut_idx_vec[Cc] = it->get_no_spl_pts();
    Cc+=1;
  }
}

// ================================================================= //

void ParamOrderClass::get_time_combi(std::vector<double> &t_time_combi_vec){
  //
  //int cerrPOSgtc=0;

  std::list<ParamClass>::iterator it;
  int Cmax;
  c_count(Cmax);

  t_time_combi_vec.resize(Cmax);

  int Cc=0;
  //
  for (it = m_ParamList.begin(); it != m_ParamList.end(); ++it){
    for (int iC = 0; iC < it->get_no_spl_pts(); ++iC) {
      t_time_combi_vec[Cc] = it->get_time_points(iC);
      Cc+=1;
    }
  }
}

// ================================================================= //

void ParamOrderClass::get_param_combi(std::vector<double> &t_param_combi_vec){
  //
  //int cerrPOSgpc=0;

  std::list<ParamClass>::iterator it;

  int Cmax;
  c_count(Cmax);

  t_param_combi_vec.resize(Cmax);

  int counter = 0;
  for(it = m_ParamList.begin(); it != m_ParamList.end(); ++it) {
    for(int j = 0; j < it -> get_no_spl_pts(); ++j) {

      t_param_combi_vec[counter] = it -> get_spl_pts_arr(j);;
      counter += 1;
    }
  }
}

void ParamOrderClass::get_lb_combi(std::vector<double> &t_lb_combi_vec){
  //
  //int cerrPOSgpc=0;

  std::list<ParamClass>::iterator it;

  int Cmax;
  c_count(Cmax);

  t_lb_combi_vec.resize(Cmax);

  int counter = 0;
  for(it = m_ParamList.begin(); it != m_ParamList.end(); ++it) {
    for(int j = 0; j < it -> get_no_spl_pts(); ++j) {

      t_lb_combi_vec[counter] = it -> get_low_bound(j);;
      counter += 1;
    }
  }
}

void ParamOrderClass::get_up_combi(std::vector<double> &t_ub_combi_vec){
  //
  //int cerrPOSgpc=0;

  std::list<ParamClass>::iterator it;

  int Cmax;
  c_count(Cmax);

  t_ub_combi_vec.resize(Cmax);

  int counter = 0;
  for(it = m_ParamList.begin(); it != m_ParamList.end(); ++it) {
    for(int j = 0; j < it -> get_no_spl_pts(); ++j) {

      t_ub_combi_vec[counter] = it -> get_up_bound(j);;
      counter += 1;
    }
  }
}
