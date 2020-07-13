/********************************************************************
//
//  file: param_interpolation.cpp
//
//  contents:
//    - sort function which passes time, timevector
//      and parametervector to spline
//    - Cubic Spline interpolation for variable parameters
//    - Similar to Hermite Spline 
//    - periodic Spline
//
*********************************************************************/

#include "param_interpolation.hpp"

void params_sort (realtype &t, std::vector<double> &params, std::vector<double> &par_vec, std::vector<double> &time_vec, std::vector<int> &param_idx_cuts);
double CatmullRomSpline(realtype &t, std::vector<double> &time_vec,
                         std::vector<double> &par_vec);

void params_sort (
  realtype &t, std::vector<double> &params, std::vector<double> &par_vec, std::vector<double> &time_vec, std::vector<int> &param_idx_cuts){
  // Get number of parameter
  int no_par = param_idx_cuts.size();

  params.resize(no_par);

  std::vector<double> tmp_time_vec;
  std::vector<double> tmp_par_vec;
  int tmp_no_vals;
  int idx_count = 0;

  for (int i = 0; i < no_par; ++i) {
    tmp_no_vals = param_idx_cuts[i];
    if (tmp_no_vals==1){
      params[i] = par_vec[idx_count];
      ++idx_count;
    }else{
      tmp_par_vec.resize(tmp_no_vals);
      tmp_time_vec.resize(tmp_no_vals);
      for (int j = 0; j < tmp_no_vals; ++j) {
        tmp_par_vec[j] = par_vec[idx_count];
        tmp_time_vec[j] = time_vec[idx_count];
        ++idx_count;
      }
      params[i] = CatmullRomSpline(t, tmp_time_vec, tmp_par_vec);
    }
  }
}

// ========================================================================
// ========================================================================
// ========================================================================

double CatmullRomSpline( // (labled with ! need check)
    realtype &t, std::vector<double> &time_vec, std::vector<double> &par_vec
){
  int idx0, idx1, idx2, idx3;
  double t0, t1, t2, t3;
  double y0, y1, y2, y3;

  idx0 = 0;
  idx1 = 0;
  idx2 = 0;
  idx3 = 0;
  t0 = t1 = t2 = t3 = 0.;
  y0 = y1 = y2 = y3 = 0.;
  //if(t > time_vec.back()) { //!
   //t = t  - time_vec.back();
  //}
  //for (signed int i = 0; i <= time_vec.size(); ++i) { // i < time_vec.size !
  for(size_t i = 0; i <= time_vec.size(); i++)  {
    if (i == time_vec.size()) {

      idx0 = time_vec.size() - 2;
      t0 = time_vec[idx0];
      y0 = par_vec[idx0];

      idx1 = time_vec.size() - 1;
      t1 = time_vec[idx1];
      y1 = par_vec[idx1];

      idx2 = time_vec.size()  - time_vec.size();
      t2 = time_vec[idx2];
      y2 = par_vec[idx2];

      idx3 = time_vec.size() + 1 - time_vec.size();
      t3 = time_vec[idx3];
      y3 = par_vec[idx3];
      break;
    }else if (t>=time_vec[i] && t<time_vec[i+1]){
      if (i==0) {
        idx0 = time_vec.size()-1;
        t0 = time_vec[idx0]-time_vec.back();
      } else {
        idx0 = i-1;
        t0 = time_vec[idx0];
      }
      y0 = par_vec[idx0];
      idx1 = i;
      t1 = time_vec[idx1];
      y1 = par_vec[idx1];
      if ( i == time_vec.size()-1 ) {
        idx2 = 0;
        t2 = time_vec[idx2]+time_vec.back();
      } else {
        idx2 = i+1;
        t2 = time_vec[idx2];
      }
      y2 = par_vec[idx2];
      if ( i == time_vec.size()-2 ) {
        idx3 = 0;
        t3 = time_vec[idx3]+time_vec.back();
      } else if ( i == time_vec.size()-1 ) {
        idx3 = 1;
        t3 = time_vec[idx3]+time_vec.back();
      } else {
        idx3 = i+2;
        t3 = time_vec[idx3];
      }
      y3 = par_vec[idx3];
      break;
    }
  } // search for the beginning of the interpolation intervall

  double x = (t -t1) / (t2 -t1);
  double m1 = (y2 -y0) / (t2 -t0);
  double m2 = (y3 -y1) / (t3 -t1);

  double res = (
    (2.*pow(x,3) -3.*pow(x,2) +1.) *y1
    + (pow(x,3) -2.*pow(x,2) +x) *(t2-t1) *m1
    + (-2.*pow(x,3) +3.*pow(x,2)) *y2
    + (pow(x,3) -pow(x,2)) *(t2-t1) *m2
  );
  return res;
}
