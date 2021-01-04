#ifndef INTERPOLATION
#define INTERPOLATION
//*********************************************************************
//
//  file: param_interpolation.hpp
//
//  contents:
//    - Cubic Spline interpolation for variable parameters
//    - Similar to Hermite Spline
//
//**********************************************************************
//
#include "header.hpp"

//
// ================================================================= //
void params_sort (
  realtype &t, std::vector<double> &params, std::vector<double> &par_vec, std::vector<double> &time_vec, std::vector<int> &param_idx_cuts
);
// ================================================================= //
double CatmullRomSpline(
  realtype &t, std::vector<double> &time_vec, std::vector<double> &par_vec
);
// ================================================================= //

#endif // INTERPOLATION
