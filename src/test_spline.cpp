/*
R package paropt
Copyright (C) 2021 Konrad Krämer

This file is part of R package paropt


paroptis free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with paropt
If not see: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html#SEC4
*/

#include <testthat.h>
#include "param_interpolation.hpp"



context("spline") {

  std::vector<double> parameter = {2., 4., 8., 16., 32.};
  std::vector<double> time = {1., 2., 3., 4., 5.};

  std::vector<double> time_output = {1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.};


  std::vector<double> true_res = { 2.000000,  2.861111,  4.000000,  5.60000,  8.000000, 11.277778, 16.000000, 26.36111, 32.000000};



  test_that("spline") {
    for(int i = 0; i < time_output.size(); i++) {
      double temp_time = static_cast<realtype>(time_output[i]);
      double temp0 = CatmullRomSpline(temp_time, time, parameter);
      double temp1 = round(true_res[i]*10)/10;
      temp0 = round(temp0*10)/10;
        expect_true(temp1 == temp0);
    }
  }
}








context("spline2") {

  std::vector<double> parameter = {23.78077, 12.13902, 16.55716, 20.57256, 19.39014, 14.65934, 14.10898, 20.27352, 11.73665, 12.27773};
  std::vector<double> time = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

  std::vector<double> time_output = {0.0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5,  6.0,  6.5,  7.0,  7.5,  8.0,  8.5,  9.0,  9.5, 10.0};


  std::vector<double> true_res = { 23.5, 23.7, 23.780766, 19.850645, 12.139024, 13.515224, 16.557159, 19.134152, 20.572556, 20.637245, 19.390141, 17.111683, 14.659339, 13.143470, 14.108976, 17.981299, 20.273520, 17.080506, 11.736654,  10.8000, 12.277733};



  test_that("spline2") {
    for(int i = 0; i < time_output.size(); i++) {
      double temp_time = static_cast<realtype>(time_output[i]);
      double temp0 = CatmullRomSpline(temp_time, time, parameter);
      double temp1 = round(true_res[i]*10)/10;
      temp0 = round(temp0*10)/10;
      double diff = std::abs(temp0 - temp1);
        expect_true(diff <= 1);
    }
  }
}












context("spline3") {

  std::vector<double> parameter = {15.042140,  7.556236, 12.508205,  6.637941,  3.791856,  3.308695,  6.336402, 11.480164, 17.184272, 12.215001, 20.494376, 10.235471, 16.808018};
  std::vector<double> time = {0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24.};

  std::vector<double> time_output = {0.0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5,  6.0,  6.5,  7.0,  7.5,  8.0,  8.5,  9.0,  9.5, 10.0};


  std::vector<double> true_res = { 15.042140, 15.513836,  12.43,  9.134323,  7.556236,
      8.864886, 10.514551, 11.923050, 12.508205, 11.870119, 10.338031,  8.423464,  6.637941, 5.380095,  4.597005,
      4.122862,  3.791856,  3.480764,  3.236724,  3.149459,  3.308695};



  test_that("spline3") {
    for(int i = 0; i < time_output.size(); i++) {
      double temp_time = static_cast<realtype>(time_output[i]);
      double temp0 = CatmullRomSpline(temp_time, time, parameter);
      double temp1 = round(true_res[i]*10)/10;
      temp0 = round(temp0*10)/10;
      double diff = std::abs(temp0 - temp1);
        expect_true(diff <= 1);
    }
  }
}
