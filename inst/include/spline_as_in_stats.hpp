#ifndef SPLINESTATS
#define SPLINESTATS

#include "header.hpp"


static void natural_spline(int n, double *x, double *y, double *b, double *c, double *d);


static void spline_coef(int n, double *x, double *y,
	                      double *b, double *c, double *d);


static void spline_eval(int nu, double *u, double *v,
	          int n, double *x, double *y, double *b, double *c, double *d);

double wrapper_spline(double t, std::vector<double> &time_vec, std::vector<double> &par_vec);
#endif // SPLINESTATS
