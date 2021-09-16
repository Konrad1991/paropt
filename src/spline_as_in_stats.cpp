/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998--2021  The R Core Team
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*	Spline Interpolation
 *	--------------------
 *	C code to perform spline fitting and interpolation.
 *	There is code here for:
 *
 *	1. Natural splines.
 *
 *	2. Periodic splines
 *
 *	3. Splines with end-conditions determined by fitting
 *	   cubics in the start and end intervals (Forsythe et al).
 *
 *
 *	Computational Techniques
 *	------------------------
 *	A special LU decomposition for symmetric tridiagonal matrices
 *	is used for all computations, except for periodic splines where
 *	Choleski is more efficient.
 */


 /*
modified for paropt in order to use it from C++
 */


//#include "spline_as_in_stats.hpp"
#include <vector>
#include <iostream>
#include <memory>

/*
 *	Natural Splines
 *	---------------
 *	Here the end-conditions are determined by setting the second
 *	derivative of the spline at the end-points to equal to zero.
 *
 *	There are n-2 unknowns (y[i]'' at x[2], ..., x[n-1]) and n-2
 *	equations to determine them.  Either Choleski or Gaussian
 *	elimination could be used.
 */


static void natural_spline(int n, double *x, double *y, double *b, double *c, double *d) {

if(n < 2) {
   return;
}

x--; y--; b--; c--; d--;

if(n < 3) {
    double t = (y[2] - y[1]);
    b[1] = t / (x[2]-x[1]);
    b[2] = b[1];
    c[1] = c[2] = d[1] = d[2] = 0.0;
    return;
}

const int nm1 = n-1;
int i;

/* Set up the tridiagonal system */
/* b = diagonal, d = offdiagonal, c = right hand side */
d[1] = x[2] - x[1];
c[2] = (y[2] - y[1])/d[1];

for( i=2 ; i<n ; i++) {
     d[i] = x[i+1] - x[i];
     b[i] = 2.0 * (d[i-1] + d[i]);
     c[i+1] = (y[i+1] - y[i])/d[i];
     c[i] = c[i+1] - c[i];
}

/* Gaussian elimination */
for(i=3 ; i<n ; i++) {
	double t = d[i-1]/b[i-1];
	b[i] = b[i] - t*d[i-1];
	c[i] = c[i] - t*c[i-1];
}

/* Backward substitution */
c[nm1] = c[nm1]/b[nm1];
for(i=n-2 ; i>1 ; i--) {
  c[i] = (c[i]-d[i]*c[i+1])/b[i];
}


/* End conditions */
c[1] = c[n] = 0.0;

/* Get cubic coefficients */
b[1] = (y[2] - y[1])/d[1] - d[i] * c[2];
c[1] = 0.0;
d[1] = c[2]/d[1];
b[n] = (y[n] - y[nm1])/d[nm1] + d[nm1] * c[nm1];
for(i=2 ; i<n ; i++) {
	b[i] = (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2.0*c[i]);
	d[i] = (c[i+1]-c[i])/d[i];
	c[i] = 3.0*c[i];
}
c[n] = 0.0;
d[n] = 0.0;

  return;
}


// 1. spline coef
/* These were/are the public interfaces */
static void spline_coef(int n, double *x, double *y,
	                      double *b, double *c, double *d) {

natural_spline(n, x, y, b, c, d);
}

// 2. spline eval

static void spline_eval(int nu, double *u, double *v,
	          int n, double *x, double *y, double *b, double *c, double *d) {

/* Evaluate  v[l] := spline(u[l], ...),	    l = 1,..,nu, i.e. 0:(nu-1)
 * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
*/
const int n_1 = n - 1;
int i, l;
double dx;


for(l = 0; l < nu; l++) v[l] = u[l];

  for(l = 0, i = 0; l < nu; l++) {
	   double ul = v[l];
	    if(ul < x[i] || (i < n_1 && x[i+1] < ul)) {
	       /* reset i  such that  x[i] <= ul <= x[i+1] : */
	        i = 0;
	        int j = n;
	         do {
		             int k = (i+j) / 2;
		             if(ul < x[k]) j = k; else i = k;
	         } while(j > i+1);
	    }
	dx = ul - x[i];

	/* for natural splines extrapolate linearly left */
	double tmp = (ul < x[0]) ? 0.0 : d[i];
	v[l] = y[i] + dx*(b[i] + dx*(c[i] + dx*tmp));
  }

}


double wrapper_spline(double t, std::vector<double> &time_vec, std::vector<double> &par_vec) {

  double* tvec = time_vec.data();
  double* pvec = par_vec.data();

  int n = time_vec.size();

  double* b;
  double* c;
  double* d;
  double* yout;
  b = new double[n];
  c = new double[n];
  d = new double[n];
  yout = new double[n];

  spline_coef(n, tvec, pvec, b, c, d);

  spline_eval(n, &t, yout, n, tvec, pvec, b, c, d);

  double ret = yout[0];

  delete[] b;
  delete[] c;
  delete[] d;
  delete[] yout;

  return ret;
}
