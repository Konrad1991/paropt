#ifndef TYPES
#define TYPES

#include "header.hpp"
typedef int (*OS)(double &t, std::vector<double> &params, std::vector<double> &states);

typedef void (*OS2)(double t, double* params, int size_params, double* y_, double* ydot_, int NEQ);

#endif // TYPES
