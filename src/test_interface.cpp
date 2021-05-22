/* !!revision!!
What is this?
*/


#include "header.hpp"


// [[Rcpp::interfaces(r, cpp)]]
int test_interfaces(int a) {
  return a*2;
}
