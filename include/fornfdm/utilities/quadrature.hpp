// quadrature.hpp
//
// reverse of the Fornberg algorithm do calculate 
// numerical integral weights instead of finite difference weights. 
// note: the # of repeated integrals is descending.
//
// JAF 8/9/2026

#ifndef FORNFDM_UTILS_QUADRATURE_H
#define FORNFDM_UTILS_QUADRATURE_H

std::pair<double*, double*> quadrature(
  const double* node_start, const double* node_end,  
  double x_left, double x_right, 
  std::size_t min_m, std::size_t max_m,
  double* dest_left, double* dest_right 
)
{



  return {} // returns a std::pair<double*,double*> of 2 "one past the end" output iterators
}

#endif