// Quadrature.hpp
//
// Algorithm to compute the weights in a numerical integration. 
// similar concepts to Fornber algorithm. 
//
// JAF 3/13/2026 

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "print.hpp"

template<typename Container>
std::vector<double> QuadratureWeights(double x_left, double x_right, Container nodes, std::size_t M = 1)
{
  // calculates numerical integration weights at {a0, a1, ..., a N-1} in nodes 
  // for the M th repeated integral 

  // number of nodes 
  std::size_t N = nodes.size(); 

  // 2D array. 
  // M + N entries for nodes a0, M + N - 1 for a1, ..., M+N-k for node k 
  // currently just a (M+N) X N matrix 
  std::vector<std::vector<double>> weights_mat_right{ M-1+N, std::vector<double>(N, 0.0) }; 
  std::vector<std::vector<double>> weights_mat_left{ M-1+N, std::vector<double>(N, 0.0) }; 

  // notation: kth weight using n nodes for mth repeated integral
  // I(m,n,k)

  // FIRST COLUMN 
  // initialize I(1,1,0), I(2,1,0), ..., I(M + N,1,0) 
  // the mth repeated integral for just 1 node
  std::size_t fact = 1; 
  for(std::size_t m=0; m < M+N-1; ++m){
    fact *= (m+1); 
    weights_mat_right.at(m).at(0) = std::pow(x_right - nodes.at(0), m+1) / fact; 
    weights_mat_left.at(m).at(0) = std::pow(x_left - nodes.at(0), m+1) / fact; 
  }
  // print_mat(weights_mat_right, "right mat"); 

  double c1 = 1.0; 
  double c2; 
  double c3; 
  // 2ND UNTIL LAST COLUMN 
  for(std::size_t n=1; n < N; ++n)
  {
    std::cout << "USING " << n+1 << "NODES =======" << std::endl; 
    c2 = 1.0; 
    // loop through previous columns 
    for(std::size_t k=0; k < n; ++k)
    {
      // std::cout << "node " << k << ": " << std::endl; 
      c3 = nodes.at(n)-nodes.at(k); 
      // std::cout << "a good" << std::endl; 
      c2 *= c3; 
      
      // newest column depends on OLD values of prev column ...
      if(k == n - 1)
      {
        // newest column has 
        for(std::size_t m=0; m < M - 1 + N - n; ++m)
        {
          // (3.9) 
          weights_mat_right.at(m).at(n) = ((x_right-nodes.at(n-1))*weights_mat_right.at(m).at(n-1) - (m+1)*weights_mat_right.at(m+1).at(n-1)) * (c1/c2); 
          // std::cout << "b good" << std::endl; 

          weights_mat_left.at(m).at(n) = ((x_left-nodes.at(n-1))*weights_mat_left.at(m).at(n-1) - (m+1)*weights_mat_left.at(m+1).at(n-1)) * (c1/c2); 
          // std::cout << "c good" << std::endl;           
        } 
      }; 
      // updates m = 1,2,..., M + N - n th repeated integral 
      // std::size_t m = (M-1+N) - n; // m < ; ++m
      for(std::size_t m=0; m < (M-1+N)-n; ++m)
      {
        // std::cout << m << ", "; 
        // (3.8) I(m,n,k)(x) = [ (x-a_k)*I(m+1,n-1,k)(x) - (m-1)*I(m+1,n-1,k)(x) ] / [a_k - a_n]
        weights_mat_right.at(m).at(k) = - ((x_right-nodes.at(n)) * weights_mat_right.at(m).at(k) - (m+1) * weights_mat_right.at(m+1).at(k)) / c3; 
        weights_mat_left.at(m).at(k) = - ((x_left-nodes.at(n)) * weights_mat_left.at(m).at(k) - (m+1) * weights_mat_left.at(m+1).at(k)) / c3; 
      }; 
      // std::cout << "\n";
      // std::cout << std::endl; 
    }
    c1 = c2; 
  }
  // return weights_mat.at(M-1);
  for(int m=0; m<M-1+N; ++m){
    std::transform(
      weights_mat_right.at(m).begin(), weights_mat_right.at(m).end(),
      weights_mat_left.at(m).begin(), 
      weights_mat_right.at(m).begin(), 
      [](double r, double l){ return r-l; }
    );
  }

  fornfdm::utils::print_mat(weights_mat_right, "diff mat"); 
  // return right - left 
  std::vector<double> result = std::move(weights_mat_right.at(M-1)); 
  return result;  
}

#endif