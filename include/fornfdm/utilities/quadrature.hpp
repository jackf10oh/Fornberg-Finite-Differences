// quadrature.hpp
//
// reverse of the Fornberg algorithm to calculate 
// numerical integral weights instead of finite difference weights. 
//
// note: typically the function returns a 1 past the end
// of where the result is written to.
// if the function is given inputs that would exceed the std::array's 
// it will return the "dest" iterator without advancing it any. 
//
// JAF 8/9/2026

#include<array>
#include<vector>
#include<iterator>
#include<cstdint>
#include<type_traits>

#ifndef FORNFDM_UTILS_QUADRATURE_H
#define FORNFDM_UTILS_QUADRATURE_H

namespace fornfdm{
  namespace utils{
 
template<typename ForwardIn, typename OutputIt>
OutputIt allocating_quadrature(
  ForwardIn node_start, ForwardIn node_end,
  const typename std::iterator_traits<ForwardIn>::value_type& x_left,  
  const typename std::iterator_traits<ForwardIn>::value_type& x_right,  
  std::size_t M_repeats, 
  OutputIt dest 
)
{
  using Scalar = typename std::iterator_traits<ForwardIn>::value_type; 
  std::size_t num_nodes = std::distance(node_start, node_end);
  std::vector<Scalar> left_vec(num_nodes * (M_repeats+num_nodes-1) + 1); 
  std::vector<Scalar> right_vec(num_nodes * (M_repeats+num_nodes-1) + 1);

  // first column gets set to (left_x - a(0)) ^ i / (i!) for i=0,...,M_repeats
  left_vec[0] = 1;
  right_vec[0] = 1;
  std::size_t fact = 1; 
  Scalar left_exp = 1.0;
  Scalar right_exp = 1.0;
  for(auto row_idx = 1; row_idx < M_repeats+num_nodes; ++row_idx)
  {
    fact = fact * row_idx;
    left_exp = left_exp * (x_left - *node_start);
    right_exp = right_exp * (x_right - *node_start);
    left_vec[row_idx * num_nodes] = left_exp / fact;
    right_vec[row_idx * num_nodes] = right_exp / fact;
  }

  Scalar c1, c2, c3 = 1.0; 
  std::size_t current_M = std::max(M_repeats, M_repeats + num_nodes - 1);
  for(auto node=std::next(node_start); node!=node_end; ++node)
  {
    std::size_t count = 0;
    c2 = 1.0;
    for(auto old_node = node_start; old_node!=node; ++old_node, ++count)
    {
      c1 = *node - *old_node;
      c2 *= c1;

      if(std::next(old_node) == node)
      {
        std::size_t m = current_M;
        do
        {
          --m;
          // (3.9) -----------------------
          std::size_t idx = m * num_nodes + count + 1;
          left_vec[idx] = c3 * ((x_left - *old_node) * left_vec[idx-1] - m * left_vec[idx+num_nodes-1]) / c2;
          right_vec[idx] = c3 * ((x_right - *old_node) * right_vec[idx-1] - m * right_vec[idx+num_nodes-1]) / c2;
        } while (m != 0);
      }
      std::size_t m = current_M;
      do
      {
        --m;
        // (3.8) -----------------------
        std::size_t idx = m * num_nodes + count;
        left_vec[idx] = (m * left_vec[idx + num_nodes] - (x_left - *node)* left_vec[idx]) / c1;
        right_vec[idx] = (m * right_vec[idx + num_nodes] - (x_right - *node) * right_vec[idx]) / c1;
      } while (m != 0);
    }
    c3 = c2;
    if(current_M > M_repeats)
    {
      --current_M;
    }
  }
  // return the actual difference of both vectors. 
  auto end = std::next(right_vec.cbegin(), num_nodes * (M_repeats+1));
  return std::transform(
    right_vec.cbegin(), end,
    left_vec.cbegin(), 
    dest, 
    [](const Scalar& r, const Scalar& l){ return r-l; }
  );
}

template<typename ForwardIn, typename OutputIt, std::size_t N=100>
OutputIt quadrature(
  ForwardIn node_start, ForwardIn node_end,
  const typename std::iterator_traits<ForwardIn>::value_type& x_left,  
  const typename std::iterator_traits<ForwardIn>::value_type& x_right,  
  std::size_t M_repeats, 
  OutputIt dest 
)
{
  using Scalar = typename std::iterator_traits<ForwardIn>::value_type; 

  std::size_t num_nodes = std::distance(node_start, node_end);
  if((num_nodes * (M_repeats+num_nodes-1) + 1) > N){
    return allocating_quadrature<ForwardIn,OutputIt>(
      node_start,node_end,
      x_left, x_right,
      M_repeats,
      dest
    );
  }

  std::array<Scalar,N> left_arr{}; 
  std::array<Scalar,N> right_arr{};

  // first column gets set to (left_x - a(0)) ^ i / (i!) for i=0,...,M_repeats
  left_arr[0] = 1;
  right_arr[0] = 1;
  std::size_t fact = 1; 
  Scalar left_exp = 1.0;
  Scalar right_exp = 1.0;
  for(auto row_idx = 1; row_idx < M_repeats+num_nodes; ++row_idx)
  {
    fact = fact * row_idx;
    left_exp = left_exp * (x_left - *node_start);
    right_exp = right_exp * (x_right - *node_start);
    left_arr[row_idx * num_nodes] = left_exp / fact;
    right_arr[row_idx * num_nodes] = right_exp / fact;
  }

  Scalar c1, c2, c3 = 1.0; 
  std::size_t current_M = std::max(M_repeats, M_repeats + num_nodes - 1);
  for(auto node=std::next(node_start); node!=node_end; ++node)
  {
    std::size_t count = 0;
    c2 = 1.0;
    for(auto old_node = node_start; old_node!=node; ++old_node, ++count)
    {
      c1 = *node - *old_node;
      c2 *= c1;

      if(std::next(old_node) == node)
      {
        std::size_t m = current_M;
        do
        {
          --m;
          // (3.9) -----------------------
          std::size_t idx = m * num_nodes + count + 1;
          left_arr[idx] = c3 * ((x_left - *old_node) * left_arr[idx-1] - m * left_arr[idx+num_nodes-1]) / c2;
          right_arr[idx] = c3 * ((x_right - *old_node) * right_arr[idx-1] - m * right_arr[idx+num_nodes-1]) / c2;
        } while (m != 0);
      }
      std::size_t m = current_M;
      do
      {
        --m;
        // (3.8) -----------------------
        std::size_t idx = m * num_nodes + count;
        left_arr[idx] = (m * left_arr[idx + num_nodes] - (x_left - *node)* left_arr[idx]) / c1;
        right_arr[idx] = (m * right_arr[idx + num_nodes] - (x_right - *node) * right_arr[idx]) / c1;
      } while (m != 0);
    }
    c3 = c2;
    if(current_M > M_repeats)
    {
      --current_M;
    }
  }
  // return the actual difference of both vectors. 
  auto end = std::next(right_arr.cbegin(), num_nodes * (M_repeats+1));
  return std::transform(
    right_arr.cbegin(), end,
    left_arr.cbegin(), 
    dest, 
    [](const Scalar& r, const Scalar& l){ return r-l; }
  );
}

  } // end namespace utils
} // end namespace fornfdm 

#endif