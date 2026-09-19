// QuadratureStackCalc.hpp
//
// the quadrature analog of FornbergStackCalc
//
// JAF 8/8/2026

#include<array>
#include<cmath>
#include<cassert>
#include "../types.hpp"

#ifndef FORNFDM_UTILS_QUADRATURESTACKCALC_H
#define FORNFDM_UTILS_QUADRATURESTACKCALC_H

namespace fornfdm{
  namespace utils{

template<std::size_t M, std::size_t N, typename ScalarType_ = fornfdm::Scalar>
class QuadratureStackCalc
{
  public:
    // Type Defs ------------ 
    using ScalarType = ScalarType_;

    // Member Data ---------------- 
    static constexpr std::size_t numNodesMax = M; 
    static constexpr std::size_t nthIntegral = N; 

  private:
    using ArrayType = std::array<ScalarType, M*(N+1)>; 
    ArrayType m_arr;          // single allocation of memory rows*cols big 
    std::size_t m_nodes_used = 0;

  public:
    // Constructors + Destructor =================
    QuadratureStackCalc()=default;
    QuadratureStackCalc(const QuadratureStackCalc& other)=default;

    template<typename Iter>
    QuadratureStackCalc(ScalarType x_left, ScalarType x_right, Iter start, Iter end)
    { 
      calculate(x_left,x_right,start,end); 
    }

    ~QuadratureStackCalc()=default;

    // Member Functions ==================================== 

    // Const getter to stored weights 
    const auto& getArray() const { return m_arr; }
    auto getNumNodesUsed() const { return m_nodes_used; }  

    // Updates m_arr to contain weights up to nth repeated integral
    template<typename Iter>
    void calculate(ScalarType x_left, ScalarType x_right, Iter start, Iter end)
    {
      using T = ArrayType::pointer;
      constexpr std::size_t N = numNodesMax * (nthIntegral + numNodesMax - 1) + 1;
      // make sure distance(start,end) <= numNodesMax 
      auto d = std::distance(start,end); 
      assert((d <= numNodesMax) && "QuadratureStackCalc error: distance(start,end) > numNodesMax");  
      m_nodes_used = d; 
      fornfdm::utils::quadrature<Iter,T,N>(start, end, x_left, x_right, N, m_arr.data());
    }
};

  } // end namespace utils
} // end namespace fornfdm

#endif // QuadratureStackCalc.hpp