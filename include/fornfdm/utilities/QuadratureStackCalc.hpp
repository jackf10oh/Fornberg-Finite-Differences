// QuadratureStackCalc.hpp
//
// reverses the Fornberg algorithm to compute weights for 
// a definite integral on [x_left, x_right]. 
// computes from nthIntegral repeated integral to 
// numIntegralsMax repeated integrals. 
// most commonly just going 1 repeated integral to 1 repeated integral.
//
// no free function since the memory usage is pretty strange
//
// JAF 8/8/2026

#ifndef FORNFDM_UTILS_QUADRATURESTACKCALC_H
#define FORNFDM_UTILS_QUADRATURESTACKCALC_H

namespace fornfdm{
  namespace utils{

template<std::size_t numNodes_, std::size_t nthIntegral_, class ScalarType_, std::size_t numIntegralsMax_=nthIntegral_>
class QuadratureStackCalc
{
  private:
    // Type Defs ------------ 
    using ScalarType = ScalarType_;

    // Member Data ---------------- 
    static constexpr std::size_t nthIntegral = nthIntegral_; 
    static constexpr std::size_t numIntegralsMax = numIntegralsMax_; 
    static constexpr std::size_t numNodes = numNodes_; 
    std::size_t m_nodes_used = 0;
    std::array<ScalarType, numNodes * (numNodes+numIntegralsMax-nthIntegral)> m_buff_left;
    std::array<ScalarType, numNodes * (numNodes+numIntegralsMax-nthIntegral)> m_buff_right;

  public:
    // Constructors + Destructor =================
    QuadratureStackCalc()=default;
    QuadratureStackCalc(const QuadratureStackCalc& other)=default;
    ~QuadratureStackCalc()=default;

    // Member Functions ------------
    const auto& getArray() const { return m_buff_right; }

    template<class InputIter, class = std::enable_if_t<std::is_base_of_v<std::random_access_iterator_tag, typename std::iterator_traits<InputIter>::iterator_category>>>
    void calculate(ScalarType x_left, ScalarType x_right, InputIter node_start)
    {
      // notation: kth weight using n nodes for mth repeated integral
      // I(m,n,k)

      // Initialize first column I(nthIntegral,1,0), I(nthIntegral+1,1,0), ..., I(numIntegralsMax+numNodes,1,0) 
      std::size_t fact = 1;
      ScalarType xr_powm = std::pow(x_right - *node_start,nthIntegral);
      ScalarType xl_powm = std::pow(x_left - *node_start,nthIntegral);
      for(std::size_t m=0; m < numIntegralsMax+numNodes-nthIntegral; ++m){
        fact *= nthIntegral + m;
        m_buff_right[m * numNodes] = xr_powm / fact; 
        m_buff_left[m * numNodes] = xl_powm / fact;
        xr_powm *= x_right - *node_start;
        xl_powm *= x_left - *node_start;
      }

      ScalarType c1 = 1.0; 
      ScalarType c2; 
      ScalarType c3; 
      // iteratre 2nd until last column 
      for(std::size_t n=1; n < numNodes; ++n)
      {
        c2 = 1.0; 
        // loop through previous columns 
        for(std::size_t k=0; k < n; ++k)
        {
          c3 = node_start[n] - node_start[k]; 
          c2 *= c3; 
          // newest column depends on OLD values of prev column ...
          if(k == n - 1)
          {
            for(std::size_t m=0; m < numIntegralsMax - nthIntegral + numNodes - n; ++m)
            {
              // (3.9) --------------------
              // I(m,n,n)(x) = [(x-a_n)*I(m,n-1,n-1) - (m+1)*I(m,n-1,n)] * (c1/c2)
              m_buff_right[m*numNodes + n] = ((x_right - node_start[n-1])*m_buff_right[m*numNodes+n-1] - (m+nthIntegral)*m_buff_right[(m+1)*numNodes+n-1])*(c1/c2);
              m_buff_left[m*numNodes + n] = ((x_left - node_start[n-1])*m_buff_left[m*numNodes+n-1] - (m+nthIntegral)*m_buff_left[(m+1)*numNodes+n-1])*(c1/c2);
            } 
          }; 
          // updates m = nthIntegral, nthIntegral+1, ..., (numIntegralsMax + numNodes - n) th repeated integral 
          for(std::size_t m=0; m < (numIntegralsMax-nthIntegral+numNodes-n); ++m)
          {
            // (3.8) -----------------
            // I(m,n,k)(x) = [m*I(m+1,n,k-1) - (x-a_n)*I(m,n,k-1)] * (a_n-a_k)
            m_buff_right[m*numNodes+k] = ((m+nthIntegral) * m_buff_right[(m+1)*numNodes+k] - (x_right-node_start[n]) * m_buff_right[m*numNodes+k]) / c3;
            m_buff_left[m*numNodes+k] = ((m+nthIntegral) * m_buff_left[(m+1)*numNodes+k] - (x_left-node_start[n]) * m_buff_left[m*numNodes+k]) / c3;
          }; 
        }
        c1 = c2; 
      }
      std::transform(
        m_buff_right.begin(), std::next(m_buff_right.begin(), (numIntegralsMax-nthIntegral+1)*numNodes),
        m_buff_left.begin(), 
        m_buff_right.begin(),
        std::minus<ScalarType>{}
      );
    }
};

  } // end namespace utils
} // end namespace fornfdm

#endif // QuadratureStackCalc.hpp