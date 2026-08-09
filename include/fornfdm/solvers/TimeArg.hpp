// TimeArg.hpp
//
// Class containing start + end times, and a step size. 
// enforces shared_ptr to match solvers interface.
//
// JAF 8/8/2026 

#include<memory>
#include<cstddef>
#include<iterator>
#include "../types.hpp"

#ifndef FORNFDM_SOLVERS_TIMEARG_H
#define FORNFDM_SOLVERS_TIMEARG_H

namespace fornfdm{
namespace solvers{

class TimeArg
{
  private:
    // Nested Types ---- 
    struct HiddenType{};
  public:
    class Builder
    {
      private:
        // Member Data ----------- 
        std::size_t m_n_steps = 0;
        fornfdm::Real m_t0;
        fornfdm::Real m_t1;
        fornfdm::Real m_dt;
      public:
        // Member Funcs --------
        std::shared_ptr<const TimeArg> build()
        {
          assert((m_t0 < m_t1) && "error t0 must be < t1");
          assert(((m_dt > 0.0) || (m_n_steps>0)) && "error dt must be > 0.0 or num_steps must be > 0");
          if(m_n_steps != 0)
          {
            return std::make_shared<const TimeArg>(m_t0,m_t1,(m_t1-m_t0)/(m_n_steps-1), HiddenType{}); 
          }
          else
          {
            return std::make_shared<const TimeArg>(m_t0,m_t1,m_dt, HiddenType{}); 
          }
        }
        Builder& setStart(fornfdm::Real t) noexcept { m_t0=t; return *this; }
        Builder& setStepSize(fornfdm::Real dt) noexcept { m_n_steps = 0; m_dt=dt; return *this; }
        Builder& setNumSteps(std::size_t n) noexcept { m_n_steps = n; return *this; }
        Builder& setStop(fornfdm::Real t) noexcept { m_t1=t; return *this; }
    };
    struct UniformSpacedIter;
    using value_type = fornfdm::Real;
    using reference = void;
    using const_reference = void;
    using const_iterator = UniformSpacedIter;
    using iterator = const_iterator;
    using difference_type = std::ptrdiff_t;
    using size_type = std::size_t;

  private:
    // Member Data ----------- 
    fornfdm::Real m_start;
    fornfdm::Real m_stepsize;
    std::size_t m_num_steps;
    fornfdm::Real m_stop;
  
  public:
    // Constructors ------
    TimeArg()=delete;
    TimeArg(fornfdm::Real t0, fornfdm::Real t1, fornfdm::Real dt, HiddenType) noexcept
      : m_start(t0), m_stepsize(dt), m_stop(t1), m_num_steps(((t1-t0)/dt)+1)
    {}
    TimeArg(const TimeArg& other)=delete;
    ~TimeArg()=default;

    // Member Funcs --------
    static Builder builder() noexcept { return Builder{}; }
    fornfdm::Real getStart() const noexcept { return m_start; }
    fornfdm::Real getStop() const noexcept { return m_stop; }
    fornfdm::Real getStepSize() const noexcept { return m_stepsize; }
    std::size_t getNumSteps() const noexcept { return m_num_steps; }
    size_type size() const noexcept { return m_num_steps; }
    const_iterator cbegin() const noexcept { return {this, 0}; }
    const_iterator cend() const noexcept 
    { 
      constexpr std::size_t maximum = static_cast<std::size_t>(std::numeric_limits<difference_type>::max());   
      assert((m_num_steps <= maximum) && "fatal error, narrowing std::size_t -> std::ptrdiff_t overflowed"); 
      return {this, static_cast<difference_type>(m_num_steps)}; 
    }

    // Iterator class --------------------------- 
    class UniformSpacedIter
    {
      friend TimeArg;
      public:
        // Type Defs ---------------
        using iterator_category = std::random_access_iterator_tag; 
        using value_type = fornfdm::Real;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = void;
      private:
        // Member Data ------------
        const TimeArg* m_arg;
        difference_type m_offset;
        // Constructors
        UniformSpacedIter(const TimeArg* arg, difference_type offset)
          : m_arg(arg), m_offset(offset)
        {}
      public:
        UniformSpacedIter()=default;
        UniformSpacedIter(const UniformSpacedIter& other)=default;
        ~UniformSpacedIter()=default;
        
        // Operators ------------
        value_type operator*() const { return m_arg->m_start + m_offset * m_arg->m_stepsize; }
        UniformSpacedIter& operator++(){ ++m_offset; return *this; }
        UniformSpacedIter operator++(int){ UniformSpacedIter tmp = *this; ++m_offset; return tmp; }
        UniformSpacedIter& operator+=(difference_type i){ m_offset += i; return *this; }
        UniformSpacedIter operator+(difference_type i) const { return {m_arg, m_offset+i}; }
        friend UniformSpacedIter operator+(difference_type i, const UniformSpacedIter& it){ return {it.m_arg, it.m_offset + i}; }
        UniformSpacedIter& operator--(){ --m_offset; return *this; }
        UniformSpacedIter operator--(int){ UniformSpacedIter tmp = *this; --m_offset; return tmp; }
        UniformSpacedIter& operator-=(difference_type i){ m_offset -= i; return *this; }
        UniformSpacedIter operator-(difference_type i) const { return {m_arg, m_offset - i}; }
        difference_type operator-(const UniformSpacedIter& it) const { return m_offset - it.m_offset; }
        value_type operator[](difference_type i) const { return m_arg->m_start + (m_offset + i) * (m_arg->m_stepsize); }
        bool operator<(const UniformSpacedIter& it) const { return m_offset < it.m_offset; }
        bool operator<=(const UniformSpacedIter& it) const { return m_offset <= it.m_offset; }
        bool operator>(const UniformSpacedIter& it) const { return m_offset > it.m_offset; }
        bool operator>=(const UniformSpacedIter& it) const { return m_offset >= it.m_offset; }
        bool operator==(const UniformSpacedIter& it) const { return m_offset == it.m_offset; }
        bool operator!=(const UniformSpacedIter& it) const { return m_offset != it.m_offset; }
    };
};

} // end namespace solvers
} // end namespace fornfdm

#endif