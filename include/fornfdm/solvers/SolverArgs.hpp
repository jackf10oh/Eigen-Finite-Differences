// SolverArgs.hpp
//
// P.O.D. class containing mesh of time, mesh of space, and initial conditions. 
// Initial conditions are specified by first N solution at first N entries of time. 
//
// JAF 3/4/2026 

#ifndef FORNFDM_SOLVERS_SOLVERARGS_H
#define FORNFDM_SOLVERS_SOLVERARGS_H

#include<cassert>
#include<memory>
#include<vector>
#include<Eigen/Core> 
#include "../types.hpp"
#include "../Mesh.hpp"

namespace fornfdm{
  namespace solvers{ 

template<class M, class Container>
struct SolverArgs
{
  // Mesh1D or MeshXD the PDE operates on 
  std::shared_ptr<M> mesh; 

  // list of times the solver marches through 
  std::shared_ptr<const Container> times; 
  
  // first N solution values. 
  // ! has to be atleast >= max_order + 1. 
  // where max_order is the highest order in the LHS time derivatives expression (texprs)  
  // defaulted to empty so it can be assigned later... 
  std::vector<fornfdm::Vector> initialConditions = {};
};

// CTAD guideline ... 
template<class M, class C>
SolverArgs(std::shared_ptr<M>, std::shared_ptr<const C>, std::vector<fornfdm::Vector>)
  ->SolverArgs<M, C>; 

template<class M, class C>
SolverArgs(std::shared_ptr<M>, std::shared_ptr<const C>)
  ->SolverArgs<M, C>; 

class TimeArg
{
  private:
    // Nested Types ---- 
    struct HiddenType{};
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
        Builder& setStart(fornfdm::Real t){ m_t0=t; return *this; }
        Builder& setStepSize(fornfdm::Real dt){ m_n_steps = 0; m_dt=dt; return *this; }
        Builder& setNumSteps(std::size_t n){ m_n_steps = n; return *this; }
        Builder& setStop(fornfdm::Real t){ m_t1=t; return *this; }
    };

    // Member Data ----------- 
    fornfdm::Real m_start;
    fornfdm::Real m_stepsize;
    std::size_t m_num_steps;
    fornfdm::Real m_stop;
  
  public:
    // Constructors ------
    TimeArg()=delete;
    TimeArg(fornfdm::Real t0, fornfdm::Real t1, fornfdm::Real dt, HiddenType)
      : m_start(t0), m_stepsize(dt), m_stop(t1), m_num_steps(((t1-t0)/dt)+1)
    {}
    TimeArg(const TimeArg& other)=delete;
    ~TimeArg()=default;

    // Member Funcs --------
    static Builder builder(){ return Builder{}; }

    fornfdm::Real getStart() const { return m_start; }
    fornfdm::Real getStop() const { return m_stop; }
    fornfdm::Real getStepSize() const { return m_stepsize; }
    fornfdm::Real getNumSteps() const { return m_num_steps; }
};

  } // end namespace solvers
} // end namespace fornfdm 

#endif /// SolverArgs.hpp 