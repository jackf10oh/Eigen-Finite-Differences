// StepContexts.hpp
//
// 1. ) POD struct for 
// current time, 
// immediate next time, 
// and (optionally) container of times.
// 
// 2. ) Wrapper class for 
// shared_ptr<> to mesh domain in space
// raw pointer to TExprExecutor
// raw pointer to RHS linear operator expression
// possibly add raw pointer to entire solver later?  
// 
// JAF 3/21/2026 

#ifndef STEPCONTEXTS_H
#define STEPCONTEXTS_H

#include<cstddef>
#include "../Types.hpp"
#include "../Mesh.hpp"

namespace fdm{
  namespace osteps{    

// namespace traits{
// struct NA{}; 
// } // end namespace traits 

template<typename C = void>
struct TimeContext
{
  double now; 
  double next; 
  std::shared_ptr<const C> container; 
  TimeContext(double t1=0.0, double t2=0.0, std::shared_ptr<const C> c=nullptr)
    : now(t1), next(t2), container(c)
  {}
  TimeContext(const TimeContext& other)=delete; 
}; 

template<typename M = fdm::Mesh, typename X=void, typename R=void, typename S=void>
class Context
{
  private:
    using MCleaned = std::remove_cv_t<std::remove_reference_t<M>>; 
    static_assert(std::is_same<fdm::Mesh, MCleaned>::value, "shared_ptr<> passed to osteps::Context must point to Mesh or const Mesh");

    // Member Data --------------------------------------------
    const std::shared_ptr<M> m_mesh;
    X* const m_executor; 
    R* const m_rhs_expr;
    S* const m_solver; 
  public:
    // Contstructor ======================================================= 
    Context(std::shared_ptr<M> m=nullptr, X* x=nullptr, R* r=nullptr, S* s=nullptr)
      : m_mesh(std::move(m)), m_executor(x), m_rhs_expr(r), m_solver(s)
    {}

    Context(const Context& other) = delete; 

    ~Context()=default; 
    
    // Member Funcs ==================================================
    // const / non const getters. 
    M* getMesh(){ return m_mesh.get(); }
    const M* getMesh() const { return m_mesh.get(); }
    std::shared_ptr<M> getMeshShared(){ return m_mesh; } 
    std::shared_ptr<const M> getMeshShared() const { return m_mesh; }; 
    X* getExecutor(){ return m_executor; }
    const X* getExecutor() const { return m_executor; }
    R* getRhsExpr(){ return m_rhs_expr; }
    const R* getRhsExpr() const { return m_rhs_expr; } 
    S* getSolver(){ return m_solver; }
    const S* getSolver() const { return m_solver; }  
}; 

// entry point functions
template<typename C=void>
auto make_time(double t1=0.0, double t2=0.0, std::shared_ptr<const C> c=nullptr){ return TimeContext(t1,t2,c); }; 

template<typename M = fdm::Mesh, typename X=void, typename R=void, typename S=void>
auto make_context(std::shared_ptr<M> m=nullptr, X* x=nullptr, R* r=nullptr, S* s=nullptr){ return Context(m,x,r,s); }

  } // end namespace osteps 
} // end namespace fdm 

#endif // StepContexts.hpp 