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

namespace OSteps{

namespace traits{
struct NA{}; 
} // end namespace traits 

template<typename C=traits::NA>
struct TimeContext
{
  double now; 
  double next; 
  std::shared_ptr<const C> container; 
  TimeContext(double t1=0.0, double t2=0.0, std::shared_ptr<const C> c=nullptr)
    : now(t1), next(t2), container(c)
  {}
}; 

template<typename M=traits::NA, typename X=traits::NA, typename R=traits::NA, typename S=traits::NA>
class Context
{
  private:
    // Member Data --------------------------------------------
    std::shared_ptr<M> mesh;
    X* const executor; 
    R* const rhs_expr;
    S* const solver; 
  public:
    // Contstructor ======================================================= 
    Context(std::shared_ptr<M> m=nullptr, X* const x=nullptr, R* const r=nullptr, S* s=nullptr)
      : mesh(std::move(m)), executor(x), rhs_expr(r), solver(s)
    {}

    // Member Funcs ==================================================
    // const / non const getters. 
    std::shared_ptr<M>& getMesh(){ return mesh; }; 
    const std::shared_ptr<M>& getMesh() const { return mesh; }; 
    X& getExecutor(){ return *executor; }
    const X& getExecutor() const { return *executor; }
    R& getRhsExpr(){ return *rhs_expr; }
    const R& getRhsExpr() const { return *rhs_expr; } 
    S& getSolver(){ return *solver; }
    const S& getSolver() const { return *solver; }  

}; 

// entry point functions
template<typename C=traits::NA>
auto make_time(double t1=0.0, double t2=0.0, std::shared_ptr<const C> c=nullptr){ return TimeContext(t1,t2,c); }; 

template<typename M=traits::NA, typename X=traits::NA, typename R=traits::NA>
auto make_context(std::shared_ptr<const M> m=nullptr, X* x=nullptr, R* r=nullptr){ return Context(m,x,r); }

} // end namespace OSteps 

#endif // StepContexts.hpp 