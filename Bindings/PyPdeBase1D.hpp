// PyPdeBase1D.hpp
//
//
//
// JAF 1/27/2026

#ifndef PYPDEBASE1D_H
#define PYPDEBASE1D_H

#include<LinOps/All.hpp>
#include<OutsideSteps/All.hpp> 
#include<TExprs/GenInterp.hpp>

template<typename Derived>
struct PDE_Base_Impl
{
  // Constructors + Destructor ---------------------------------------
  PDE_Base_Impl()=default; 
  PDE_Base_Impl(const PDE_Base_Impl& other) = delete; // non copyable!
  ~PDE_Base_Impl()=default; 
  // Member Funcs ---------------------------------------  
  auto& GetLhs(){return static_cast<Derived*>(this)->GetLhs();} 
  auto& GetRhs(){return static_cast<Derived*>(this)->GetRhs();} 
}; 

/* Concrete PDE type that inherits from GenInterpolator 
and stores a PDE_Base_Impl as member data. Upon construction 
it Binds the GenInterpolators lhs and rhs pointers 
to PDE_Base_Impl's lhs and rhs  
*/
template<typename PDE_IMPL>
class Concrete_PDE_1D : public PDE_IMPL, public texprs::GenInterp
<
  typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetLhs())>::type, 
  typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetRhs())>::type, 
  LinOps::SharedConstMesh1D, 
  Fds::BcPtr_t, 
  std::vector<Eigen::VectorXd>
>
{
  private:
    // Type Defs -----------
    using interp_base_t = typename texprs::GenInterp<typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetLhs())>::type, typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetRhs())>::type, LinOps::SharedConstMesh1D, Fds::BcPtr_t, std::vector<Eigen::VectorXd>>; 
    // Member Data ------------
  public:
    Concrete_PDE_1D()
      : PDE_IMPL(), interp_base_t(PDE_IMPL::GetLhs(), PDE_IMPL::GetRhs())
    {}
}; 

/* same ideas. just changes Mesh1D to XD, Bc to BcXD*/
template<typename PDE_IMPL>
class Concrete_PDE_XD : public PDE_IMPL, public texprs::GenInterp
<
  typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetLhs())>::type, 
  typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetRhs())>::type, 
  LinOps::SharedConstMeshXD, 
  Fds::BcXDPtr_t, 
  std::vector<Eigen::VectorXd>
>
{
  private:
    // Type Defs -----------
    using interp_base_t = typename texprs::GenInterp<typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetLhs())>::type, typename std::remove_reference<decltype(std::declval<PDE_IMPL>().GetRhs())>::type, LinOps::SharedConstMeshXD, Fds::BcXDPtr_t, std::vector<Eigen::VectorXd>>; 
    // Member Data ------------
  public:
    Concrete_PDE_XD()
      : PDE_IMPL(), interp_base_t(PDE_IMPL::GetLhs(), PDE_IMPL::GetRhs())
    {}
}; 

/* example of how to use bases. 
doubles like diffusion,convection,etc are captured by reference.
Altering them influences the stored equations. 

struct HeatPDE_impl : public PDE_Base_Impl<HeatPDE_impl>
{ 
  double diffusion = 1.0; 
  double convection = 0.0; 
  double reaction = 0.0; 

  LinOps::IOp U = LinOps::IOp(); 
  Fds::NthDerivOp Ux = Fds::NthDerivOp(1); 
  Fds::NthDerivOp Uxx = Fds::NthDerivOp(2); 

  // rhs expr in space 
  decltype(diffusion*Uxx + convection*Ux + reaction*U) rhs_expr = diffusion*Uxx + convection*Ux + reaction*U;  
  
  // lhs expr in time 
  texprs::NthTimeDeriv Ut = texprs::NthTimeDeriv(1); 

  auto& GetLhs(){ return Ut; } 
  auto& GetRhs(){ return rhs_expr; } 
}; 

using HeatPDE = Concrete_PDE_1D<HeatPDE_impl>; 

  HeatPDE pde; 
  
  pde.Args().domain_mesh_ptr = LinOps::make_Mesh1D(0.0, 10.0, 17); 
  pde.Args().time_mesh_ptr = LinOps::make_Mesh1D(0.0,5.0,21); 
  pde.Args().bcs = std::make_shared<Fds::BCPair>(Fds::make_dirichlet(0.0), Fds::make_dirichlet(0.0)); 

  LinOps::Discretization1D d; 
  BumpFunc f{.L=0.0, .R=10.0, .c = 5.0, .h = 3.0}; 
  d.set_init(pde.Args().domain_mesh_ptr, f); 
  pde.Args().ICs = std::vector<Eigen::VectorXd>{ std::move(d.values()) };
  pde.Args().time_dep_flag = false;

  pde.Reset(); 
  pde.fillStoredSolutions(); 

  print_mat(pde.StoredData(), "Solution"); 

*/ 

#endif // PyPdeBase1D.hpp 