// Wave2D.hpp
//
// Bindings for Wave equation in 2D. 
// Forcing term in Center oscilates according to sin(t) 
// and 4 walls damp reflection using Robin Boundary Conditions 
//
// JAF 2/2/2026 

#ifndef WAVE2D_H
#define WAVE2D_H

#include<LinOps/All.hpp>
#include<OutsideSteps/All.hpp> 
#include<TExprs/All.hpp>

#include<Utilities/BumpFunc.hpp> // forcing term at origin 

class origin_bump : public OSteps::OStepBaseXD<origin_bump>
{
  public:  
    // Member Data ----------------------------
    BumpFunc bump_1d = BumpFunc{.L = -1.0, .R = 1.0, .c =0.0, .h = 1.0, .focus=10}; 
  private:
    const double* m_inv_coeff_ptr = nullptr; 

  public:
    // Member Funcs ------------------------------------------------
    // Calls SetTime on any LinOps inside the expression
    template<typename TIME_ITER, typename LHS_EXECUTOR, typename RHS_EXPR, OSteps::FDStep_Type Step>
    void BeforeLinAlgebra(TIME_ITER& time_iter, LinOps::MeshXD_SPtr_t& mesh, LHS_EXECUTOR& exec, RHS_EXPR& rhs_expr)
    { 
      m_inv_coeff_ptr = &exec.inv_coeff(); 
    }

    template<OSteps::FDStep_Type STEP>
    void MatBeforeStep(double t, const LinOps::MeshXD_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat) const 
    { /* do nothing to stencil...*/}

    // Adds a 2d bump at origin to solution before step. regardless of implicit or explicit -> IMEX scheme
    template<OSteps::FDStep_Type STEP>
    void SolBeforeStep(double t, const LinOps::MeshXD_SPtr_t& mesh, OSteps::StridedRef_t Sol) const 
    {
      auto forcing = [&](double x, double y){ return bump_1d(x)*bump_1d(y)*std::sin(t); }; 
      Eigen::VectorXd disc_vals = LinOps::make_Discretization(mesh, forcing).values(); 
      disc_vals *=  (*m_inv_coeff_ptr); 
      Sol += disc_vals; 
    }

    template<OSteps::FDStep_Type STEP>
    void SolAfterStep(double t, const LinOps::MeshXD_SPtr_t& mesh, OSteps::StridedRef_t Sol) const 
    { /* do nothing to solution after step...*/} 
}; 

struct Wave2D_impl
{
  // Utt 
  using Lhs_t = TExprs::NthTimeDeriv;  
  Lhs_t Lhs = TExprs::NthTimeDeriv(2);
  
  // = Uxx + Uyy 
  LinOps::DirectionalNthDerivOp Uxx = LinOps::DirectionalNthDerivOp(2,0);  
  LinOps::DirectionalNthDerivOp Uyy = LinOps::DirectionalNthDerivOp(2,1);

  using Rhs_t = decltype(Uxx + Uyy); 
  Rhs_t Rhs = Uxx + Uyy; 

  // One outside step for forcing term... 
  origin_bump forcing_term = origin_bump{};

  // double damping = 2.0;
  using BC_t = OSteps::BCList< OSteps::BCPair<OSteps::RobinBC, OSteps::RobinBC>, OSteps::BCPair<OSteps::RobinBC,OSteps::RobinBC> >;
  BC_t boundary_condition = BC_t( 
    OSteps::BCPair(
      OSteps::RobinBC(1.0,-1.0,0.0), // left 
      OSteps::RobinBC(1.0,1.0,0.0) // right 
    ), 
    OSteps::BCPair(
      OSteps::RobinBC(1.0,-1.0,0.0), //similar 
      OSteps::RobinBC(1.0,1.0,0.0)
    )
  ); 

  using OStep_t = decltype(std::tie(forcing_term, boundary_condition));   
  OStep_t osteps = std::tie(forcing_term, boundary_condition); 

}; 

struct Wave2D : public Wave2D_impl, public TExprs::GenInterp<Wave2D_impl::Lhs_t, Wave2D_impl::Rhs_t, Wave2D_impl::OStep_t, LinOps::MeshXD_SPtr_t>
{
  using Interp_t = TExprs::GenInterp<Wave2D_impl::Lhs_t, Wave2D_impl::Rhs_t, Wave2D_impl::OStep_t, LinOps::MeshXD_SPtr_t>; 
  Wave2D()
    :Wave2D_impl(), 
    Interp_t(this->Lhs, this->Rhs, this->osteps)
  {
    // std::get<0>(this->osteps).m_inv_coeff_ptr = &this->Lhs.inv_coeff();
  };  

  void set_damping(double damping)
  {
    // osteps are std::tie(forcing_term, bcs) 
    auto& damping_bcs = std::get<1>(this->osteps); 

    // sets dim1 to new damping 
    auto& dim_1 = std::get<0>(damping_bcs.m_list); 
    dim_1.m_left.val_coeff = damping; 
    dim_1.m_right.val_coeff = damping; 

    // sets dim1 to new damping 
    auto& dim_2 = std::get<1>(damping_bcs.m_list); 
    dim_2.m_left.val_coeff = damping; 
    dim_2.m_right.val_coeff = damping; 
  }

  void set_bump_height(double height)
  {
    // osteps are std::tie(forcing_term, bcs) 
    auto& forcing_term = std::get<0>(this->osteps);
    
    forcing_term.bump_1d.h = std::sqrt(height); 
  }
}; 

#endif 


