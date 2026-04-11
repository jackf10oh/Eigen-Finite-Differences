// CvecDiffuseReac.cpp
//
// First attempt at binding PDEs for python use
//
// JAF 1/24/2026 

/*
-I"C:/eigen-install/include/eigen3" -I"C:\Users\jackf\anaconda3\Lib\site-packages\pybind11\include" -I"C:\Users\jackf\anaconda3\include"
*/
 
class ConvectionDiffusionAdvectionPDE 
{
  public:
    double a; 
    double b; 
    double c; 

  private:
    LinOps::IOp U = LinOps::IOp();  
    Fds::NthDerivOp Ux = Fds::NthDerivOp(1);  
    Fds::NthDerivOp Uxx = Fds::NthDerivOp(2);     
    texprs::NthTimeDeriv Ut = texprs::NthTimeDeriv(1); 

    decltype(a*Uxx + b*Ux + c*U) m_assembled; 

  private:
    using lhs_expr_t = decltype(Ut); 
    using rhs_expr_t = decltype(m_assembled); 
    using domain_mesh_ptr_t = std::shared_ptr<const LinOps::Mesh1D>
    using bc_ptr_t = Fds::BcPtr_t; 
    using ic_cont_t = std::vector<Eigen::VectorXd>; 
    using interp_t = texprs::GenInterp<lhs_expr_t, rhs_expr_t, domain_mesh_ptr_t, bc_ptr_t, ic_cont_t>; 
    using solver_args_t = GenSolverArgs<domain_mesh_ptr_t, bc_ptr_t, ic_cont_t>; 
    interp_t m_interpolator;  


  public:
    ConvectionDiffusionAdvectionPDE(double a_, double b_, double c_)
      : a(a_), b(b_), c(c_), 
      m_assembled(a*Uxx + b*Ux + c*U), 
      m_interpolator(Ut, m_assembled)
    {};

  const auto& StoredData() const { return m_interpolator.StoredData(); }; 

  double SolAt(double t, double x){ return m_interpolator.SolAt(t,x); }; 

  const auto& Args() const { return m_interpolator.Args(); }; 
  void SetArgs(const solver_args_t& args){ m_interpolator.SetArgs(args); }; 

  void SetDomain(domain_mesh_ptr_t d_ptr){ m_interpolator.Args().domain_mesh_ptr = d_ptr; m_interpolator.Reset(); }; 

  void SetDomain(std::shared_ptr<const LinOps::Mesh1D> t_ptr){ m_interpolator.Args().time_mesh_ptr = t_ptr; m_interpolator.Reset(); }; 

  void SetBCs(bc_ptr_t bc_ptr){ m_interpolator.Args().bcs = bc_ptr; m_interpolator.Reset(); }; 

  void SetICs(ic_cont_t cont){ m_interpolator.Args().ICs = std::move(cont); m_interpolator.Reset(); }; 

  void SetTimeDep(bool b){m_interpolator.Args().time_dep_flag=b; m_interpolator.Reset(); }; 

}; 
