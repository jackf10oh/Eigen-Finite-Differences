// basic_api.cpp
//
// JAF 8/31/2026

#include<iostream>
#include<memory>
#include<fornfdm/all.hpp>

using namespace fornfdm; 
using namespace fornfdm::linops; 

using std::endl, std::cout; 

struct len{ static constexpr double value=2.0; };

int main()
{
  // uniform spacing on [0,10], 1D mesh
  auto m = std::make_shared<Mesh>(linspaced(11,0.0,10.0), 1);
  cout << m->getAxis(0).transpose() << "\n";

  NthPartialDeriv<1,0> Ux{};    // order 1, direction 0 
  Ux.setMesh(m);                // constructs 11 x 11 stencil 
  cout << Ux.toEigen() << "\n"; // cast to Eigen::SparseMatrixBase<...>

  NthPartialDeriv<2,0> Uxx{};   // order 2, direction 0
  Uxx.setMesh(m);
  cout << Uxx.toEigen() << "\n";

  Eigen::SparseMatrix<double,Eigen::RowMajor> rhs;

  // naive calculation
  rhs = Uxx.toEigen() + 0.5 * Ux.toEigen();
  cout << rhs << "\n";

  // using new expressions!
  auto xpr = NthPartialDeriv<2,0>{} + 0.5 * NthPartialDeriv<1,0>{}; 
  xpr.setMesh(m);
  cout << xpr.toEigen() << "\n";

  // and more!
  auto m02 = std::make_shared<Mesh>(linspaced(6,0.0,5.0), /*dim=*/2);
  NthPartialDeriv<2,1> Uyy{};
  TimeDepCoeff c = [](double t){ return t*t; };
  auto xpr02 = c * Uyy; 
  xpr02.setMesh(m02); 
  xpr02.setTime(2.0);
  cout << xpr02.getStencil() << "\n";
  cout << xpr02.toEigen() << "\n";

  // getters
  std::shared_ptr<const Mesh> m03 = xpr02.getMesh(); 
  double t = xpr02.getTime();

  // node_selector tags
  NthPartialDeriv<1,0,Forward<4>> Ux_forward{};     // 4 nodes to the right 
  // struct len{ static constexpr double value=2.0; }; // before main()
  NthPartialDeriv<1,0,Periodic<len>> Ux_periodic{}; // wraps around the axis! 
};