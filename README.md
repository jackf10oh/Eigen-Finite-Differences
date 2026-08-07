<a id="readme-top"></a>
<!-- PROJECT LOGO -->
<br />
<div align="center">
<h3 align="center">Fornberg-Finite-Differences</h3>

<p align="center">
    C++ finite difference library for arbitrary grid domains and time steps.
 </p>
 
<p align="center">
 <img src="wave.png">
</p>
    

</div>

<!-- ABOUT THE PROJECT -->
<h3>About The Project</h3>
Fornberg-Finite-Differences is a comprehensive fdm library made to handle partial derivative equations in any dimension. It uses Eigen's C++ matrix math library extensively and adds to SparseMatrixBase through Eigen's plugin macro: EIGEN_SPARSEMATRIXBASE_PLUGIN. any libraries already making use of this macro can be added through the new macro EIGEN_SPARSEMATRIXBASE_PLUGIN_OTHER.

<h3>Features</h3>

 - any number of dimensions 1-4
 - any rectangular grid domain
 - any subdivision through time
 - expressions of linear operators in space
 - expressions of time derivatives of any order 
 - autonomous and time dependent coefficients 
 - Explicit Euler, Implicit Euler, Crank Nicolson solvers
 - Dirichlet, Neumann, Robin boundary conditions
 - external forcing terms
 - extensible outside_steps framework to allow arbitrary modifications between steps of solvers

<!-- GETTING STARTED -->
<h3>Prerequisites</h3>

You will need to install Eigen3  for C++. It is the only dependency of Fornberg-Finite-Differences.

[Downloading Eigen](https://libeigen.gitlab.io/#download)

<h3>Download</h3>

```
git clone https://github.com/jackf10oh/Fornberg-Finite-Differences
```

<!-- USAGE EXAMPLES -->
<h3>Usage</h3>

```cpp
#include<fornfdm/plugin.hpp> 
#include<fornfdm/all.hpp>
int main()
{
  // Domain + Times
  fornfdm::solvers::SolverArgs args{
    /*.mesh=*/ make_Mesh(fornfdm::linspaced(20,0.0,3.14159), 1), // equal spaced on [0,pi] 
    /*.times=*/ std::make_shared<const fornfdm::RealVector>(fornfdm::linspaced(100,0.0,0.5)) 
  }; 

  // Initial Conditions
  auto v = fornfdm::discretize(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { v }; 

  // LHS in time
  fornfdm::texprs::NthTimeDeriv<1> Ut; 

  // RHS in space
  fornfdm::linops::NthPartialDeriv<2,0> Uxx; 

  // Boundary Conditions
  fornfdm::osteps::Dirichlet left(0.0), right(0.0);
  fornfdm::osteps::BCPair bcs(left,right); 

  // Solving. prints to std::cout.
  fornfdm::solvers::ExplicitSolver my_solver(Ut,Uxx,std::tie(bcs));  
  my_solver.calculate(args, fornfdm::solvers::PrintSaver{});
}
```

<!-- LICENSE -->
<h3>License</h3>

Distributed under the MIT License. See `LICENSE.txt` for more information.

<!-- CONTACT -->
<h3>Contact</h3>
Jack Feds - jackf10oh@gmail.com

Project Link: [https://github.com/jackf10oh/Fornberg-Finite-Differences](https://github.com/jackf10oh/Fornberg-Finite-Differences)

<!-- ACKNOWLEDGMENTS -->
<h3>Acknowledgments</h3>

[Fornberg (1988)](https://ww2.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0935077-0/S0025-5718-1988-0935077-0.pdf)

[Eigen](https://gitlab.com/libeigen/eigen)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
<!-- Shields.io badges. You can a comprehensive list with many more badges at: https://github.com/inttter/md-badges -->
