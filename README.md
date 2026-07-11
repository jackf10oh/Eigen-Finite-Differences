<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a id="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
![Contributors](https://img.shields.io/github/contributors/jackf10oh/Fornberg-Finite-Differences.svg?style=for-the-badge)
![Forks](https://img.shields.io/github/forks/jackf10oh/Fornberg-Finite-Differences.svg?style=for-the-badge)
![Stargazers](https://img.shields.io/github/stars/jackf10oh/Fornberg-Finite-Differences.svg?style=for-the-badge)
![Issues](https://img.shields.io/github/issues/jackf10oh/Fornberg-Finite-Differences.svg?style=for-the-badge)
![project_license](https://img.shields.io/github/license/jackf10oh/Fornberg-Finite-Differences.svg?style=for-the-badge)
![LinkedIn](https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555)

<!-- PROJECT LOGO -->
<br />
<div align="center">
<h3 align="center">Fornberg-Finite-Differences</h3>
<p align="center">
    C++ fdm library for arbitrary grid domains and time steps.
 </p>
</div>

<p align="center">
 <img src="wave.png">
</p>

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
// workaround to allow expressions of different rows/cols to be added 
// before setMesh() sets their rows/cols to be equal
#define eigen_assert(x)
#include<fornfdm/all.hpp>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <fornfdm/plugin.hpp> 

#include<iostream>
#include<iomanip>
#include<fornfdm/utilities/print.hpp> 
#include<fornfdm/utilities/BumpFunc.hpp>
#include<Eigen/SparseCore> // macro plugin takes effect. 

using namespace fornfdm; 

int main()
{
  // Domain + Times  
  constexpr double pi = 3.14159265385; 
  fornfdm::solvers::SolverArgs args{
    .mesh = make_Mesh(fornfdm::linspaced(20,0.0,pi), 1), 
    .times = std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(100,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = fornfdm::discretize(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0>{}; 

  // Boundary Conditions 
  auto left = osteps::Dirichlet(0.0); 
  auto right = left;
  osteps::BCPair bcs(left,right); 

  // Solving. prints to std::cout.
  solvers::ExplicitSolver my_solver(Ut,Uxx,std::tie(bcs));  
  my_solver.calculate(args, solvers::PrintSaver{});
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
