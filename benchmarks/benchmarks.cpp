// benchmarks.cpp
//
//
//
// JAF 12/12/2025

// compilation command 
// g++ -O3 -DNDEBUG     ./benchmarks/benchmarks.cpp  -lbenchmark -lpthread     -o benchmarks_main

#include<benchmark/benchmark.h>

#include "../All.hpp"

static void BENCHMARK_NthDerivOp_set_mesh(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    std::size_t order_n = 1; 
    auto my_mesh = make_Mesh1D(0.0,double(r), r+1); 
    NthDerivOp D(order_n); 
    state.ResumeTiming(); 
    D.set_mesh(my_mesh); 
  }
};

static void BENCHMARK_ConvectionDiffusionExpression_set_mesh(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = state.range(0);

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_Mesh1D(0.0,r,n_gridpoints);

    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    double dt = 0.01;
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    state.ResumeTiming();
    fdm_scheme.set_mesh(my_mesh);
  };
}

static void BENCHMARK_ConvectionDiffusionExpression_single_explicit_step(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = state.range(0);

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_Mesh1D(0.0,r,n_gridpoints);

    Discretization1D my_vals;
    
    // Bump centered at r/2. Zero at 0.0 and r 
    auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
    
    // fill my_vals with func(x)
    my_vals.set_init(my_mesh, func); 
    
    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    double dt = 0.01;
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    // set the boundary conditions to Dirichlet 0
    auto left = std::make_shared<DirichletBC>(0.0);
    auto right = left; 
    fdm_scheme.lbc_ptr = left; 
    fdm_scheme.rbc_ptr = right; 

    // fill out the stencil matrix
    fdm_scheme.set_mesh(my_mesh);

    state.ResumeTiming();
    benchmark::DoNotOptimize(fdm_scheme.explicit_step(my_vals));

  };
}

static void BENCHMARK_ConvectionDiffusionExpression_N_explicit_steps(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = 1000;

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_Mesh1D(0.0,r,n_gridpoints);

    Discretization1D my_vals;
    
    // Bump centered at r/2. Zero at 0.0 and r 
    auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
    
    // fill my_vals with func(x)
    my_vals.set_init(my_mesh, func); 
    
    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    int NSteps =  state.range(0); 
    double T = 10.0;
    int dt = T/NSteps; 
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    // set the boundary conditions to Dirichlet 0
    auto left = std::make_shared<DirichletBC>(0.0);
    auto right = left; 
    fdm_scheme.lbc_ptr = left; 
    fdm_scheme.rbc_ptr = right; 

    // fill out the stencil matrix
    fdm_scheme.set_mesh(my_mesh);

    state.ResumeTiming();
    for(int n=0; n<NSteps; n++)
    {
      // explcit euler 
      my_vals = std::move(fdm_scheme.explicit_step(my_vals).values()); 
      // my_vals = fdm_scheme.explicit_step(my_vals); // copies an empty shared_ptr<Mesd1D> to my_vals 
    }
  };
}

static void BENCHMARK_ConvectionDiffusionExpression_N_Manual_explicit_steps(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = 1000;

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_Mesh1D(0.0,r,n_gridpoints);

    Discretization1D my_vals;
    
    // Bump centered at r/2. Zero at 0.0 and r 
    auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
    
    // fill my_vals with func(x)
    my_vals.set_init(my_mesh, func); 
    
    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    int NSteps =  state.range(0); 
    double T = 10.0;
    int dt = T/NSteps; 
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    // set the boundary conditions to Dirichlet 0
    auto left = std::make_shared<DirichletBC>(0.0);
    auto right = left; 
    fdm_scheme.lbc_ptr = left; 
    fdm_scheme.rbc_ptr = right; 

    // fill out the stencil matrix
    fdm_scheme.set_mesh(my_mesh);
    const Matrix stencil = fdm_scheme.asMatrix(); 

    state.ResumeTiming();

    // Nsteps explcit euler steps: don't defer to fdm_stencil's .explicit_step() method
    for(int n=0; n<NSteps; n++)
    {
      // use move construct new set of values from result Eigen matrix * Eigen vector
      my_vals = stencil * my_vals.values();  
      // my_vals = fdm_scheme.explicit_step(my_vals); // copies an empty shared_ptr<Mesd1D> to my_vals 
      
      // use BCs manually. 
      left->SetSolL(my_vals.values(),my_mesh); 
      right->SetSolR(my_vals.values(),my_mesh); 
    }

  };
}

static void BENCHMARK_ConvectionDiffusionExpression_single_solve_implicit(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming(); 
    auto r = 10.0; 
    int n_gridpoints = state.range(0);

    // uniform mesh from 0.0 to r with n_gridpoints
    MeshPtr_t my_mesh = make_Mesh1D(0.0,r,n_gridpoints);

    Discretization1D my_vals;
    
    // Bump centered at r/2. Zero at 0.0 and r 
    auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
    
    // fill my_vals with func(x)
    my_vals.set_init(my_mesh, func); 
    
    // create an fdm_scheme for convection diffusion equation 
    using D = NthDerivOp; 
    double dt = 0.01;
    auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));

    // set the boundary conditions to Dirichlet 0
    auto left = std::make_shared<DirichletBC>(0.0);
    auto right = left; 
    fdm_scheme.lbc_ptr = left; 
    fdm_scheme.rbc_ptr = right; 

    // fill out the stencil matrix
    fdm_scheme.set_mesh(my_mesh);

    state.ResumeTiming();
    benchmark::DoNotOptimize(fdm_scheme.solve_implicit(my_vals));

  };
}

BENCHMARK(BENCHMARK_NthDerivOp_set_mesh)->Arg(1e2)->Arg(1e3)->Arg(1e4)->Arg(1e5);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_set_mesh)->Arg(1e2)->Arg(1e3)->Arg(1e4);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_single_explicit_step)->Arg(1e2)->Arg(1e3)->Arg(1e4);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_N_explicit_steps)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_N_Manual_explicit_steps)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_ConvectionDiffusionExpression_single_solve_implicit)->Arg(1e2)->Arg(1e3)->Arg(1e4);

BENCHMARK_MAIN();