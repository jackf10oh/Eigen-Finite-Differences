// benchmarks.cpp
//
//
//
// JAF 12/12/2025

// compilation command 
// g++ -O3 -DNDEBUG     ./benchmarks/benchmarks.cpp  -lbenchmark -lpthread     -o benchmarks_main

#include<benchmark/benchmark.h>

#include<fornfdm/all.hpp>
using namespace fornfdm;

static void BENCHMARK_setMesh_order_1_direction_0(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1); 
    linops::NthPartialDeriv<1,0> D; 
    state.ResumeTiming(); 
    D.setMesh(mesh); 
  }
};

static void BENCHMARK_setMesh_order_2_direction_0(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1); 
    linops::NthPartialDeriv<2,0> D; 
    state.ResumeTiming(); 
    D.setMesh(mesh); 
  }
};

static void BENCHMARK_setMesh_order_3_direction_0(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1); 
    linops::NthPartialDeriv<3,0> D; 
    state.ResumeTiming(); 
    D.setMesh(mesh); 
  }
};

static void BENCHMARK_setMesh_saxby(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1); 
    auto D = 4.0 * linops::NthPartialDeriv<1,0>{} + linops::NthPartialDeriv<2,0>{}; 
    state.ResumeTiming(); 
    D.setMesh(mesh); 
  }
}; 

static void BENCHMARK_setMesh_sum(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1); 
    auto D = linops::NthPartialDeriv<1,0>{} + linops::NthPartialDeriv<2,0>{}; 
    state.ResumeTiming(); 
    D.setMesh(mesh); 
  }
}; 

static void BENCHMARK_explicit_euler_1d(benchmark::State &state)
{
  for(auto _ : state)
  {
    state.PauseTiming(); 
    double pi = 3.14159265385;
    int n_gridpoints = state.range(0);
    int m_timesteps = state.range(1);

    // uniform mesh from 0.0 to r with n_gridpoints
    solvers::SolverArgs args{
      .mesh = make_Mesh(linspaced(n_gridpoints,0.0,pi),1),
      .times = std::make_shared<const Eigen::VectorXd>(linspaced(m_timesteps,0.0,pi))
    };
    args.initialConditions = { discretize(args.mesh, [](fornfdm::Scalar x){ return std::sin(x); }) };

    // Left hand side in time
    auto Ut = texprs::NthTimeDeriv<1>{};

    // Right hand side in space 
    auto Uxx = linops::NthPartialDeriv<2,0>{};

    // set the boundary conditions to Dirichlet 0
    auto left = osteps::Dirichlet(0.0);
    auto right = left; 
    osteps::BCPair bcs{left, right};

    solvers::ExplicitSolver solver(Ut, Uxx, std::tie(bcs));

    state.ResumeTiming();
    benchmark::DoNotOptimize(solver.calculate(std::move(args),solvers::LastSaver{}));
  }
}

static void BENCHMARK_implicit_euler_1d(benchmark::State &state)
{
  for(auto _ : state)
  {
    state.PauseTiming(); 
    double pi = 3.14159265385;
    int n_gridpoints = state.range(0);
    int m_timesteps = state.range(1);

    // uniform mesh from 0.0 to r with n_gridpoints
    solvers::SolverArgs args{
      .mesh = make_Mesh(linspaced(n_gridpoints,0.0,pi),1),
      .times = std::make_shared<const Eigen::VectorXd>(linspaced(m_timesteps,0.0,pi))
    };
    args.initialConditions = { discretize(args.mesh, [](fornfdm::Scalar x){ return std::sin(x); }) };

    // Left hand side in time
    auto Ut = texprs::NthTimeDeriv<1>{};

    // Right hand side in space 
    auto Uxx = linops::NthPartialDeriv<2,0>{};

    // set the boundary conditions to Dirichlet 0
    auto left = osteps::Dirichlet(0.0);
    auto right = left; 
    osteps::BCPair bcs{left, right};

    solvers::ImplicitSolver solver(Ut, Uxx, std::tie(bcs));

    state.ResumeTiming();
    benchmark::DoNotOptimize(solver.calculate(std::move(args),solvers::LastSaver{}));
  }
}

static void BENCHMARK_crank_nicolson_1d(benchmark::State &state)
{
  for(auto _ : state)
  {
    state.PauseTiming(); 
    double pi = 3.14159265385;
    int n_gridpoints = state.range(0);
    int m_timesteps = state.range(1);

    // uniform mesh from 0.0 to r with n_gridpoints
    solvers::SolverArgs args{
      .mesh = make_Mesh(linspaced(n_gridpoints,0.0,pi),1),
      .times = std::make_shared<const Eigen::VectorXd>(linspaced(m_timesteps,0.0,pi))
    };
    args.initialConditions = { discretize(args.mesh, [](fornfdm::Scalar x){ return std::sin(x); }) };

    // Left hand side in time
    auto Ut = texprs::NthTimeDeriv<1>{};

    // Right hand side in space 
    auto Uxx = linops::NthPartialDeriv<2,0>{};

    // set the boundary conditions to Dirichlet 0
    auto left = osteps::Dirichlet(0.0);
    auto right = left; 
    osteps::BCPair bcs{left, right};

    solvers::CrankNicolsonSolver solver(Ut, Uxx, std::tie(bcs));

    state.ResumeTiming();
    benchmark::DoNotOptimize(solver.calculate(std::move(args),solvers::LastSaver{}));
  }
}

BENCHMARK(BENCHMARK_setMesh_order_1_direction_0)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_setMesh_order_2_direction_0)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_setMesh_order_3_direction_0)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_setMesh_saxby)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_setMesh_sum)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(BENCHMARK_explicit_euler_1d)->Args({100,100})->Args({100,200})->Args({100,400})->Args({100,800});
BENCHMARK(BENCHMARK_implicit_euler_1d)->Args({100,100})->Args({100,200})->Args({100,400})->Args({100,800});
BENCHMARK(BENCHMARK_crank_nicolson_1d)->Args({100,100})->Args({100,200})->Args({100,400})->Args({100,800});

BENCHMARK_MAIN();