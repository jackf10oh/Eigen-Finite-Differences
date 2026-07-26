// benchmarks.cpp
//
// JAF 12/12/2025

#include<benchmark/benchmark.h>
// #define FORNFDM_CUSTOM_SCALAR float
#include<fornfdm/all.hpp>
using namespace fornfdm;

static void fornberg_algo(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int n_nodes = state.range(0);  
    int order_m = state.range(1);
    fornfdm::Vector v = fornfdm::linspaced(n_nodes, 0.0, n_nodes - 1);  
    fornfdm::Vector u(n_nodes * (order_m + 1));  
    state.ResumeTiming(); 
    fornfdm::utils::fornberg(
      v.cbegin(), std::next(v.cbegin(),n_nodes), 
      0.0, 
      order_m, 
      u.begin()
    );
  }
};

static void setMesh_order_1_direction_0(benchmark::State &state)
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

static void setMesh_order_2_direction_0(benchmark::State &state)
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

static void setMesh_order_3_direction_0(benchmark::State &state)
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

static void setMesh_saxby(benchmark::State &state)
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

static void setMesh_sum(benchmark::State &state)
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

static void setTime_assignment(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1);
    linops::TimeDepCoeff func = [](fornfdm::Real t, fornfdm::Scalar x){ return t + x*x; };  
    auto D = func * linops::NthPartialDeriv<1,0>{} + linops::NthPartialDeriv<2,0>{}; 
    D.setMesh(mesh); 
    state.ResumeTiming(); 
    D.setTime(10.0); 
    fornfdm::CSRMatrix stencil = D;
    benchmark::DoNotOptimize(stencil);
  }
}; 

static void evalTime_assignment(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r+1,0.0,double(r)), 1);
    linops::TimeDepCoeff func = [](fornfdm::Real t, fornfdm::Scalar x){ return t + x*x; };  
    auto D = func * linops::NthPartialDeriv<1,0>{} + linops::NthPartialDeriv<2,0>{}; 
    D.setMesh(mesh); 
    state.ResumeTiming(); 
    fornfdm::CSRMatrix stencil = D.evalTime(10.0);
    benchmark::DoNotOptimize(stencil);
  }
}; 

template<std::size_t M>
static void Executor_getRhsExpression_assignment(benchmark::State &state)
{
  for(auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);
    texprs::NthTimeDeriv<1> Ut; 
    auto exec = texprs::make_Executor<M+1>(Ut);
    for(int i=0; i< M+1; ++i)
    {
      fornfdm::Scalar t = i;
      exec.pushTime(t);
      fornfdm::Vector s = fornfdm::linspaced(r, (r-1) * i, (r-1) * (i+1));
      exec.pushSolution(s);
    }
    fornfdm::Vector dest(r); // size r. don't want allocation in the benchmark. 
    state.ResumeTiming(); 
    dest = exec.getRhsExpression();
    benchmark::DoNotOptimize(dest);
  }
}

template<std::size_t M>
static void Executor_rotateStoredSolutions(benchmark::State &state)
{
  for(auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);
    texprs::NthTimeDeriv<1> Ut; 
    auto exec = texprs::make_Executor<M+1>(Ut);
    for(int i=0; i< M+1; ++i)
    {
      fornfdm::Scalar t = i;
      exec.pushTime(t);
      fornfdm::Vector s = fornfdm::linspaced(r, (r-1) * i, (r-1) * (i+1));
      exec.pushSolution(s);
    }
    state.ResumeTiming(); 
    exec.rotateStoredSolutions(1);
  }
}

static void explicit_euler_1d(benchmark::State &state)
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
      .times = std::make_shared<const fornfdm::RealVector>(linspaced(m_timesteps,0.0,pi))
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

static void implicit_euler_1d(benchmark::State &state)
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
      .times = std::make_shared<const fornfdm::RealVector>(linspaced(m_timesteps,0.0,pi))
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

static void crank_nicolson_1d(benchmark::State &state)
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
      .times = std::make_shared<const fornfdm::RealVector>(linspaced(m_timesteps,0.0,pi))
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

BENCHMARK(fornberg_algo)->Args({1,2})->Args({1,3})->Args({2,3})->Args({1,4})->Args({2,4})->Args({3,4});
BENCHMARK(setMesh_order_1_direction_0)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(setMesh_order_2_direction_0)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(setMesh_order_3_direction_0)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(setMesh_saxby)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(setMesh_sum)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(setTime_assignment)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(evalTime_assignment)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(Executor_getRhsExpression_assignment<1>)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(Executor_getRhsExpression_assignment<2>)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(Executor_getRhsExpression_assignment<4>)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(Executor_rotateStoredSolutions<1>)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(Executor_rotateStoredSolutions<2>)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(Executor_rotateStoredSolutions<4>)->Arg(100)->Arg(200)->Arg(400)->Arg(800);
BENCHMARK(explicit_euler_1d)->Args({100,100})->Args({100,200})->Args({100,400})->Args({100,800});
BENCHMARK(implicit_euler_1d)->Args({100,100})->Args({100,200})->Args({100,400})->Args({100,800});
BENCHMARK(crank_nicolson_1d)->Args({100,100})->Args({100,200})->Args({100,400})->Args({100,800});

BENCHMARK_MAIN();