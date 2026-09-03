// benchmarks.cpp
//
// JAF 12/12/2025

#include<benchmark/benchmark.h>
#include<fornfdm/all.hpp>
using namespace fornfdm;

template<std::size_t nthOrder>
static void setMesh_single(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    linops::NthPartialDeriv<nthOrder,0> deriv; 
    state.ResumeTiming(); 
    deriv.setMesh(mesh); 
    benchmark::DoNotOptimize( deriv );
  }
};

static void setMesh_expresion(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    auto deriv = linops::NthPartialDeriv<2,0>{} + 0.5 * linops::NthPartialDeriv<1,0>{}; 
    state.ResumeTiming(); 
    deriv.setMesh(mesh); 
    benchmark::DoNotOptimize( deriv );
  }
};

static void naive_assign(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    linops::NthPartialDeriv<1,0> u_x;
    linops::NthPartialDeriv<2,0> u_xx; 
    u_x.setMesh(mesh); 
    u_xx.setMesh(mesh);

    CSRMatrix dest;
    auto s = mesh->getAxis(0).size();
    dest.resize(s, s);
    dest.reserve(3 * s);
    state.ResumeTiming(); 
    dest = u_xx.getStencil() + 0.5 * u_x.getStencil();
    benchmark::DoNotOptimize( dest );
  }
};

static void expression_assign(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    auto xpr = linops::NthPartialDeriv<2,0>{} + 0.5 * linops::NthPartialDeriv<1,0>{};
    xpr.setMesh(mesh);

    CSRMatrix dest;
    auto s = mesh->getAxis(0).size();
    dest.resize(s, s);
    dest.reserve(3 * s);
    state.ResumeTiming(); 
    dest = xpr.getStencil();
    benchmark::DoNotOptimize( dest );
  }
};

static void naive_scale(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    linops::NthPartialDeriv<1,0> u_x;
    linops::NthPartialDeriv<2,0> u_xx; 
    u_x.setMesh(mesh); 
    u_xx.setMesh(mesh);

    CSRMatrix dest;
    auto s = mesh->getAxis(0).size();
    dest.resize(s, s);
    dest.reserve(3 * s);
    state.ResumeTiming(); 
    dest = 2.0 * (u_xx.getStencil() + 0.5 * u_x.getStencil());
    benchmark::DoNotOptimize( dest );
  }
};

static void expression_scale(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    auto xpr = linops::NthPartialDeriv<2,0>{} + 0.5 * linops::NthPartialDeriv<1,0>{};
    xpr.setMesh(mesh);

    CSRMatrix dest;
    auto s = mesh->getAxis(0).size();
    dest.resize(s, s);
    dest.reserve(3 * s);
    state.ResumeTiming(); 
    dest = 2.0 * xpr.getStencil();
    benchmark::DoNotOptimize( dest );
  }
};

static void naive_apply(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    linops::NthPartialDeriv<1,0> u_x;
    linops::NthPartialDeriv<2,0> u_xx; 
    u_x.setMesh(mesh); 
    u_xx.setMesh(mesh);

    auto s = mesh->getAxis(0).size();
    Vector u = mesh->getAxis(0); 
    Vector v(s); 
    state.ResumeTiming(); 
    v = 2.0 * (u_xx.getStencil() + 0.5 * u_x.getStencil()) * u;
    benchmark::DoNotOptimize( v );
  }
};

static void expression_apply(benchmark::State &state)
{
  for (auto _ : state)
  {
    state.PauseTiming();
    int r = state.range(0);  
    auto mesh = make_Mesh(linspaced(r,0.0,double(r)), 1); 
    auto xpr = linops::NthPartialDeriv<2,0>{} + 0.5 * linops::NthPartialDeriv<1,0>{};
    xpr.setMesh(mesh);

    auto s = mesh->getAxis(0).size();
    Vector u = mesh->getAxis(0); 
    Vector v(s); 
    state.ResumeTiming(); 
    v = 2.0 * xpr.getStencil() * u;
    benchmark::DoNotOptimize( v );
  }
};

// linops 
BENCHMARK(setMesh_single<1>)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(setMesh_single<2>)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(setMesh_expresion)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(naive_assign)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(expression_assign)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(naive_scale)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(expression_scale)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(naive_apply)->Arg(20)->Arg(40)->Arg(80)->Arg(160);
BENCHMARK(expression_apply)->Arg(20)->Arg(40)->Arg(80)->Arg(160);

BENCHMARK_MAIN();