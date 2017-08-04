//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file einspline3d_benchmark.h
 * @brief define einspline3d_benchmark and random_position_generator
 */
#ifndef QMCPLUSPLUS_EINSPLINE3D_BENCHMARK_H
#define QMCPLUSPLUS_EINSPLINE3D_BENCHMARK_H
#include <omp.h>
#include <random/random.hpp>
#include <Utilities/Timer.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <simd/allocator.hpp>
//#include <spline2/einspline_engine.hpp>
#include <spline2/MultiBspline.hpp>

namespace qmcplusplus
{

template<typename T>
struct einspline3d_benchmark
{
#if (__cplusplus< 201103L)
  typedef MultiBspline<T> engine_type;
  typedef typename engine_type::real_type real_type;
  typedef TinyVector<real_type,3> pos_type;
  static const int nullptr=0;
#else
  using engine_type=MultiBspline<T>;
  using real_type=typename engine_type::real_type;
  using pos_type=TinyVector<real_type,3>;
#endif

  bool own_engine;
  engine_type* einspliner;

  pos_type start;
  pos_type end;
#if (__cplusplus< 201103L)
  vector<real_type> psi;
  vector<real_type> lap;
  vector<real_type> grad;
  vector<real_type> hess;
#else
  aligned_vector<real_type> psi;
  aligned_vector<real_type> lap;
  aligned_vector<real_type> grad;
  aligned_vector<real_type> hess;
#endif

  einspline3d_benchmark():own_engine(false),einspliner(nullptr),start(0.0),end(1.0) { }

  ~einspline3d_benchmark()
  {
    if(own_engine && einspliner!=nullptr) delete einspliner;
  }

  /** copy constructor to share einspliner and everything else is generated */
  einspline3d_benchmark(einspline3d_benchmark<T>& in): start(in.start), end(in.end)
  {
    own_engine=false;
    einspliner=in.einspliner;
    int num_splines=einspliner->num_splines();
    psi.resize(num_splines);
    lap.resize(num_splines*3);
    grad.resize(num_splines*3);
    hess.resize(num_splines*6);
  }


  template<typename VT>
  void assign(int i, VT& data)
  {
    einspliner->set(i,data);
  }

  void set(int nx, int ny, int nz, int num_splines, bool initialize=true)
  {
    TinyVector<int,3> ng(nx,ny,nz);
    if(einspliner == nullptr)
    {
      own_engine=true;
      einspliner=new engine_type;
      einspliner->create(start,end,ng,PERIODIC,num_splines);
      if(initialize)
      {
#if (__cplusplus< 201103L)
        BoostRandom<real_type> rng; 
#else
        RandomGenerator<real_type> rng; 
#endif
        Array<real_type,3> data(nx,ny,nz);
        rng.generate_uniform(data.data(),data.size());
        for(int i=0; i<num_splines; ++i)
          einspliner->set(i,data);
      }
    }
    psi.resize(num_splines);
    grad.resize(num_splines*3);
    lap.resize(num_splines*3);
    hess.resize(num_splines*6);
  }

  inline TinyVector<double,3> test_all(const vector<pos_type>& Vpos,
                                       const vector<pos_type>& VGLpos, const vector<pos_type>& VGHpos)
  {
    TinyVector<double,3> timers(0.0);
    Timer clock;
    test_v(Vpos);
    timers[0]+=clock.elapsed();
    clock.restart();
    test_vgl(VGLpos);
    timers[1]+=clock.elapsed();
    clock.restart();
    test_vgh(VGHpos);
    timers[2]+=clock.elapsed();
    return timers;
  }

  inline void test_v(const vector<pos_type>& coord) //const
  {
    for(int i=0; i<coord.size(); ++i)
      einspliner->evaluate(coord[i],psi);
  }

  inline void test_vgl(const vector<pos_type>& coord)// const
  {
    for(int i=0; i<coord.size(); ++i)
      einspliner->evaluate_vgl(coord[i],psi,grad,lap);
  }

  inline void test_vgh(const vector<pos_type>& coord) //const
  {
    for(int i=0; i<coord.size(); ++i)
      einspliner->evaluate_vgh(coord[i],psi,grad,hess);
  }

};

template<typename T>
struct random_position_generator
{
  typedef TinyVector<T,3> pos_type;
#if (__cplusplus< 201103L)
  typedef uint32_t uint_type;
  BoostRandom<T> myrand;
#else
  typedef typename RandomGenerator<T>::uint_type uint_type;
  RandomGenerator<T> myrand;
#endif
  vector<pos_type> Vpos, VGLpos, VGHpos;
  random_position_generator(int n,  uint_type seed):myrand(seed)
  {
    Vpos.resize(n);
    VGLpos.resize(n);
    VGHpos.resize(n);
  }
  inline void randomize()
  {
    int n3=Vpos.size()*3;
    myrand.generate_uniform(&Vpos[0][0],n3);
    myrand.generate_uniform(&VGLpos[0][0],n3);
    myrand.generate_uniform(&VGHpos[0][0],n3);
    //for(int i=0; i<Vpos.size(); ++i)
    //  Vpos[i]=pos_type(myrand(),myrand(),myrand());
    //for(int i=0; i<VGLpos.size(); ++i)
    //  VGLpos[i]=pos_type(myrand(),myrand(),myrand());
    //for(int i=0; i<VGHpos.size(); ++i)
    //  VGHpos[i]=pos_type(myrand(),myrand(),myrand());
  }
  inline void randomize_v()
  {
    myrand.generate_uniform(&Vpos[0][0],Vpos.size()*3);
  }
};
}
#endif

