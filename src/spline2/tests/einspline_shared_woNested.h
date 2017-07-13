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
#include <spline2/MultiBspline.hpp>
using namespace std;

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

  bool bInitialized;
  int numSplinesPerThread;
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

  einspline3d_benchmark():
    own_engine(false), einspliner(nullptr), bInitialized(false), numSplinesPerThread(0),start(0.0),end(1.0)
  { }

  /** copy constructor to share einspliner and everything else is generated */
  einspline3d_benchmark(const einspline3d_benchmark<T>& in): 
    start(in.start), end(in.end), bInitialized(in.bInitialized)
  {
    own_engine=false;
    this->numSplinesPerThread = in.numSplinesPerThread;
    this->einspliner=in.einspliner;

    this->psi.resize(numSplinesPerThread);
    this->lap.resize(numSplinesPerThread*3);
    this->grad.resize(numSplinesPerThread*3);
    this->hess.resize(numSplinesPerThread*6);
  }

  /** copy constructor to share einspliner and everything else is generated */
  einspline3d_benchmark<T>& operator = (const einspline3d_benchmark<T> in)
  {
    this->numSplinesPerThread = in.numSplinesPerThread;
    this->einspliner=in.einspliner;

    this->psi.resize(numSplinesPerThread);
    this->lap.resize(numSplinesPerThread*3);
    this->grad.resize(numSplinesPerThread*3);
    this->hess.resize(numSplinesPerThread*6);

    return *this;
  }

  void set(int nx, int ny, int nz, int num_splines, int nTeamThreads)
  {
    TinyVector<int,3> ng(nx,ny,nz);
    if(einspliner == nullptr)
    {
      own_engine=true;
      einspliner=new engine_type;
      einspliner->create(start,end,ng,PERIODIC,num_splines,nTeamThreads);
      Array<real_type,3> data(nx,ny,nz);
#if (__cplusplus< 201103L)
      BoostRandom<real_type> rng; 
#else
      RandomGenerator<real_type> rng; 
#endif
      rng.generate_uniform(data.data(),data.size());
      for(int i=0; i<num_splines; ++i)
        einspliner->set(i,data);
    }
    numSplinesPerThread = num_splines;
    psi.resize(numSplinesPerThread);
    grad.resize(numSplinesPerThread*3);
    lap.resize(numSplinesPerThread*3);
    hess.resize(numSplinesPerThread*6);
    bInitialized = true;
  }

  inline TinyVector<double,3> test_all(const vector<pos_type>& Vpos,
                                       const vector<pos_type>& VGLpos, const vector<pos_type>& VGHpos,
                                       const int iIDInTeam )
  {
    TinyVector<double,3> timers(0.0);
    Timer clock;
    test_v_omp(Vpos,iIDInTeam);
    timers[0]+=clock.elapsed();
    clock.restart();
    test_vgl_omp(VGLpos,iIDInTeam);
    timers[1]+=clock.elapsed();
    clock.restart();
    test_vgh_omp(VGHpos,iIDInTeam);
    timers[2]+=clock.elapsed();
    return timers;
  }

  inline void test_v_omp(const vector<pos_type>& coord, const int iIDInTeam) //const
  {
      for(int i=0; i<coord.size(); ++i) {
        einspliner->evaluate(coord[i],psi, iIDInTeam );
      }
  }

  inline void test_vgl_omp(const vector<pos_type>& coord, const int iIDInTeam)// const
  {
      for(int i=0; i<coord.size(); ++i)  {
        einspliner->evaluate_vgl(coord[i],psi,grad,lap, iIDInTeam );
      }
  }

  inline void test_vgh_omp(const vector<pos_type>& coord, const int iIDInTeam) //const
  {
      for(int i=0; i<coord.size(); ++i) {
        einspliner->evaluate_vgh(coord[i],psi,grad,hess, iIDInTeam );
      }
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
  }
};
}
#endif

