//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include <Configuration.h>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <spline/einspline_engine.hpp>
namespace qmcplusplus
{
  template<typename ENGT>
    struct einspline3d_benchmark
    {
      typedef typename einspline_engine<ENGT>::real_type real_type;
      typedef typename einspline_engine<ENGT>::value_type value_type;
      typedef TinyVector<real_type,3> pos_type;
      einspline_engine<ENGT> einspliner;
      pos_type start;
      pos_type end;
      vector<pos_type> pos;
      Vector<value_type> psi;
      Vector<value_type> lap;
      Vector<TinyVector<value_type,3> > grad;
      Vector<Tensor<value_type,3> > hess;

      einspline3d_benchmark():start(0.0),end(1.0){}

      einspline3d_benchmark(int nx, int ny, int nz, int num_splines)
        : start(0.0),end(1.0)
      {
        set(nx,ny,nz);
      }

      void set(int nx, int ny, int nz, int num_splines)
      {
        TinyVector<int,3> ng(nx,ny,nz);
        einspliner.create(start,end,ng,PERIODIC,num_splines);
        Array<value_type,3> data(nx,ny,nz);
        for(int i=0; i<num_splines; ++i)
        {
          for(int j=0; j<data.size();++j) data(j) = Random();
          einspliner.set(i,data);
        }
      }

      /** generate random positions and containers
       */
      void set_samples(int n)
      {
        pos.resize(n);
        RandomGenerator_t myrand;
        for(int i=0; i<n; ++i) pos[i]=pos_type(myrand(),myrand(),myrand());

        psi.resize(n);
        grad.resize(n);
        lap.resize(n);
        hess.resize(n);
      }

      inline TinyVector<double,3> test_all(int niters)
      {
        TinyVector<double,3> timers;
        Timer clock;
        for(int i=0; i<niters; ++i)
        {
          clock.restart();
          test_v();
          timers[0]+=clock.elapsed();
          clock.restart();
          test_vgl();
          timers[1]+=clock.elapsed();
          clock.restart();
          test_vgh();
          timers[2]+=clock.elapsed();
        }
        return timers;
      }

      inline void test_v()
      {
        for(int i=0; i<pos.size(); ++i) einspliner.evaluate(pos[i],psi);
      }

      inline void test_vgl()
      {
        for(int i=0; i<pos.size(); ++i) einspliner.evaluate(pos[i],psi,grad,lap);
      }

      inline void test_vgh()
      {
        for(int i=0; i<pos.size(); ++i) einspliner.evaluate(pos[i],psi,grad,hess);
      }

    };
}

int main(int argc, char** argv)
{
  using namespace qmcplusplus;

  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome(argc,argv,mycomm->rank());
  qmcplusplus::Random.init(0,1,11);

  int nx=73,ny=91,nz=29;
  int num_splines=21;
  int niters=100;
  int nsamples=32;

  if(mycomm->rank()==0)
  {
    cout << "#einspline benchmark grid = " << nx << " " << ny << " " << nz
      << " num_splines = " << num_splines << " num_samples = " << nsamples << " niters=" << niters << endl;
    cout << "#datatype     value      vgl      vgh " << endl;
  }

  typedef TinyVector<double,3> timer_type;
  timer_type d_timer_t,s_timer_t,z_timer_t,c_timer_t;

#pragma omp parallel
  {
    einspline3d_benchmark<multi_UBspline_3d_d> d_bench;
    d_bench.set(nx,ny,nz,num_splines);
    d_bench.set_samples(nsamples);
    timer_type d_timer=d_bench.test_all(niters);

    einspline3d_benchmark<multi_UBspline_3d_s> s_bench;
    s_bench.set(nx,ny,nz,num_splines);
    s_bench.set_samples(nsamples);
    timer_type s_timer=s_bench.test_all(niters);

    einspline3d_benchmark<multi_UBspline_3d_z> z_bench;
    z_bench.set(nx,ny,nz,num_splines);
    z_bench.set_samples(nsamples);
    timer_type z_timer=z_bench.test_all(niters);

    einspline3d_benchmark<multi_UBspline_3d_c> c_bench;
    c_bench.set(nx,ny,nz,num_splines);
    c_bench.set_samples(nsamples);
    timer_type c_timer=c_bench.test_all(niters);

#pragma omp critical
    {
      d_timer_t += d_timer;
      s_timer_t += s_timer;
      z_timer_t += z_timer;
      c_timer_t += c_timer;
    }
  }

  double nops=num_splines*nsamples*niters*omp_get_max_threads();
  cout.precision(6);
  cout.setf(std::ios::scientific, std::ios::floatfield);
  cout << "einspline::double          " << nops/d_timer_t << endl;
  cout << "einspline::float           " << nops/s_timer_t << endl;
  cout << "einspline::complex<double> " << nops/z_timer_t << endl;
  cout << "einspline::complex<float>  " << nops/c_timer_t << endl;

  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/
