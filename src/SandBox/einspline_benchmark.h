//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file einspline3d_benchmark.h
 * @brief define einspline3d_benchmark and random_position_generator
 */
#ifndef QMCPLUSPLUS_EINSPLINE3D_BENCHMARK_H
#define QMCPLUSPLUS_EINSPLINE3D_BENCHMARK_H
#include <Configuration.h>
#include <Utilities/RandomGenerator.h>
#include <Utilities/Timer.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <simd/simd.hpp>
#include <spline/einspline_engine.hpp>
//#include <SandBox/einspline_util.h>
//#include <einspline/multi_bspline_copy.h>
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
  Vector<value_type> psi;
  Vector<value_type> lap;
  Vector<TinyVector<value_type,3> > grad;
  Vector<Tensor<value_type,3> > hess;

  einspline3d_benchmark():start(0.0),end(1.0) { }

  einspline3d_benchmark(int nx, int ny, int nz, int num_splines)
    : start(0.0),end(1.0)
  {
    set(nx,ny,nz);
  }

  einspline3d_benchmark(einspline3d_benchmark<ENGT>& in)
  {
    std::cout << "Calling copy constructor " << std::endl;
    this->einspliner.spliner=copy_multi_UBspline_3d_d(in.einspliner.spliner);
    //this->einspliner.spliner=(ENGT*)(malloc(sizeof(ENGT)));
    //this->einspliner.spliner=new ENGT;
    //copy_einspline(in.einspliner.spliner, this->einspliner.spliner);
  }

  /** return the spline engine */
  ENGT* spliner()
  {
    return einspliner.spliner;
  }

  template<typename VT>
  void assign(int i, VT& data)
  {
    einspliner.set(i,data);
  }

  void set(int nx, int ny, int nz, int num_splines, bool initialize=true)
  {
    TinyVector<int,3> ng(nx,ny,nz);
    einspliner.create(start,end,ng,PERIODIC,num_splines);
    if(initialize)
    {
      Array<value_type,3> data(nx,ny,nz);
      for(int i=0; i<num_splines; ++i)
      {
        for(int j=0; j<data.size(); ++j)
          data(j) = Random();
        einspliner.set(i,data);
      }
    }
    psi.resize(num_splines);
    grad.resize(num_splines);
    lap.resize(num_splines);
    hess.resize(num_splines);
  }

  inline TinyVector<double,3> test_all(const std::vector<pos_type>& Vpos,
                                       const std::vector<pos_type>& VGLpos, const std::vector<pos_type>& VGHpos)
  {
    TinyVector<double,3> timers;
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

  inline void test_v(const std::vector<pos_type>& coord) //const
  {
    for(int i=0; i<coord.size(); ++i)
      einspliner.evaluate(coord[i],psi);
  }

  inline void test_vgl(const std::vector<pos_type>& coord)// const
  {
    for(int i=0; i<coord.size(); ++i)
      einspliner.evaluate_vgl(coord[i],psi,grad,lap);
  }

  inline void test_vgh(const std::vector<pos_type>& coord) //const
  {
    for(int i=0; i<coord.size(); ++i)
      einspliner.evaluate_vgh(coord[i],psi,grad,hess);
  }

};

template<typename T>
struct random_position_generator
{
  typedef TinyVector<T,3> pos_type;
  typedef RandomGenerator_t::uint_type uint_type;
  RandomGenerator_t myrand;
  std::vector<pos_type> Vpos, VGLpos, VGHpos;
  random_position_generator(int n,  uint_type seed):myrand(seed)
  {
    Vpos.resize(n);
    VGLpos.resize(n);
    VGHpos.resize(n);
  }
  inline void randomize()
  {
    for(int i=0; i<Vpos.size(); ++i)
      Vpos[i]=pos_type(myrand(),myrand(),myrand());
    for(int i=0; i<VGLpos.size(); ++i)
      VGLpos[i]=pos_type(myrand(),myrand(),myrand());
    for(int i=0; i<VGHpos.size(); ++i)
      VGHpos[i]=pos_type(myrand(),myrand(),myrand());
  }
};
}
#endif
