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
    
    



/** @file einspline_validation.cpp
 *
 * Test code to check the correctness of the einspline
 */
#include <SandBox/SplineTestBase.h>
#include <Message/Communicate.h>
#include <Estimators/accumulators.h>
#include <mpi/collectives.h>
#include <getopt.h>
namespace qmcplusplus
{

template<typename ENGT>
void SplineTest<ENGT>
::test(std::vector<TinyVector<typename SplineTest<ENGT>::real_type,3> >& coord)
{
  typedef TinyVector<real_type,3> pos_type;
  typedef TinyVector<value_type,3> grad_type;
  typedef Tensor<value_type,3> hess_type;
  Vector<value_type> psi(num_splines);
  Vector<value_type> lap(num_splines);
  Vector<TinyVector<value_type,3> > grad(num_splines);
  Vector<Tensor<value_type,3> > hess(num_splines);
  for(int k=0; k<hess.size(); ++k)
    hess[k]=hess_type(0.0);
  accumulator_set<real_type> dv,dvgl_v,dvgl_g,dvgl_l,dvgh_v,dvgh_g,dvgh_h;
  accumulator_set<real_type> dv_copy;
  einspline3d_benchmark<ENGT> d_bench;
  d_bench.set(nx,ny,nz,num_splines,false);
  for(int i=0; i<num_splines; ++i)
    d_bench.assign(i,inData);
  for(int i=0; i<nsamples; ++i)
  {
    value_type v,l;
    grad_type g;
    hess_type h;
    //evaluate analytic solution
    infunc.v(coord[i],v);
    d_bench.einspliner.evaluate(coord[i],psi);
    dv(diff(v,psi[0]));
    infunc.vgl(coord[i],v,g,l);
    d_bench.einspliner.evaluate_vgl(coord[i],psi,grad,lap);
    dvgl_v(diff(v,psi[0]));
    dvgl_g(diff(g,grad[0]));
    dvgl_l(diff(l,lap[0]));
    //for(int k=0; k<hess.size(); ++k)
    //  hess[k]=hess_type(0.0);
    infunc.vgh(coord[i],v,g,h);
    d_bench.einspliner.evaluate_vgh(coord[i],psi,grad,hess);
    dvgh_v(diff(v,psi[0]));
    dvgh_g(diff(g,grad[0]));
    dvgh_h(diff(h,hess[0]));
    for(int k=1; k<hess.size(); ++k)
      if(diff(hess[0],hess[k]) > std::numeric_limits<real_type>::epsilon())
        app_log() << "Check multi evaluation = " << hess[0] << "\n"  << hess[k] << std::endl;
  }
  app_log() << "Average difference = " << dv.mean() << std::endl;
  app_log() << "Average difference VGL = " << dvgl_v.mean() << " " << dvgl_g.mean() << " " << dvgl_l.mean()  <<  std::endl;
  app_log() << "Average difference VGH = " << dvgh_v.mean() << " " << dvgh_g.mean() << " " << dvgh_h.mean()  <<  std::endl;
}
}

int main(int argc, char** argv)
{
  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("einspline",mycomm->rank());
  Random.init(0,1,-1);
  SplineTestBase param(argc,argv);
  typedef TinyVector<double,3> pos_type;
  std::vector<pos_type> coord(param.nsamples);
  for(int i=0; i<coord.size(); ++i)
    coord[i]=pos_type(Random(),Random(),Random()/2.0);
  std::vector<TinyVector<float,3> > coord_s(coord.size());
  for(int i=0; i<coord.size(); ++i)
    convert(coord[i],coord_s[i]);
  std::cout << "\n Testing double " << std::endl;
  SplineTest<multi_UBspline_3d_d> test_d(param);
  test_d.test(coord);
  //cout << "\n Testing std::complex<double> " << std::endl;
  //SplineTest<multi_UBspline_3d_z> test_z(param);
  //test_z.test(coord);
  //cout << "\n Testing float " << std::endl;
  //SplineTest<multi_UBspline_3d_s> test_s(param);
  //test_s.test(coord_s);
  //cout << "\n Testing std::complex<float> " << std::endl;
  //SplineTest<multi_UBspline_3d_c> test_c(param);
  //test_c.test(coord_s);
  return 0;
}
