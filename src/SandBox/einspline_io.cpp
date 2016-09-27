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
#include <einspline/multi_bspline_copy.h>
#include <Message/Communicate.h>
#include "Numerics/HDFNumericAttrib.h"
#include <io/hdf_archive.h>
#include <mpi/collectives.h>
#include <spline/einspline_util.hpp>
#include <getopt.h>

//testing subgrid
qmcplusplus::TinyVector<double,3> lower(0.0,0.0,0.0);
qmcplusplus::TinyVector<double,3> upper(1.0,1.0,1.0);

namespace qmcplusplus
{
template<typename ENGT>
void SplineTest<ENGT>
::test(std::vector<TinyVector<typename SplineTest<ENGT>::real_type,3> >& coord)
{
  //read data from a file
  int ngx=160;
  int ngy=160;
  int ngz=160;
  Array<double,3> big(ngx,ngy,ngz);
  Array<float,3> bigf(ngx,ngy,ngz);
  {
    hdf_archive h5f;
    h5f.open("sample.h5");
    h5f.read(big,"big_40");
  }
  //convert double to float
  simd::copy(bigf.data(), big.data(), big.size());
  typedef TinyVector<real_type,3> pos_type;
  pos_type start(0.0);
  pos_type end(1.0);
  TinyVector<int,3> ng(ngx,ngy,ngz);
  //create/set multi_UBspline_3d as usual on the dense grid
  einspline_engine<ENGT> einspliner;
  einspliner.create(start,end,ng,PERIODIC,num_splines);
  Timer clock;
  clock.restart();
  for(int i=0; i<num_splines; ++i)
    einspliner.set(i,bigf);
  double t_spline=clock.elapsed();
//#pragma omp parallel for
//      for(int i=0; i<num_splines; ++i)
//      {
//        einspliner.set(i,big);
//      }
  einspline_engine<ENGT> einspliner2;
  einspliner2.create(start,end,ng,PERIODIC,num_splines);
  //dump the spline to a hdf5
  clock.restart();
  {
    hdf_archive h5f;
    h5f.create("spline.h5");
    h5f.write(einspliner,"bspline");
  }
  double t_write=clock.elapsed();
  clock.restart();
  //initialize einspline2 from the file
  {
    hdf_archive h5f;
    h5f.open("spline.h5");
    h5f.read(einspliner2,"bspline");
  }
  double t_read=clock.elapsed();
  Vector<value_type> psi(num_splines);
  Vector<real_type> psi_c(num_splines);
  accumulator_set<real_type> dv;
  for(int i=0; i<nsamples; ++i)
  {
    value_type v_dense, v_mixed;
    einspliner.evaluate(coord[i],psi);
    einspliner2.evaluate(coord[i],psi_c);
    for(int k=0; k<num_splines; ++k)
      dv(diff(psi[k],psi_c[k]));
  }
  std::cout << "Average diff " << dv.mean() << std::endl;
  app_log() << "Timing spline=" << t_spline << " write = " << t_write << " read = " << t_read
            << "  saving = " << t_spline/t_read<< std::endl;
}
}
int main(int argc, char** argv)
{
  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("einspline",mycomm->rank());
  Random.init(0,1,11);
  SplineTestBase param(argc,argv);
  typedef TinyVector<double,3> pos_type;
  std::vector<pos_type> coord(param.nsamples);
  for(int i=0; i<coord.size(); ++i)
    coord[i]=pos_type(Random(),Random(),Random());
//  for(int i=0; i<coord.size(); ++i)
//    coord[i]=pos_type((upper[0]-lower[0])*Random()+lower[0]
//        , (upper[1]-lower[1])*Random()+lower[1]
//        , (upper[2]-lower[2])*Random()+lower[2]
//        );
//
  std::vector<TinyVector<float,3> > coord_s(coord.size());
  for(int i=0; i<coord.size(); ++i)
    convert(coord[i],coord_s[i]);
  std::cout << "\nTesting einspline transformation" << std::endl;
  SplineTest<multi_UBspline_3d_s> test_d(param);
  test_d.test(coord_s);
  return 0;
}
