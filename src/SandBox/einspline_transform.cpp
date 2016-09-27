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
#include <spline/einspline_util.hpp>
#include <Message/Communicate.h>
#include "Numerics/HDFNumericAttrib.h"
#include <io/hdf_archive.h>
#include <mpi/collectives.h>
#include <getopt.h>

//testing subgrid
qmcplusplus::TinyVector<double,3> lower(0.25,0.25,0.25);
qmcplusplus::TinyVector<double,3> upper(0.75,0.75,0.75);

namespace qmcplusplus
{
template<typename PT>
inline bool inbox(const PT& pos)
{
  return (pos[0]>lower[0] && pos[0]<upper[0]
          &&  pos[1]>lower[1] && pos[1]<upper[0]
          &&  pos[2]>lower[2] && pos[2]<upper[0]
         );
}

template<typename ENGT>
void SplineTest<ENGT>
::test(std::vector<TinyVector<typename SplineTest<ENGT>::real_type,3> >& coord)
{
  //read data from a file
  int ngx=160;
  int ngy=160;
  int ngz=160;
  Array<double,3> big(ngx,ngy,ngz),small(ngx/2,ngy/2,ngz/2);
  {
    hdf_archive h5f;
    h5f.open("sample.h5");
    h5f.read(big,"big_40");
    h5f.read(small,"small_40");
  }
  typedef TinyVector<real_type,3> pos_type;
  pos_type start(0.0);
  pos_type end(1.0);
  TinyVector<int,3> ng(ngx,ngy,ngz), ngh(ngx/2,ngy/2,ngz/2);
  //create/set UBspline_3d_X with the dense grid
  typename bspline_engine_traits<ENGT>::SingleSplineType* dense=0;
  dense=einspline::create(dense,start,end,ng,PERIODIC);
  einspline::set(dense,big.data());
  //create/set UBspline_3d_X with the coarse grid
  typename bspline_engine_traits<ENGT>::SingleSplineType* coarse=0;
  coarse=einspline::create(coarse,start,end,ngh,PERIODIC);
  einspline::set(coarse,small.data());
  //create/set multi_UBspline_3d as usual on the dense grid
  einspline_engine<ENGT> einspliner;
  einspliner.create(start,end,ng,PERIODIC,num_splines);
  for(int i=0; i<num_splines; ++i)
    einspliner.set(i,big);
  //create multi_UBspline_3d as usual on a truncated grid based on the dense grid
  einspline_engine<multi_UBspline_3d_d> newspline;
  GridConvert<real_type> gconv;
  gconv.create(newspline.spliner,dense,lower,upper,num_splines);
  for(int i=0; i<num_splines; ++i)
    einspline::set(newspline.spliner,i,dense,gconv.Offset,gconv.N);
  std::cout << "Original dense grid " << std::endl;
  print_spliner(einspliner.spliner);
  std::cout << "Truncated dense grid " << std::endl;
  print_spliner(newspline.spliner);
  Vector<value_type> psi(num_splines);
  Vector<real_type> psi_c(num_splines);
  accumulator_set<real_type> dv,ds,din,dout;
  int n_in=0;
  int n_out=0;
  for(int i=0; i<nsamples; ++i)
  {
    value_type v_dense, v_mixed;
    einspliner.evaluate(coord[i],psi);
    v_dense = einspline::evaluate(dense,coord[i]);
    ds(diff(v_dense,psi[0]));
    if(inbox(coord[i]))
    {
      n_in++;
      newspline.evaluate(coord[i],psi_c);
      v_mixed=psi_c[0];
      din(diff(v_mixed,v_dense));
    }
    else
    {
      n_out++;
      v_mixed = einspline::evaluate(coarse,coord[i]);
      dout(diff(v_mixed,v_dense));
    }
    dv(diff(v_dense,v_mixed));
    if(verbose)
    {
      std::cout << "At " << coord[i] << std::endl;
      std::cout << "Dense Grid " << v_dense
           <<  "\nMixed Grid " << v_mixed
           << std::endl << std::endl;
    }
  }
  app_log() << "Number of coords in the box       " << n_in << std::endl;
  app_log() << "Number of coords outside the box  " << n_out << std::endl;
  app_log() << "Average difference org-mixed grid " << dv.mean()
            << "\n\tInside  " << din.mean()
            << "\n\tOutside " << dout.mean()
            << "\n\tmulti-spline:single-spline ="  << ds.mean()
            << std::endl;
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
  std::vector<TinyVector<double,3> > coord_s(coord.size());
  for(int i=0; i<coord.size(); ++i)
    convert(coord[i],coord_s[i]);
  std::cout << "\nTesting einspline transformation" << std::endl;
  SplineTest<multi_UBspline_3d_d> test_d(param);
  test_d.test(coord);
  return 0;
}
