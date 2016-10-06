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
    
    



#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "OhmmsData/OhmmsElementBase.h"
#include "Optimize/VarList.h"
#include "Numerics/TricubicBsplineSet.h"
#include "Numerics/OneDimGridBase.h"
#include "Utilities/Clock.h"
using namespace qmcplusplus;

struct TestFunc
{

  typedef std::complex<double> value_type;
  TinyVector<double,3> K;
  double d2factor;

  TestFunc(int nk0=1, int nk1=1, int nk2=1)
  {
    const double twopi = 8.0*std::atan(1.0);
    K=TinyVector<double,3>
      (twopi*static_cast<double>(nk0), twopi*static_cast<double>(nk1), twopi*static_cast<double>(nk2));
    d2factor = -dot(K,K);
  }
  //inline double operator()(const TinyVector<double,3>& pos) {
  //  return sin(k0*pos[0])*sin(k1*pos[1])*sin(k2*pos[2]);
  //}
  //inline double operator()(double x, double y, double z) {
  //  return sin(k0*x)*sin(k1*y)*sin(k2*z);
  //}
  //
  template<class PV>
  inline value_type f(const PV& pos)
  {
    double phi(dot(pos,K));
    return value_type(std::cos(phi),std::sin(phi));
  }

  template<class PV>
  inline TinyVector<value_type,3> df(const PV& pos)
  {
    double phi(dot(pos,K));
    double c=std::cos(phi);
    double s=std::sin(phi);
    return TinyVector<value_type,3>
           (value_type(-s*K[0],c*K[0]),
            value_type(-s*K[1],c*K[1]),
            value_type(-s*K[2],c*K[2]));
  }

  template<class PV>
  inline std::complex<double> d2f(const PV& pos)
  {
    return d2factor*f(pos);
  }

};

struct ComboFunc
{

  typedef std::complex<double> value_type;
  std::vector<value_type> C;
  std::vector<TestFunc*> F;

  ComboFunc() {}
  ~ComboFunc()
  {
    for(int i=0; i<F.size(); i++)
      delete F[i];
  }

  void push_back(value_type c, TestFunc* fn)
  {
    C.push_back(c);
    F.push_back(fn);
  }

  template<class PV>
  inline value_type f(const PV& pos)
  {
    value_type res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->f(pos);
    return res;
  }

  template<class PV>
  inline TinyVector<value_type,3> df(const PV& pos)
  {
    TinyVector<value_type,3> res(0.0,0.0,0.0);
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->df(pos);
    return res;
  }

  template<class PV>
  inline value_type d2f(const PV& pos)
  {
    value_type res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->d2f(pos);
    return res;
  }

};

int main(int argc, char** argv)
{
  double ri = 0.0;
  double rf = 1.0;
  std::vector<int> npts(3);
  npts[0]=101;
  npts[1]=101;
  npts[2]=101;
  double xcut=0.23;
  double ycut=0.67;
  const int nk0=1;
  const int nk1=1;
  const int nk2=1;
  typedef std::complex<double> value_type;
  typedef LinearGrid<double> GridType;
  GridType gridX, gridY, gridZ;
  gridX.set(ri,rf,npts[0]);
  gridY.set(ri,rf,npts[1]);
  gridZ.set(ri,rf,npts[2]);
  //Create an analytic function for assignment
  ComboFunc infunc;
  infunc.push_back(value_type(0.5,0.0),new TestFunc(1,1,1));
  infunc.push_back(value_type(0.3,0.0),new TestFunc(1,1,2));
  infunc.push_back(value_type(0.01,0.0),new TestFunc(5,3,2));
  infunc.push_back(value_type(0.01,0.0),new TestFunc(5,7,1));
  //infunc.push_back(0.1,new TestFunc(1,2,1));
  //infunc.push_back(0.01,new TestFunc(2,1,1));
  //infunc.push_back(0.01,new TestFunc(2,2,1));
  //infunc.push_back(0.001,new TestFunc(2,1,2));
  //infunc.push_back(0.001,new TestFunc(2,2,2));
  //infunc.push_back(0.001,new TestFunc(5,5,5));
  //infunc.push_back(-0.3,new TestFunc(7,2,3));
  //infunc.push_back(0.01,new TestFunc(7,7,7));
  //infunc.push_back(0.001,new TestFunc(5,5,5));
  //Write to an array
  Array<value_type,3> inData(npts[0],npts[1],npts[2]);
  value_type *it=inData.data();
  for(int ix=0; ix<npts[0]; ix++)
  {
    double x(gridX(ix));
    for(int iy=0; iy<npts[1]; iy++)
    {
      double y(gridY(iy));
      for(int iz=0; iz<npts[2]; iz++)
      {
        TinyVector<double,3> pos(x,y,gridZ(iz));
        (*it)=infunc.f(pos);
        ++it;
      }
    }
  }
  TricubicBspline<value_type> aorb;
  aorb.setGrid(ri,rf,ri,rf,ri,rf,npts[0]-1,npts[1]-1,npts[2]-1);
  aorb.Init(inData);
  value_type lap,val,lap0,val0,g;
  TinyVector<value_type,3> grad,grad0;
  double x0=xcut;
  double y0=ycut;
  double z=-0.5, dz=0.01;
  std::ofstream fout("bs.value.dat");
  std::ofstream dfout("bs.grad.dat");
  std::ofstream d2fout("bs.lap.dat");
  fout.setf(std::ios::scientific, std::ios::floatfield);
  dfout.setf(std::ios::scientific, std::ios::floatfield);
  d2fout.setf(std::ios::scientific, std::ios::floatfield);
  while(z<1.5)
  {
    TinyVector<double,3> pos(x0,y0,z);
    val=aorb.evaluate(pos,grad,lap);
    val0=infunc.f(pos);
    grad0=infunc.df(pos);
    lap0=infunc.d2f(pos);
    fout << z
         << std::setw(20) << val0.real() << std::setw(20) << val.real()
         << std::setw(20) << val0.imag() << std::setw(20) << val.imag()
         << std::setw(20) << val.real()-val0.real() << std::setw(20) << val.imag()-val0.imag()
         << std::endl;
    dfout << z
          << std::setw(20) << grad0[0].real() << std::setw(20) << grad[0].real()
          << std::setw(20) << grad0[0].imag() << std::setw(20) << grad[0].imag()
          << std::setw(20) << grad0[1].real() << std::setw(20) << grad[1].real()
          << std::setw(20) << grad0[1].imag() << std::setw(20) << grad[1].imag()
          << std::setw(20) << grad0[2].real() << std::setw(20) << grad[2].real()
          << std::setw(20) << grad0[2].imag() << std::setw(20) << grad[2].imag()
          << std::endl;
    d2fout << z
           << std::setw(20) << lap0.real() << std::setw(20) <<  lap.real()
           << std::setw(20) << lap0.imag() << std::setw(20) <<  lap.imag()
           << std::setw(20) << (lap.real()-lap0.real())
           << std::setw(20) << (lap.imag()-lap0.imag())
           << std::endl;
    z+=dz;
  }
  return 0;
}
