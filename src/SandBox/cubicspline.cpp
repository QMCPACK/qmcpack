//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/CubicBspline.h"
#include "Numerics/CubicSpline.h"
#include "Numerics/CubicSplineEngine.h"
#include "Configuration.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommunicateGroup.h"

struct TestFunc1
{
  double k0;
  double minusksq;
  TestFunc1(int nk0=1, double L=1.0)
  {
    const double twopi = 8.0*std::atan(1.0);
    k0=twopi*static_cast<double>(nk0)/L;
    minusksq = -k0*k0;
  }
  //inline double f(double x) { return std::cos(k0*x); }
  //inline double df(double x) { return -k0*std::sin(k0*x); }
  //inline double d2f(double x) { return minusksq*std::cos(k0*x); }
  inline double f(double x)
  {
    return std::sin(k0*x);
  }
  inline double df(double x)
  {
    return k0*std::cos(k0*x);
  }
  inline double d2f(double x)
  {
    return minusksq*std::sin(k0*x);
  }
};

struct ComboFunc1
{

  std::vector<double> C;
  std::vector<TestFunc1*> F;
  double Y0;
  double Xmax;

  ComboFunc1(double dl):Xmax(dl),Y0(0.0) {}
  ~ComboFunc1()
  {
    for(int i=0; i<F.size(); i++)
      delete F[i];
  }

  void push_back(double c, TestFunc1* fn)
  {
    C.push_back(c);
    F.push_back(fn);
    Y0+= c*fn->f(0.0);
  }

  inline double f(double x)
  {
    double res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->f(x);
    return res-Y0;
  }

  inline double df(double x)
  {
    double res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->df(x);
    return res;
  }

  inline double d2f(double x)
  {
    double res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->d2f(x);
    return res;
  }

  template<class NFunc>
  void  compare(NFunc& aorb, const std::string& fname)
  {
    std::ofstream fout(fname.c_str());
    fout.setf(std::ios::scientific, std::ios::floatfield);
    fout.precision(12);
    double lap,val,lap0,val0,grad,grad0;
    double dx=0.01;
    double x0=-0.107*Xmax;
    double xL=1.307*Xmax;
    while(x0<xL)
    {
      val=aorb.splint(x0, grad, lap);
      val0=f(x0);
      grad0=df(x0);
      lap0=d2f(x0);
      fout << x0 << std::setw(20) << val0 << std::setw(20) << val
           <<std::setw(20) <<  (val-val0)
           <<std::setw(20) <<  (grad-grad0)/grad0
           <<std::setw(20) <<  (lap-lap0)/lap0 << std::endl;
      //val=aorb.splint(x0);
      //val0=f(x0);
      //fout << x0 << std::setw(20) << val0 << std::setw(20) << val
      //  <<std::setw(20) <<  (val-val0)
      //  << std::endl;
      x0+=dx;
    }
  }
};

using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv);
  Random.init(0,1,-1);
  double ri = 0.0;
  double rf = 3.145;
  double xcut=0.23;
  const int nk0=1;
  int npts=101;
  if(argc>1)
  {
    npts=atoi(argv[1]);
  }
  std::cout << "Ri= " << ri << " Rf= " << rf << "  " << npts << std::endl;
  typedef LinearGrid<double> GridType;
  GridType gridX;
  gridX.set(ri,rf,npts);
  double dL=rf-ri;
  //Create an analytic function for assignment
  ComboFunc1 infunc(dL);
  infunc.push_back(0.5,new TestFunc1(1,dL));
  infunc.push_back(-0.3,new TestFunc1(3,dL));
  infunc.push_back(0.01,new TestFunc1(5,dL));
  infunc.push_back(-0.01,new TestFunc1(9,dL));
  std::vector<double> inData(npts);
  for(int ix=0; ix<npts; ix++)
  {
    inData[ix]=infunc.f(gridX(ix));
  }
  OneDimCubicSplinePBC<double,double> aorb(&gridX,inData);
  aorb.spline();
  infunc.compare(aorb,"cs.dat");
  bool closedEnd=true;
  CubicBspline<double,0,PBC_CONSTRAINTS> borb;
  borb.Init(ri,rf,inData,closedEnd);
  infunc.compare(borb,"bs.dat");
  CubicSpline<double,0,PBC_CONSTRAINTS> corb;
  corb.Init(ri,rf,inData,closedEnd);
  infunc.compare(corb,"csnew.dat");
  int niter=100000000;
  Timer myTimer;
  myTimer.restart();
  double sum=0.0,grad,lap;
  for(int i=0; i<niter; i++)
  {
    sum += aorb.splint(Random(), grad, lap);
  }
  std::cout << "Cubic spline " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += borb.splint(Random(), grad, lap);
  }
  std::cout << "Cubic Bspline " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += corb.splint(Random(), grad, lap);
  }
  std::cout << "New Cubic Spline " << myTimer.elapsed() << " " << sum << std::endl;
  return 0;
}
