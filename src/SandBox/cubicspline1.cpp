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

/** Implements a screened Function \f$u[r]=(1-z(r/rcut))/(1+B*z(r/rcut)\f$
 *
 * Short-range functor introduced by Wagner and Mitas, cond-mat/0610088
 */
struct TestFunc1
{
  typedef double real_type;
  ///input B
  real_type B0;
  ///input Rcut
  real_type Rcut;
  ///1/Rcut
  real_type OneOverRc;
  ///id
  std::string ID;
  ///constructor
  explicit TestFunc1(real_type b, real_type rc=7.5)
  {
    reset(b,rc);
  }
  inline void reset()
  {
    OneOverRc=1.0/Rcut;
  }
  void reset(real_type b, real_type rc)
  {
    B0=b;
    Rcut=rc;
    reset();
  }

  inline real_type f(real_type r)
  {
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return (1-z)/(1+B0*z);
  }
  inline real_type df(real_type r)
  {
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return -(1+B0)/(1+B0*z)*(1+B0*z)*OneOverRc*12*x*(1-2.0*x+x*x);
  }
  inline real_type d2f(real_type r)
  {
    return 1.0;
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
  }

  inline double f(double x)
  {
    double res=0;
    for(int i=0; i<C.size(); i++)
      res += C[i]*F[i]->f(x);
    return res;
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
    double x0=0.0013;
    double xL=0.99*Xmax;
    while(x0<xL)
    {
      val=aorb.splint(x0, grad, lap);
      //val=aorb.splint(x0);
      val0=f(x0);
      grad0=df(x0);
      lap0=d2f(x0);
      fout << x0 << std::setw(20) << val0 << std::setw(20) << val
           <<std::setw(20) <<  (val-val0)
           //<<std::setw(20) <<  (grad-grad0)/grad0
           //<<std::setw(20) <<  (lap-lap0)/lap0 << std::endl;
           //<<std::setw(20) <<  grad0
           <<std::setw(20) <<  grad
           //<<std::setw(20) <<  lap0
           <<std::setw(20) <<  lap
           << std::endl;
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
  double rf = 5;
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
  infunc.push_back(-0.137958,new TestFunc1(4.55682,rf));
  infunc.push_back(-1.09405,new TestFunc1(1.83679,rf));
  infunc.push_back(1.06578,new TestFunc1(4.59201,rf));
  infunc.push_back(0.96602,new TestFunc1(0.226152,rf));
  std::vector<double> inData(npts);
  for(int ix=0; ix<npts; ix++)
  {
    inData[ix]=infunc.f(gridX(ix));
  }
  OneDimCubicSpline<double,double> aorb(&gridX,inData);
  //double yp0=(infunc.f(0.001)-infunc.f(0.0))*1000.;
  aorb.spline(0,infunc.df(ri),npts-1,0.0);
  infunc.compare(aorb,"cs.dat");
  bool closedEnd=true;
  CubicBspline<double,0,FIRSTDERIV_CONSTRAINTS> borb;
  borb.Init(ri,rf,inData,closedEnd,infunc.df(ri),0.0);
  infunc.compare(borb,"bs.dat");
  CubicSpline<double,0,FIRSTDERIV_CONSTRAINTS> corb;
  corb.Init(ri,rf,inData,closedEnd,infunc.df(ri),0.0);
  infunc.compare(corb,"csnew.dat");
  int niter=100000000;
  Timer myTimer;
  double sum=0.0,grad,lap;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += aorb.splint(Random(), grad, lap);
  }
  std::cout << "Cubic spline " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += aorb.splint(Random());
  }
  std::cout << "Cubic spline (only value) " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += borb.splint(Random());
  }
  std::cout << "Cubic Bspline (only value) " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += corb.splint(Random(), grad, lap);
  }
  std::cout << "New Cubic Spline " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += corb.splint(Random());
  }
  std::cout << "New Cubic Spline (only value) " << myTimer.elapsed() << " " << sum << std::endl;
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
    sum += infunc.f(Random());
  }
  std::cout << "Analytic function (only value) " << myTimer.elapsed() << " " << sum << std::endl;
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += infunc.df(Random());
  }
  std::cout << "Analytic function (only gradient) " << myTimer.elapsed() << " " << sum << std::endl;
  /*
  sum=0.0;
  myTimer.restart();
  for(int i=0; i<niter; i++)
  {
    sum += corb.splint(Random(), grad, lap);
  }
  std::cout << "New Cubic Spline " << myTimer.elapsed() << " " << sum << std::endl;
  */
  return 0;
}
