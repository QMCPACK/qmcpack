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
    
    



/**@file radfunc.cpp
 * @brief Implement a code to debug radial functors for Jastrow orbitals.
 *
 * Currently, the test uses pade and wm functors.
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>
#include "Utilities/OhmmsInfo.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "QMCWaveFunctions/Jastrow/WMFunctor.h"
#include "QMCWaveFunctions/Jastrow/DerivWMFunctor.h"
#include "QMCWaveFunctions/Jastrow/DerivPadeFunctors.h"
using namespace qmcplusplus;


template<typename FT, typename DFT>
struct RadFunctorTest
{
  FT&  u;
  DFT& ub;
  FT& up;
  FT& um;
  RadFunctorTest(FT& u_in, DFT& ub_in, FT& up_in, FT& um_in)
    : u(u_in), ub(ub_in), up(up_in), um(um_in) {}

  void run(double r, double delta)
  {
    double v=u.f(r);
    double dv=u.df(r);
    double ub_p=up.evaluate(r);
    double ub_m=um.evaluate(r);
    std::cout << "checking derivative functor ";
    std::cout <<  (ub_p-ub_m)/(2*delta)-ub.f(r) << std::endl;
    double du,d2udr2;
    std::cout << "checking Functor::evaulate ";
    double e=u.evaluate(r,du,d2udr2);
    double dv_p=u.f(r+delta);
    double dv_m=u.f(r-delta);
    std::cout << "  dudr = " << (dv_p-dv_m)/(2*delta)-du;
    std::cout << " d2udr2 = " << (dv_p+dv_m-2.0*v)/(delta*delta)-d2udr2 << std::endl;
    std::cout << "checking DFunctor::evaulate ";
    e=ub.evaluate(r,du,d2udr2);
    v=ub.f(r);
    dv=ub.df(r);
    dv_p=ub.f(r+delta);
    dv_m=ub.f(r-delta);
    std::cout << " dudr = " << (dv_p-dv_m)/(2*delta)-du;
    std::cout << " d2udr2 = " << (dv_p+dv_m-2.0*v)/(delta*delta)-d2udr2 << std::endl;
  }
};

int main(int argc, char** argv)
{
  //histograms<boost::mt19937>();
  //histograms<boost::lagged_fibonacci607>();
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo Welcome(argc,argv,OHMMS::Controller->rank());
  const double delta=0.001;
  double rc=5.0;
  double r=3.5;
  double b=3.0;
  if(argc>3)
  {
    b=atof(argv[1]);
    r=atof(argv[2]);
    rc=atof(argv[3]);
  }
  else
  {
    std::cout << "Using default values" << std::endl;
    std::cout << "Usage : radfunc B distance cutoff-distance " << std::endl;
  }
  std::cout << "rc= " << rc << " distance= " << r << std::endl;
  std::cout << "Printing differences: small numbers are good." << std::endl;
  //test WMFunctors
  {
    typedef WMFunctor<double> RadFunctor;
    typedef DWMDBFunctor<double> DerivRadFunctor;
    RadFunctor u(b,rc);
    RadFunctor up(b+delta,rc);
    RadFunctor um(b-delta,rc);
    DerivRadFunctor ub(b,rc);
    std::cout << std::endl << "Testing WM functors " << std::endl;
    RadFunctorTest<RadFunctor,DerivRadFunctor> test(u,ub,up,um);
    test.run(r,delta);
  }
  //test PadeFunctors
  {
    typedef PadeFunctor<double> RadFunctor;
    typedef DPadeDBFunctor<double> DerivRadFunctor;
    RadFunctor u(-0.5,b);
    RadFunctor up(-0.5,b+delta);
    RadFunctor um(-0.5,b-delta);
    DerivRadFunctor ub(-0.5,b);
    std::cout << std::endl << "Testing pade functors " << std::endl;
    RadFunctorTest<RadFunctor,DerivRadFunctor> test(u,ub,up,um);
    test.run(r,delta);
  }
  OHMMS::Controller->finalize();
  return 0;
}

