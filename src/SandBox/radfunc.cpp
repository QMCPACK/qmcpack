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
using namespace std;
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
    cout << "checking derivative functor ";
    cout <<  (ub_p-ub_m)/(2*delta)-ub.f(r) << endl;
    double du,d2udr2;
    cout << "checking Functor::evaulate ";
    double e=u.evaluate(r,du,d2udr2);
    double dv_p=u.f(r+delta);
    double dv_m=u.f(r-delta);
    cout << "  dudr = " << (dv_p-dv_m)/(2*delta)-du;
    cout << " d2udr2 = " << (dv_p+dv_m-2.0*v)/(delta*delta)-d2udr2 << endl;
    cout << "checking DFunctor::evaulate ";
    e=ub.evaluate(r,du,d2udr2);
    v=ub.f(r);
    dv=ub.df(r);
    dv_p=ub.f(r+delta);
    dv_m=ub.f(r-delta);
    cout << " dudr = " << (dv_p-dv_m)/(2*delta)-du;
    cout << " d2udr2 = " << (dv_p+dv_m-2.0*v)/(delta*delta)-d2udr2 << endl;
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
    cout << "Using default values" << endl;
    cout << "Usage : radfunc B distance cutoff-distance " << endl;
  }
  cout << "rc= " << rc << " distance= " << r << endl;
  cout << "Printing differences: small numbers are good." << endl;
  //test WMFunctors
  {
    typedef WMFunctor<double> RadFunctor;
    typedef DWMDBFunctor<double> DerivRadFunctor;
    RadFunctor u(b,rc);
    RadFunctor up(b+delta,rc);
    RadFunctor um(b-delta,rc);
    DerivRadFunctor ub(b,rc);
    cout << endl << "Testing WM functors " << endl;
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
    cout << endl << "Testing pade functors " << endl;
    RadFunctorTest<RadFunctor,DerivRadFunctor> test(u,ub,up,um);
    test.run(r,delta);
  }
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/
