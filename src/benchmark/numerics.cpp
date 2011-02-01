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
/**@file numeric.cpp
 * @brief Test for matrix operations
 */
#include <Configuration.h>
#include "Numerics/DeterminantOperators.h"
#include <Numerics/MatrixOperators.h>
#include <Numerics/OhmmsBlas.h>
#include <Utilities/RandomGenerator.h>
#include <Utilities/Timer.h>
using namespace qmcplusplus;

void test_numerics(int n, int niters)
{
  typedef TinyVector<double,3> pos_type;
  Matrix<double> A(n,n);
  Vector<double> v(n),v_out(n);
  Vector<pos_type> p(n),p_out(n);

  for(int j=0; j<n; ++j) v[j]=Random();
  for(int j=0; j<n; ++j) p[j]=pos_type(Random(),Random(),Random());

  Timer clock;
  for(int i=0; i<niters; ++i)
  {
    for(int j=0; j<n; ++j) v_out[j]=dot(A[j],v.data(),n);
    for(int j=0; j<n; ++j) p_out[j]=dot(A[j],p.data(),n);
  }
  double dt_dot=clock.elapsed();

  clock.restart();
  for(int i=0; i<niters; ++i)
  {
    MatrixOperators::product(A,v,v_out.data());
    MatrixOperators::product(A,p.data(),p_out.data());
  }
  double dt_gemm=clock.elapsed();

  cout << " n= " << n << " dot/gemm = " << dt_dot/dt_gemm
    << " dot = " << dt_dot/static_cast<double>(niters) << " gemm = " << dt_gemm/static_cast<double>(niters) << endl;
}

int main(int argc, char** argv)
{
  Random.init(0,1,11);

  int niters=1<<20;
  int n=4;

  for(int m=n; m<512; m+= 8)
  {
    int i=niters/n/n;
    test_numerics(m,i);
  }

  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/
