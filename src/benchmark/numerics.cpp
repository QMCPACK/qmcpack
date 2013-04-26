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
#include <simd/simd.hpp>
using namespace qmcplusplus;

void compare_blas_1_2(int n, int m, int niters)
{
  typedef TinyVector<double,3> pos_type;
  Matrix<double> A(n,m);
  Vector<double> v(m),v_out(m), v2(m),v2_out(m);
  Vector<pos_type> p(m),p_out(m);
  Vector<TinyVector<double,4> > q(m),q_out(m);
  for(int j=0; j<n; ++j)
    v[j]=Random();
  for(int j=0; j<n; ++j)
    v2[j]=Random();
  for(int j=0; j<n; ++j)
    p[j]=pos_type(Random(),Random(),Random());
  Timer clock;
  for(int i=0; i<niters; ++i)
  {
    for(int j=0; j<n; ++j)
      v_out[j] =simd::dot(A[j],v.data(),m);
  }
  double dt_dot=clock.elapsed();
  clock.restart();
  for(int i=0; i<niters; ++i)
  {
    for(int j=0; j<n; ++j)
      v_out[j] =simd::dot(A[j],v.data(),m);
  }
  double dt_blas=clock.elapsed();
  clock.restart();
  for(int i=0; i<niters; ++i)
  {
    MatrixOperators::product(A,v,v_out.data());
  }
  double dt_gemm=clock.elapsed();
  clock.restart();
  for(int i=0; i<niters; ++i)
  {
    for(int j=0; j<n; ++j)
      v_out[j] =simd::dot(A[j],v.data(),m);
    for(int j=0; j<n; ++j)
      p_out[j] =simd::dot(A[j],p.data(),m);
  }
  double dt_v_blas1=clock.elapsed();
  clock.restart();
  for(int i=0; i<niters; ++i)
  {
    MatrixOperators::product(A,v.data(),v_out.data());
    MatrixOperators::product(A,p.data(),p_out.data());
  }
  double dt_v_blas2=clock.elapsed();
  clock.restart();
  for(int i=0; i<niters; ++i)
    MatrixOperators::product(A,q.data(),q_out.data());
  double dt_q_blas2=clock.elapsed();
  double f=1.0/static_cast<double>(niters);
  cout << n << " " << m << " " << dt_dot*f << " " << dt_gemm*f << " " << dt_blas*f << " "
       << " " << dt_v_blas1*f << " " << dt_v_blas2*f <<  " " << dt_q_blas2*f << endl;
}

inline int ops(int n, int m)
{
  return n*m;
}
int main(int argc, char** argv)
{
  Random.init(0,1,11);
  int niters=1<<25;
  int n=4;
  cout << "# n m dot_v gemm_v blas_v dot_g blase2_g blas2_gl " << endl;
  for(int m=n; m<600; m*= 2)
  {
    int i=niters/ops(m,m);
    compare_blas_1_2(m,m,i);
  }
  for(int m=n; m<600; m*= 2)
  {
    int i=niters/ops(m,m*8);
    compare_blas_1_2(m,m*8,i);
  }
  for(int m=n; m<600; m*= 2)
  {
    int i=niters/ops(m,m*11);
    compare_blas_1_2(m,m*11,i);
  }
  for(int m=n; m<600; m*= 2)
  {
    int i=niters/ops(m,m*13);
    compare_blas_1_2(m,m*13,i);
  }
  for(int m=n; m<600; m*= 2)
  {
    int i=niters/ops(m,m*16);
    compare_blas_1_2(m,m*16,i);
  }
  for(int m=n; m<600; m*= 2)
  {
    int i=niters/ops(m,m*32);
    compare_blas_1_2(m,m*32,i);
  }
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/
