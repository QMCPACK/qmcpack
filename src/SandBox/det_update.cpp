//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Jeongnim Kim
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
/**@file det_update.cpp
 * @brief Test code for multideterminant tree
 */
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Numerics/MatrixOperators.h"
#include "SandBox/determiant_operators.h"
#include <QMCWaveFunctions/Fermion/excitation_node.h>
#include <QMCWaveFunctions/Fermion/ci_builder.h>
#include "Utilities/IteratorUtility.h"
using namespace qmcplusplus;
using namespace std;

extern "C"
{
  void dger_(const int& m, const int& n, const double& alpha
             , const double* x, const int& incx, const double* y, const int& incy
             , double* a, const int& lda);
}

template<typename T>
inline void det_ger(T* restrict pinv,  const T* restrict tv, int m, int rowchanged, T c_ratio)
{
  T ratio_inv=1.0/c_ratio;
  double temp[m], rcopy[m];
  dgemv('T', m, m, ratio_inv, pinv, m, tv, 1, T(), temp, 1);
  temp[rowchanged]=1.0-ratio_inv;
  memcpy(rcopy,pinv+m*rowchanged,m*sizeof(T));
  dger_(m,m,-1.0,rcopy,1,temp,1,pinv,m);
}

int main(int argc, char** argv)
{
  int cmax=32;
  int M=512;
  typedef Matrix<double> mat_t;
  vector<mat_t*> dets(8);
  mat_t psi_big(M+cmax,M), psi_test(M,M), dpsi(M,M);
  for(int i=0; i<psi_big.size(); ++i)
    psi_big(i)=Random();
  for(int i=0; i<dets.size(); ++i)
    dets[i]=new mat_t(M,M);
  Timer myclock;
  double dt0=0.0,dt1=0.0,dt2=0.0;
  for(int iter=0; iter<100; ++iter)
  {
    std::copy(psi_big[0],psi_big[M],dets[0]->data());
    *dets[1] = *dets[0];
    int unoccupied=iter%cmax;
    double det_0=invert_matrix(*dets[0],true);
    psi_test=*dets[0];
    double r=BLAS::dot(M,(*dets[0])[unoccupied],psi_big[M]);
    myclock.restart();
    det_row_update(dets[0]->data(),psi_big[M],M,unoccupied,r);
    dt0 += myclock.elapsed();
    myclock.restart();
    det_ger(psi_test.data(),psi_big[M],M,unoccupied,r);
    dt1 += myclock.elapsed();
    for(int j=0; j<M; ++j)
      dets[1]->operator()(j,unoccupied)=psi_big(M,j);
    myclock.restart();
    double det_1=invert_matrix(*dets[1],true);
    dt2 += myclock.elapsed();
    //difference
    dpsi=psi_test-(*dets[0]);
    cout << setw(20) << BLAS::norm2(psi_test.size(),psi_test.data())
         << setw(20) << BLAS::norm2(dpsi.size(),dpsi.data()) << endl;
    dpsi=psi_test-(*dets[1]);
    cout << setw(20) << BLAS::norm2(psi_test.size(),psi_test.data())
         << setw(20) << BLAS::norm2(dpsi.size(),dpsi.data()) << endl;
  }
  cout << "Timing old = " << dt0 << " dger=" << dt1 << " direct=" << dt2
       << "\n dt0/dt1=" << dt0/dt1 << " dt2/dt1=" << dt2/dt1 << endl;
  delete_iter(dets.begin(),dets.end());
  return 0;
}
