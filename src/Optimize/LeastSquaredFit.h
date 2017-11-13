//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_LEASTSQUARED_FITTING_H
#define QMCPLUSPLUS_LEASTSQUARED_FITTING_H

#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"

template<typename VT, typename MT>
void LeastSquaredFitLU(VT &y, VT &sigma, MT &F, VT &a, VT &errors)
{
  int N = F.rows();
  int M = F.cols();
  a.resize (M);
  errors.resize(M);
  MT A(N,M);
  for (int i=0; i<N; ++i)
    for (int j=0; j<M; ++j)
      A(i,j) = F(i,j) / sigma[i];
  VT b(N);
  for (int i=0; i<N; ++i)
    b[i] = y[i]/sigma[i];
  MT alpha(M,M);
  alpha = 0.0;
  for (int j=0; j<M; ++j)
    for (int k=0; k<M; ++k)
      for (int i=0; i<N; ++i)
        alpha(k,j) += A(i,j) * A(i,k);
  VT beta(M);
  beta = 0.0;
  for (int k=0; k<M; ++k)
    for (int i=0; i<N; ++i)
      beta[k] += b[i]*A(i,k);
  qmcplusplus::invert_matrix(alpha,false);
  BLAS::gemv(M,M,alpha.data(),beta.data(),a.data());
  for (int i=0; i<M; ++i)
    errors[i] = std::sqrt(alpha(i,i));
}
#endif
