//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef LINEAR_FIT_H
#define LINEAR_FIT_H

#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus
{

template<typename T>
inline void LinearFit (std::vector<T> &y, Matrix<T> &A, std::vector<T> &coefs)
{
  int N = A.size(0);
  int M = A.size(1);
  if (y.size() != N)
  {
    app_error() << "Differernt number of rows in basis functions that in data "
                << "points in LinearFit.  Exitting.\n";
    abort();
  }
  // Construct alpha matrix
  Matrix<T> alpha(M,M);
  alpha = 0.0;
  for (int j=0; j<M; j++)
    for (int k=0; k<M; k++)
      for (int i=0; i<N; i++)
        alpha(k,j) += A(i,j) * A(i,k);
  // Next, construct beta vector
  std::vector<T> beta(M,0.0);
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta[k] += y[i]*A(i,k);
  invert_matrix(alpha, false);
  coefs.resize (M, 0.0);
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      coefs[i] += alpha(i,j)*beta[j];
}
}

#endif
