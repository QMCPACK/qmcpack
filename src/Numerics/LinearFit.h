#ifndef LINEAR_FIT_H
#define LINEAR_FIT_H

#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus
{

template<typename T>
inline void LinearFit (vector<T> &y, Matrix<T> &A, vector<T> &coefs)
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
  vector<T> beta(M,0.0);
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
