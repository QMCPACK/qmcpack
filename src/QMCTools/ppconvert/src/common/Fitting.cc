//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "Fitting.h"
#include "MatrixOps.h"

inline Array<double,1> operator*(Array<double,2> &A, Array<double,1> &x)
{
  assert (A.cols() == x.size());
  Array<double,1> Ax(A.rows());
  Ax = 0.0;
  for (int i=0; i<A.rows();i++)
    for (int j=0; j<A.cols(); j++)
      Ax(i) += A(i,j) * x(j);
  return Ax;
}

inline Array<double,1> diag (Array<double,2> &A)
{
  assert (A.rows() == A.cols());
  Array<double,1> d(A.rows());
  for (int i=0; i<A.rows(); i++)
    d(i) = A(i,i);
  return (d);
}


void LinFitLU (Array<double,1> &y, Array<double,1> &sigma,    // inputs
	       Array<double,2> &F,                     // input
	       Array<double,1> &a, Array<double,1> &errors) // outputs
{
  int N = F.rows();
  int M = F.cols();

  if (y.rows() != F.rows()) {
    cerr << "Differernt number of rows in basis functions that in data "
	 << "points in LinFit.  Exitting.\n";
    abort();
  }

  assert (y.rows() == sigma.rows());
  a.resize (F.cols());
  errors.resize(F.cols());

  // First, construct A matrix
  Array<double,2> A(F.rows(), F.cols());
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)
      A(i,j) = F(i,j) / sigma(i);
  
  // Next, construct b vector
  Array<double,1> b(y.rows());
  for (int i=0; i<N; i++)
    b(i) = y(i)/sigma(i);
  
  // Next, construct alpha matrix
  Array<double,2> alpha(M,M);
  alpha = 0.0;
  for (int j=0; j<M; j++)
    for (int k=0; k<M; k++)
      for (int i=0; i<N; i++)
	alpha(k,j) += A(i,j) * A(i,k);
  
  // Next, construct beta vector
  Array<double,1> beta(M);
  beta = 0.0;
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta(k) += b(i)*A(i,k);

  // Now, invert alpha
  Array<double,2> C(M,M);
  C = Inverse(alpha);
  a = C * beta;
  for (int i=0; i<M; i++)
    errors(i) = sqrt(C(i,i));
  //  errors = sqrt(diag(C));  
}




void LinFitSVD (Array<double,1> &y, Array<double,1> &sigma,  // inputs
		Array<double,2> &F,                          // input
		Array<double,1> &a, Array<double,1> &errors, // outputs
		double tolerance)
{
  int N = F.rows();
  int M = F.cols();

  if (y.rows() != F.rows()) {
    cerr << "Differernt number of rows in basis functions that in data "
	 << "points in LinFit.  Exitting.\n";
    abort();
  }

  assert (y.rows() == sigma.rows());
  a.resize (F.cols());
  errors.resize(F.cols());

  // First, construct A matrix
  Array<double,2> A(F.rows(), F.cols());
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)
      A(i,j) = F(i,j) / sigma(i);
  
  // Next, construct b vector
  Array<double,1> b(y.rows());
  for (int i=0; i<N; i++)
    b(i) = y(i)/sigma(i);
  
  // Now, compute SVD
  Array<double,2> U,V;
  Array<double,1> S;
  SVdecomp (A, U, S, V);

  //  cerr << "S = " << S << endl;
  // Zero out near-singular values
  double Smax=S(0);
  for (int i=1; i<S.size(); i++)
    Smax = max (S(i),Smax);
  Array<double,1> Sinv(S.size());
  for (int i=0; i<S.size(); i++)
    Sinv(i) = (S(i) < (tolerance*Smax)) ? 0.0 : (1.0/S(i));
  
  a = 0.0;
  for (int k=0; k<a.size(); k++) 
    for (int i=0; i<U.cols(); i++) {
      double coef = 0.0;
      for (int j=0; j < U.rows(); j++)
	coef += U(j,i) * b(j);
      coef *= Sinv(i);
      a(k) += coef * V(k,i);
    }
  
  errors = 0.0;
  for (int j=0; j<errors.size(); j++) {
    for (int i=0; i<V.cols(); i++)
      errors(j) += V(j,i)*V(j,i)*(Sinv(i)*Sinv(i));
    errors(j) = sqrt(errors(j));
  }
}

// This is a version of the above that allows certain basis functions
// to be pinned and not adjusted
double
LinFitSVD (Array<double,1> &y, Array<double,1> &sigma,  // inputs
	   Array<double,2> &F, Array<bool,  1> &adjust,  // inputs
	   Array<double,1> &t, Array<double,1> &errors, // outputs
	   double tolerance)
{
  int N = F.rows(); // Number of data points
  int M = F.cols(); // Number of basis functions

  if (y.rows() != F.rows()) {
    cerr << "Different number of rows in basis functions that in data "
	 << "points in LinFit.  Exitting.\n";
    abort();
  }

  assert (y.rows() == sigma.rows());
  assert (t.size() == F.cols());
  errors.resize(F.cols());

  // First, construct A matrix
  Array<double,2> A(F.rows(), F.cols());
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)
      A(i,j) = F(i,j) / sigma(i);
  
  // Next, construct b vector
  Array<double,1> b(y.rows());
  for (int i=0; i<N; i++)
    b(i) = y(i)/sigma(i);
  
  // Now reduce for constraints
  int P = M;
  for (int i=0; i<adjust.size(); i++) 
    if (!adjust(i))
      P--;

  // The c is for "constrained"
  Array<double,2> Ac(N,P);
  Array<double,1> bc(N), tc(P);

  // Build constrained Ac and bc
  int j=0;
  for (int col=0; col<M; col++) {
    if (adjust(col)) {
      // Copy column a A to Ac
      for (int row=0; row<N; row++) 
	Ac(row,j) = A(row,col);
      j++;
    }
    else {
      // Otherwise, subtract t(col)*A(:,col) from bc
      for (int row=0; row<N; row++)
      	b(row) -= A(row,col)*t(col);
    }
  }
  for (int row=0; row<N; row++)
    bc(row) = b(row);

  // Now, compute SVD
  Array<double,2> U (P,P), V(P,P);
  Array<double,1> S(P), Sinv(P);
  SVdecomp (Ac, U, S, V);

  // Zero out near-singular values
  double Smax=S(0);
  for (int i=1; i<S.size(); i++)
    Smax = max (S(i),Smax);
  for (int i=0; i<S.size(); i++)
    Sinv(i) = (S(i) < (tolerance*Smax)) ? 0.0 : (1.0/S(i));
  
  tc = 0.0;
  // Compute t_n, removing singular values
  for (int k=0; k<tc.size(); k++) 
    for (int i=0; i<U.cols(); i++) {
      double coef = 0.0;
      for (int j=0; j < U.rows(); j++)
	coef += U(j,i) * bc(j);
      coef *= Sinv(i);
      tc(k) += coef * V(k,i);
    }
  
  Array<double,1> cerrors(P);
  
  cerrors = 0.0;
  for (int j=0; j<cerrors.size(); j++) {
    for (int i=0; i<V.cols(); i++)
      cerrors(j) += V(j,i)*V(j,i)*(Sinv(i)*Sinv(i));
    cerrors(j) = sqrt(errors(j));
  }

  // Now copy tc values into t
  j=0;
  for (int i=0; i<M; i++)
    if (adjust(i)) {
      t(i)      = tc(j);
      errors(i) = cerrors(j);
      j++;
    }
}
