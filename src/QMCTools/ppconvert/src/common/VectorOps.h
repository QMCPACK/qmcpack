//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

// extern "C"{
// #ifdef USE_MKL
//   #include <mkl_cblas.h>
// #else  
//   #include <cblas.h> 
// #endif
// }
#include "Blitz.h"
#include "../config.h"

typedef TinyVector<int,3> Int3;
typedef Array<std::complex<double>,1> zVec;
typedef Array<cVec3,1> zVecVec;
typedef Array<cMat3,1> zMatVec;

#ifdef USE_CBLAS
extern "C"{
  #ifdef USE_MKL_CBLAS
    #include <mkl_cblas.h>
  #else
    #include <cblas.h>
  #endif
}
#else

#define F77_DZNRM2 F77_FUNC(dznrm2,DZNRM2)
#define F77_ZDSCAL F77_FUNC(zdscal,ZDSCAL)
#define F77_ZDOTC  F77_FUNC(zdotc,ZDOTC)
#define F77_ZGEMV  F77_FUNC(zgemv,ZGEMV)
#define F77_ZAXPY  F77_FUNC(zaxpy,ZAXPY)

extern "C" double F77_DZNRM2(const int *N, const void *X, const int *INC);
extern "C" void   F77_ZDSCAL(const int *N, double *ALPHA, const void *X, 
			     const int *INC);
// extern "C" void F77_ZDOTC (std::complex<double> *z, const int *N, 
// 			   const void *X, const int *INCX, 
// 			   const void *Y, const int *INCY);
extern "C" std::complex<double> F77_ZDOTC (const int *N, 
 				      const void *X, const int *INCX, 
 				      const void *Y, const int *INCY);
extern "C" void   F77_ZGEMV (char *TRANS, const int *M, const int *N, 
			     std::complex<double> *alpha, const void *A, 
			     const int *LDA, const void *X, 
			     const int *INCX, std::complex<double> *beta, 
			     const void *Y, const int *INCY);

extern "C" void   F77_ZAXPY (const int *N, std::complex<double> *ALPHA,
			     void *X, int *INCX, void *Y, int *INCY);
#endif


#ifdef USE_CBLAS
inline void Normalize (zVec &c)
{
  double norm = cblas_dznrm2(c.size(), c.data(), 1);
  norm = 1.0/norm;
  cblas_zdscal (c.size(), norm, c.data(), 1);
}

inline double norm (const zVec &c)
{ return cblas_dznrm2(c.size(), c.data(), 1); }

inline std::complex<double> conjdot(zVec &cA, zVec &cB)
{
  std::complex<double> z;
  cblas_zdotc_sub(cA.size(), cA.data(), 1, cB.data(), 1, &z);
  return z;
}

inline void 
Orthogonalize (const Array<std::complex<double>,2> &A, zVec &x)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  std::complex<double> zero(0.0, 0.0);
  std::complex<double> one (1.0, 0.0);
  std::complex<double> minusone (-1.0, 0.0);
  Array<std::complex<double>,1> S(m);
  

  cblas_zgemv(CblasColMajor, CblasConjTrans, n, m, &one,
	      A.data(), n, x.data(), 1, &zero, S.data(), 1);
  cblas_zgemv(CblasRowMajor, CblasTrans, m, n, &minusone,
	      A.data(), n, S.data(), 1, &one, x.data(), 1);
  
}

#else

inline void Normalize (zVec &c)
{
  const int inc=1;
  int n = c.size();
  double norm = F77_DZNRM2(&n, c.data(), &inc);
  norm = 1.0/norm;
  F77_ZDSCAL(&n, &norm, c.data(), &inc);
}


inline double norm (const zVec &c)
{
  double norm;
  int n = c.size();
  const int inc=1;
  return F77_DZNRM2(&n, c.data(), &inc);
}


inline std::complex<double> conjdot(zVec &cA, zVec &cB)
{
  const int n = cA.size();
  const int incA = 1;
  const int incB = 1;
  return F77_ZDOTC (&n, cA.data(), &incA, cB.data(), &incB);
  
//   std::complex<double> z;
//   F77_ZDOTC (&z, &n, cA.data(), &incA, cB.data(), &incB);
//   return z;
}

inline void 
Orthogonalize (const Array<std::complex<double>,2> &A, zVec &x)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  std::complex<double> zero(0.0, 0.0);
  std::complex<double> one (1.0, 0.0);
  std::complex<double> minusone (-1.0, 0.0);
  Array<std::complex<double>,1> S(m);
  
  // Calculate overlaps
  // Calling with column major and ConjTrans is equivalent to
  // conjugate of untransposed row major
  char trans='C';
  const int inc=1;

  F77_ZGEMV(&trans, &n, &m, &one, A.data(), &n, x.data(), &inc,
	    &zero, S.data(), &inc);

//   cblas_zgemv(CblasColMajor, CblasConjTrans, n, m, &one,
// 	      A.data(), n, x.data(), 1, &zero, S.data(), 1);

//   for (int i=0; i<m; i++) {
//     fprintf (stderr, "S[%d] = %18.14f + %18.14fi\n",
// 	     real(S[i]), imag(S[i]));
    
  // Now, subtract off components * overlaps
  trans='T';
  F77_ZGEMV(&trans, &m, &n, &minusone, A.data(), &n,
	    S.data(), &inc, &one, x.data(), &inc);
//   cblas_zgemv(CblasRowMajor, CblasTrans, m, n, &minusone,
//  	      A.data(), n, S.data(), 1, &one, x.data(), 1);

}

#endif 

inline double realconjdot(zVec &cA, zVec &cB)
{
  return conjdot(cA,cB).real();
}

    


inline double mag (std::complex<double> x)
{
  return (x.real()*x.real() + x.imag()*x.imag());
}

inline void 
Orthogonalize2 (Array<std::complex<double>,2> &A, zVec &x, int lastBand)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;
  Array<std::complex<double>,1> S(m);

  for (int row=0; row<=lastBand; row++) {
    Ar.reference (A(row,Range::all()));
    S(row) = conjdot (Ar, x);
  }
  for (int row=0; row<=lastBand; row++) 
    x -= S(row) * A(row,Range::all());
//   for (int row=0; row<=lastBand; row++) {
//     Ar.reference (A(row,Range::all()));
//     S(row) = conjdot (Ar, x);
//     if (mag(S(row)) > 1.0e-14) {
//       std::cerr << "row = " << row << " lastband = " << lastBand << std::endl;
//       std::cerr << "Error in Orthogonalize2!, s = " << S(row) << std::endl;
//       double norm = realconjdot (Ar, Ar);
//       std::cerr << "norm = " << norm << std::endl;
//     }
//   }
}

inline void
OrthogExcluding(const Array<std::complex<double>,2> &A, zVec &x,
		int excluding)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;

#ifndef __INTEL_COMPILER
  std::complex<double> S[m];
  for (int row=0; row<m; row++) {
    Ar.reference (A(row,Range::all()));
    S[row] = conjdot(Ar, x);
  }
  for (int row=0; row<m; row++) 
    if (row != excluding)
      x -= S[row] * A(row,Range::all());
#else

  double Sre[m], Sim[m];

  for (int row=0; row<m; row++) {
    Ar.reference (A(row,Range::all()));
    std::complex<double> S = conjdot(Ar, x);
    Sre[row] = real(S);
    Sim[row] = imag(S);
  }
  for (int row=0; row<m; row++) 
    if (row != excluding)
      x -= std::complex<double>(Sre[row],Sim[row]) * A(row,Range::all());
#endif
}


inline void
OrthogLower(const Array<std::complex<double>,2> &A, zVec &x,
	    int currBand)
{
  int m = currBand;
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;
  std::complex<double> S;

  for (int row=0; row<m; row++) {
    Ar.reference (A(row,Range::all()));
    S = conjdot (Ar, x);
    x -= S * A(row,Range::all());
  }
}


inline void
GramSchmidt (Array<std::complex<double>,2> &A)
{
  zVec a, b;
  for (int i=0; i<A.rows(); i++) {
    a.reference (A(i,Range::all()));
    Normalize(a);
    for (int j=i+1; j<A.rows(); j++) {
      b.reference(A(j,Range::all()));
      b = b - (conjdot(a,b)*a);
      Normalize(b);
    }
  }
}

inline void
Orthogonalize (Array<std::complex<double>,2> &A)
{
  zVec x, y;
  Array<std::complex<double>,1> S(A.rows());
  for (int iter=0; iter < 40; iter++) {
    for (int i=0; i<A.rows();i++) {
      x.reference (A(i, Range::all()));
      for (int j=i+1; j<A.rows(); j++) {
	y.reference (A(j,Range::all()));
	S(j) = conjdot (y, x);
      }
      for (int j=i+1; j<A.rows(); j++) {
      y.reference (A(j,Range::all()));
      x -= S(j) * y;
      }
      Normalize (x);
    }
  }
}


inline void CheckOrthog (const Array<std::complex<double>,2> &A,
			 zVec &x)
{
  zVec Ai;
  double normInv = 1.0/norm(x);
  for (int i=0; i<A.rows(); i++) {
    Ai.reference (A(i,Range::all()));
    if (normInv*mag(conjdot(Ai, x)) > 1.0e-13) {
      std::cerr << "CheckOrthog failed for i=" << i << ".\n";
      exit(1);
    }
  }
}

inline void
zaxpy (std::complex<double> alpha, const Array<std::complex<double>,1> &x,
       const std::complex<double> y, Array<std::complex<double>,1> &axpy)
{

}


// inline Array<std::complex<double>,3>&
// operator*= (Array<std::complex<double>,3> &A, const Array<std::complex<double>,3> &B)
// {
  

// }


#endif
