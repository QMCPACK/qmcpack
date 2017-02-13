#ifndef AFQMC_DENSEMATRIXOPERATORS_H
#define AFQMC_DENSEMATRIXOPERATORS_H

#include<iterator>
#include<complex>

//#include <AFQMC/Matrix/DenseMatrix.h>
#include "AFQMC/config.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/Blasf.h"

namespace qmcplusplus
{

namespace DenseMatrixOperators
{

inline bool isHermitian(int N, std::complex<double>* A, int LDA)
{
  for(int i=0; i<N; i++) 
   for(int j=i+1; j<N; j++) 
    if( A[i*LDA+j] != myconj(A[j*LDA+i]) )
     return false; 
  return true;
}

inline bool isHermitian(int N, double* A, int LDA)
{
  for(int i=0; i<N; i++)
   for(int j=i+1; j<N; j++)
    if( A[i*LDA+j] != A[j*LDA+i] )
     return false;
  return true;
}

inline bool isHermitian(Matrix<std::complex<double> >& A)
{
  if(A.rows() != A.cols()) return false;
  for(int i=0; i<A.rows(); i++)
   for(int j=i+1; j<A.cols(); j++)
    if( A(i,j) != myconj(A(j,i)) )
     return false;
  return true;
}

inline bool isSymmetric(int N, std::complex<double>* A, int LDA)
{
  for(int i=0; i<N; i++) 
   for(int j=i+1; j<N; j++) 
    if( A[i*LDA+j] != A[j*LDA+i] )
     return false; 
  return true;
}

inline bool isSymmetric(Matrix<std::complex<double> >& A)
{
  if(A.rows() != A.cols()) return false;
  for(int i=0; i<A.rows(); i++)
   for(int j=i+1; j<A.cols(); j++)
    if( A(i,j) != A(j,i) )
     return false;
  return true;
}

template<typename T>
inline void transpose(int N, T* A, int LDA ) {
  for (int i=0; i<N; i++)
    for (int j=0; j<i; j++)
      std::swap(A[i*LDA+j],A[j*LDA+i]);
}

bool exponentiateHermitianMatrix(int N, std::complex<double>* A, int LDA, std::complex<double>* expA, int LDEXPA); 

bool symEigenSysAll(int N, std::complex<double>* A, int LDA, double* eigVal, std::complex<double>* eigVec, int LDV); 
bool symEigenSysAll(int N, double* A, int LDA, double* eigVal, double* eigVec, int LDV); 

bool symEigenSysSelect(int N, double* A, int LDA, int neig, double* eigVal, bool getEigV, double* eigVec, int LDV); 
bool symEigenSysSelect(int N, std::complex<double>* A, int LDA, int neig, double* eigVal, bool getEigV, std::complex<double>* eigVec, int LDV); 

bool genHermitianEigenSysSelect(int N, std::complex<double>* A, int LDA, std::complex<double>* B, int LDB, int neig, double* eigVal, bool getEigV, std::complex<double>* eigVec, int LDV, int* ifail);

inline void product(const int M, const int N, const int K, const std::complex<double>* A, const int LDA, const std::complex<double>* B, const int LDB, std::complex<double>* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'N';
  const std::complex<double> one=1.0;
  const std::complex<double> zero=0.0;  

  // C = A*B -> fortran -> C' = B'*A', 
  BLAS::gemm(transa,transb, N, M, K,
             one, B, LDB, A, LDA,
             zero, C, LDC);  
}

inline void product(const int M, const int N, const int K, const std::complex<double> one, const std::complex<double>* A, const int LDA, const std::complex<double>* B, const int LDB, const std::complex<double> zero, std::complex<double>* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'N';

  // C = A*B -> fortran -> C' = B'*A', 
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);  
}

inline void product(const int M, const int N, const int K, const double* A, const int LDA, const double* B, const int LDB, double* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'N';
  const double one=1.0;
  const double zero=0.0;  

  // C = A*B -> fortran -> C' = B'*A',
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);  
}

inline void product(const int M, const int N, const int K, const double one, const double* A, const int LDA, const double* B, const int LDB, const double zero, double* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'N';

  // C = A*B -> fortran -> C' = B'*A', 
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);  
} 

inline void product_AhB(const int M, const int N, const int K, const double one, const double* A, const int LDA, const double* B, const int LDB, const double zero, double* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'T';

  // C = A'*B -> fortran -> C' = B'*A, 
     BLAS::gemm(transa,transb, N, M, K,
            one, B, LDB, A, LDA,
            zero, C, LDC);
  
}

inline void product_AhB(const int M, const int N, const int K, const std::complex<double> one, const std::complex<double>* A, const int LDA, const std::complex<double>* B, const int LDB, const std::complex<double> zero, std::complex<double>* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'C';

  // C = A^H*B -> fortran -> C' = B'*conjg(A), 
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);

}

inline void product_ABh(const int M, const int N, const int K, const std::complex<double> one, const std::complex<double>* A, const int LDA, const std::complex<double>* B, const int LDB, const std::complex<double> zero, std::complex<double>* C, const int LDC )
{
  const char transa = 'C';
  const char transb = 'N';

  // C = A*B^H -> fortran -> C' = conjg(B)*A', 
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);

}

inline std::complex<double> 
Determinant(std::complex<double>* restrict x, int n, int* restrict pivot)
{
  std::complex<double> detvalue(1.0);
  int status;
  zgetrf(n,n,x,n,pivot,status);
  for(int i=0,ip=1; i<n; i++, ip++)
  {
    if(pivot[i]==ip)
      detvalue *= x[i*n+i];
    else
      detvalue *= -x[i*n+i];
  }
  return detvalue;
}

/*
void product(const int M, const int N, const int K, const std::complex<double>,  const std::complex<double>* A, const int LDA, const std::complex<double>* B, const int LDB, const std::complex<double>, std::complex<double>* C, const int LDC );

void product(const int M, const int N, const int K, const std::complex<double>* A, const int LDA, const std::complex<double>* B, const int LDB, std::complex<double>* C, const int LDC );

void product(const int M, const int N, const int K, const double* A, const int LDA, const double* B, const int LDB, double* C, const int LDC );

void product(const int M, const int N, const int K, const double, const double* A, const int LDA, const double* B, const int LDB, const double, double* C, const int LDC );
*/

void GeneralizedGramSchmidt(std::complex<double>* A, int LDA, int nR, int nC);

} // namespace DenseMatrixOperators

} // namespace qmcplusplus

#endif

