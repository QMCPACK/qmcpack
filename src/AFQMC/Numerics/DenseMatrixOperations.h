#ifndef AFQMC_DENSEMATRIXOPERATORS_H
#define AFQMC_DENSEMATRIXOPERATORS_H

#include<iterator>
#include<complex>

//#include <AFQMC/Matrix/DenseMatrix.h>
#include "AFQMC/config.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/Blasf.h"
#include "Numerics/DeterminantOperators.h"

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

template<typename T>
inline void product(const int M, const int N, const int K, const T one, const T* A, const int LDA, const T* B, const int LDB, const T zero, T* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'N';

  // C = A*B -> fortran -> C' = B'*A', 
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);  
} 

template<typename T>
inline void product(const int M, const int N, const int K, const T one, const T* A, const int LDA, const std::complex<T>* B, const int LDB, const T zero, std::complex<T>* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'N';
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  // C = A*B -> fortran -> C' = B'*A', 
  BLAS::gemm(transa,transb, 2*N, M, K,
            one, B_, 2*LDB, A, LDA,
            zero, C_, 2*LDC);
}

template<typename T>
inline void product_AtB(const int M, const int N, const int K, const T one, const T* A, const int LDA, const std::complex<T>* B, const int LDB, const T zero, std::complex<T>* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'T';
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  // C = A'*B -> fortran -> C' = B'*A, 
  BLAS::gemm(transa,transb, 2*N, M, K,
            one, B_, 2*LDB, A, LDA,
            zero, C_, 2*LDC);
} 

template<typename T>
inline void product_AtB(const int M, const int N, const int K, const T one, const T* A, const int LDA, const T* B, const int LDB, const T zero, T* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'T';

  // C = A'*B -> fortran -> C' = B'*A, 
  BLAS::gemm(transa,transb, N, M, K,
            one, B, LDB, A, LDA,
            zero, C, LDC);
}

template<typename T>
inline void product_AtBt(const int M, const int N, const int K, const T one, const T* A, const int LDA, const T* B, const int LDB, const T zero, T* C, const int LDC )
{
  const char transa = 'T';
  const char transb = 'T';

  // C = A'*B' -> fortran -> C' = B*A, 
  BLAS::gemm(transa,transb, N, M, K,
            one, B, LDB, A, LDA,
            zero, C, LDC);
}

template<typename T>
inline void product_ABt(const int M, const int N, const int K, const T one, const T* A, const int LDA, const T* B, const int LDB, const T zero, T* C, const int LDC )
{
  const char transa = 'T';
  const char transb = 'N';

  // C = A*B' -> fortran -> C' = B*A', 
     BLAS::gemm(transa,transb, N, M, K,
            one, B, LDB, A, LDA,
            zero, C, LDC);
  
}

template<typename T>
inline void product_AhB(const int M, const int N, const int K, const T one, const T* A, const int LDA, const T* B, const int LDB, const T zero, T* C, const int LDC )
{
  const char transa = 'N';
  const char transb = 'C';

  // C = A'*B -> fortran -> C' = B'*A, 
     BLAS::gemm(transa,transb, N, M, K,
            one, B, LDB, A, LDA,
            zero, C, LDC);
  
}

template<typename T>
inline void product_ABh(const int M, const int N, const int K, const T one, const T* A, const int LDA, const T* B, const int LDB, const T zero, T* C, const int LDC )
{
  const char transa = 'C';
  const char transb = 'N';

  // C = A*B^H -> fortran -> C' = conjg(B)*A', 
  BLAS::gemm(transa,transb, N, M, K,
          one, B, LDB, A, LDA,
          zero, C, LDC);

}

template<typename T>
inline void product_Ax(const int N, const int K, const T one, const T* A, const int lda, const T* x, const T zero, T* restrict yptr)
{
  const char transa = 'T';
  BLAS::gemv(transa, K, N, one, A, lda, x, 1, zero, yptr, 1);
}

template<typename T>
inline void product_Ax(const int M, const int K, const T one, const T* A, const int lda, const std::complex<T>* x, const T zero, std::complex<T>* restrict yptr)
{
// add specialized routine if it exists 
#if defined(HAVE_MKL)
  const std::complex<T> one_ = one;
  const std::complex<T> zero_ = zero;
  const char transa = 'T';

  BLAS::gemv(transa, K, M, one_, A, lda, x, 1, zero_, yptr, 1);

#else 
  const T* x_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(x));
  T* yptr_ = reinterpret_cast<T*>(yptr);
  product(M,2,K,one,A,lda,x_,2,zero,yptr_,2);
#endif
}

template<typename T>
inline void product_Atx(const int N, const int K, const T one, const T* A, const int lda, const T* x, const T zero, T* restrict yptr)
{
  const char transa = 'N';
  BLAS::gemv(transa, K, N, one, A, lda, x, 1, zero, yptr, 1);
}

template<typename T>
inline void product_Atx(const int M, const int K, const T one, const T* A, const int lda, const std::complex<T>* x, const T zero, std::complex<T>* restrict yptr)
{
// add specialized routine if it exists 
#if defined(HAVE_MKL)
  const std::complex<T> one_ = one;
  const std::complex<T> zero_ = zero;
  const char transa = 'N';

  BLAS::gemv(transa, K, M, one_, A, lda, x, 1, zero_, yptr, 1);

#else
  const T* x_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(x));
  T* yptr_ = reinterpret_cast<T*>(yptr);
  product_AtB(K,2,M,one,A,lda,x_,2,zero,yptr_,2);
#endif
}
/*
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
*/


struct SmallDet
{

  template<typename T>
  inline static
  T D2x2(T a11, T a12,
               T a21, T a22)
  {
    return a11*a22-a21*a12;
  }

  template<typename T>
  inline static
  T D3x3(T a11, T a12, T a13,
        T a21, T a22, T a23,
        T a31, T a32, T a33)
  {
    return (a11*(a22*a33-a32*a23)-a21*(a12*a33-a32*a13)+a31*(a12*a23-a22*a13));
  }

  template<typename T>
  inline static
  T I3x3(T a11, T a12, T a13,
        T a21, T a22, T a23,
        T a31, T a32, T a33, T* M)
  {
    T det = (a11*(a22*a33-a32*a23)-a21*(a12*a33-a32*a13)+a31*(a12*a23-a22*a13));
    M[0] = (a22*a33-a32*a23)/det; 
    M[1] = (a13*a32-a12*a33)/det;
    M[2] = (a12*a23-a13*a22)/det;
    M[3] = (a23*a31-a21*a33)/det;
    M[4] = (a11*a33-a13*a31)/det;
    M[5] = (a13*a21-a11*a23)/det;
    M[6] = (a21*a32-a22*a31)/det;
    M[7] = (a12*a31-a11*a32)/det;
    M[8] = (a11*a22-a12*a21)/det;
    return det; 
  }

  template<typename T>
  inline static
  T D4x4(T a11, T a12, T a13, T a14,
        T a21, T a22, T a23, T a24,
        T a31, T a32, T a33, T a34,
        T a41, T a42, T a43, T a44)
  {
    return (a11*(a22*(a33*a44-a43*a34)-a32*(a23*a44-a43*a24)+a42*(a23*a34-a33*a24))-a21*(a12*(a33*a44-a43*a34)-a32*(a13*a44-a43*a14)+a42*(a13*a34-a33*a14))+a31*(a12*(a23*a44-a43*a24)-a22*(a13*a44-a43*a14)+a42*(a13*a24-a23*a14))-a41*(a12*(a23*a34-a33*a24)-a22*(a13*a34-a33*a14)+a32*(a13*a24-a23*a14)));
  }

  template<typename T>
  inline static
  T D5x5(T a11, T a12, T a13, T a14, T a15,
        T a21, T a22, T a23, T a24, T a25,
        T a31, T a32, T a33, T a34, T a35,
        T a41, T a42, T a43, T a44, T a45,
        T a51, T a52, T a53, T a54, T a55)
  {
    return (a11*(a22*(a33*(a44*a55-a54*a45)-a43*(a34*a55-a54*a35)+a53*(a34*a45-a44*a35))-a32*(a23*(a44*a55-a54*a45)-a43*(a24*a55-a54*a25)+a53*(a24*a45-a44*a25))+a42*(a23*(a34*a55-a54*a35)-a33*(a24*a55-a54*a25)+a53*(a24*a35-a34*a25))-a52*(a23*(a34*a45-a44*a35)-a33*(a24*a45-a44*a25)+a43*(a24*a35-a34*a25)))-a21*(a12*(a33*(a44*a55-a54*a45)-a43*(a34*a55-a54*a35)+a53*(a34*a45-a44*a35))-a32*(a13*(a44*a55-a54*a45)-a43*(a14*a55-a54*a15)+a53*(a14*a45-a44*a15))+a42*(a13*(a34*a55-a54*a35)-a33*(a14*a55-a54*a15)+a53*(a14*a35-a34*a15))-a52*(a13*(a34*a45-a44*a35)-a33*(a14*a45-a44*a15)+a43*(a14*a35-a34*a15)))+a31*(a12*(a23*(a44*a55-a54*a45)-a43*(a24*a55-a54*a25)+a53*(a24*a45-a44*a25))-a22*(a13*(a44*a55-a54*a45)-a43*(a14*a55-a54*a15)+a53*(a14*a45-a44*a15))+a42*(a13*(a24*a55-a54*a25)-a23*(a14*a55-a54*a15)+a53*(a14*a25-a24*a15))-a52*(a13*(a24*a45-a44*a25)-a23*(a14*a45-a44*a15)+a43*(a14*a25-a24*a15)))-a41*(a12*(a23*(a34*a55-a54*a35)-a33*(a24*a55-a54*a25)+a53*(a24*a35-a34*a25))-a22*(a13*(a34*a55-a54*a35)-a33*(a14*a55-a54*a15)+a53*(a14*a35-a34*a15))+a32*(a13*(a24*a55-a54*a25)-a23*(a14*a55-a54*a15)+a53*(a14*a25-a24*a15))-a52*(a13*(a24*a35-a34*a25)-a23*(a14*a35-a34*a15)+a33*(a14*a25-a24*a15)))+a51*(a12*(a23*(a34*a45-a44*a35)-a33*(a24*a45-a44*a25)+a43*(a24*a35-a34*a25))-a22*(a13*(a34*a45-a44*a35)-a33*(a14*a45-a44*a15)+a43*(a14*a35-a34*a15))+a32*(a13*(a24*a45-a44*a25)-a23*(a14*a45-a44*a15)+a43*(a14*a25-a24*a15))-a42*(a13*(a24*a35-a34*a25)-a23*(a14*a35-a34*a15)+a33*(a14*a25-a24*a15))));
  }
};

template<typename T>
inline T DeterminantSmall(T* M, int n, int* restrict pivot)
{
  switch(n)
  {
  case 0:
    return 1.0;
  case 1:
    return M[0]; 
    break;
  case 2:
  {
    return M[0]*M[3] - M[1]*M[2]; 
    break;
  }
  case 3:
  {
    return SmallDet::D3x3(M[0], M[1], M[2],
                        M[3], M[4], M[5],  
                        M[6], M[7], M[8]); 
    break;
  }
  case 4:
  {
    return SmallDet::D4x4(M[0], M[1], M[2], M[3],
                        M[4], M[5], M[6], M[7],
                        M[8], M[9], M[10], M[11], 
                        M[12], M[13], M[14], M[15]); 
    break;
    break;
  }
  case 5:
  {
    return SmallDet::D5x5(M[0], M[1], M[2], M[3], M[4],
                        M[5], M[6], M[7], M[8], M[9], 
                        M[10], M[11], M[12], M[13], M[14], 
                        M[15], M[16], M[17], M[18], M[19], 
                        M[20], M[21], M[22], M[23], M[24]); 
    break;
    break;
  }
  default:
  {
    return Determinant(M,n,n,pivot); 
  }
  }
  return 0.0;
}

template<typename T>
inline T InvertSmall(T* M, int n, T* restrict work, int* restrict pivot)
{
  switch(n)
  {
  case 0:
  {
    return T(1);
    break; 
  }
  case 1:
  {
    T det = M[0];
    M[0] = T(1)/det;
    return det;
    break;
  }
  case 2:
  {
    T one = T(1);
    T det = M[0]*M[3] - M[1]*M[2]; 
    std::swap(M[0],M[3]);  
    M[0]*=one/det;  
    M[3]*=one/det;  
    M[1]*=-one/det;  
    M[2]*=-one/det;  
    return det;
    break;
  }
  case 3:
  {
    return SmallDet::I3x3(M[0], M[1], M[2],
                        M[3], M[4], M[5],
                        M[6], M[7], M[8],M);
    break;
  }
  }
  return Invert(M,n,n,work,pivot);
}

void GeneralizedGramSchmidt(std::complex<double>* A, int LDA, int nR, int nC);

} // namespace DenseMatrixOperators

} // namespace qmcplusplus

#endif

