#ifndef AFQMC_SPARSEMATRIXOPERATORS_H
#define AFQMC_SPARSEMATRIXOPERATORS_H

#include<iterator>
#include<tuple>

#include "AFQMC/config.h"

#include "AFQMC/Numerics/sparse.h"

namespace qmcplusplus
{

namespace SparseMatrixOperators
{

// Performs a product between a sparse matrix stored in format s2D and a dense 
// matrix stored in c format
//   N: number of rows in B/C 
//   M: number of columns in B/C
//   LDB: leading dimension of B
//
//   For eack term in A, aik 
//   C(i,:) += aik * B(k.:) 
template<typename T>
inline void product_SD(const IndexType K,
             const s2D<T>* A, const int nterms,
             ComplexType* B, const IndexType LDB,
             ComplexType* C, IndexType LDC )
{
  T aik=0;
  IndexType ii=0,kk=0;
  ComplexType* lit;
  ComplexType* rit;
  for(int cnt1=0; cnt1<nterms; cnt1++) {
    std::tie(ii,kk,aik) = *(A++);
    lit=C+ii*LDC;
    rit=B+kk*LDB;
    for(int cnt2=0; cnt2<K; cnt2++)
      *(lit++) += aik*(*(rit++));
  }

}

template<class T, typename T1>
inline void product_SpMatV(int N, 
             const T& A,
             const T1* B,
             T1* C )
{
  const char trans = 'N';
  SPBLAS::csrmv( trans, N, N, T1(1.0), "GxxCxx", A.values() , A.column_data(),  A.row_index(),  A.row_index()+1, B, T1(0.0), C );
}

template<typename T> 
inline void product_SpMatV(const int M, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const T* B,
             const T beta,
             T* C )
{
  const char trans = 'N';
  SPBLAS::csrmv( trans, M, K, alpha, "GxxCxx", val , col,  row,  row+1, B, beta, C );
}

template<typename T>
inline void product_SpMatV(const int M, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* rowb,
             const int* rowe,
             const T* B,
             const T beta,
             T* C )
{
  const char trans = 'N';
  SPBLAS::csrmv( trans, M, K, alpha, "GxxCxx", val , col,  rowb,  rowe, B, beta, C );
}

template<class T, typename T1>
inline void product_SpMatTV(const int M, const int K,
             const T1 alpha,
             const T& A,
             const T1* B,
             const T1 beta,
             T1* C )
{
  const char trans = 'T';
  SPBLAS::csrmv( trans, M, K, alpha, "GxxCxx", A.values() , A.column_data(),  A.row_index(),  A.row_index()+1, B, beta, C );
}

template<typename T>
inline void product_SpMatTV(const int M, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const T* B,
             const T beta,
             T* C )
{
  const char trans = 'T';
  SPBLAS::csrmv( trans, M, K, alpha, "GxxCxx", val , col,  row,  row+1, B, beta, C );
}

template<class T, typename T1>
inline void product_SpMatM(const int M, const int N, const int K,
             const T1 alpha,
             const T& A,
             const T1* B, const int ldb,
             const T1 beta,
             T1* C, int ldc )
{
  char trans = 'N';
  SPBLAS::csrmm( trans, M, N, K, alpha, "GxxCxx", A.values() , A.column_data(),  A.row_index(),  A.row_index()+1, B, ldb, beta, C, ldc);
}

template<typename T>
inline void product_SpMatM(const int M, const int N, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const T* B, const int ldb,
             const T beta,
             T* C, int ldc )
{
  char trans = 'N';
  SPBLAS::csrmm( trans, M, N, K, alpha, "GxxCxx", val, col,  row,  row+1, B, ldb, beta, C, ldc);
}

template<typename T>
inline void product_SpMatM(const int M, const int N, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* rowb,
             const int* rowe,
             const T* B, const int ldb,
             const T beta,
             T* C, int ldc )
{
  char trans = 'N';
  SPBLAS::csrmm( trans, M, N, K, alpha, "GxxCxx", val, col,  rowb,  rowe, B, ldb, beta, C, ldc);
}

template<class T, typename T1>
inline void product_SpMatTM(const int M, const int N, const int K,
             const T1 alpha,
             const T& A,
             const T1* B, const int ldb,
             const T1 beta,
             T1* C, int ldc )
{
  char trans = 'T';
  SPBLAS::csrmm( trans, M, N, K, alpha, "GxxCxx", A.values() , A.column_data(),  A.row_index(),  A.row_index()+1, B, ldb, beta, C, ldc);
}

template<typename T>
inline void product_SpMatTM(const int M, const int N, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const T* B, const int ldb,
             const T beta,
             T* C, int ldc )
{
  char trans = 'T';
  SPBLAS::csrmm( trans, M, N, K, alpha, "GxxCxx", val, col,  row,  row+1, B, ldb, beta, C, ldc);
}

template<typename T>
inline void product_SpMatM(const int M, const int N, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const std::complex<T>* B, const int ldb,
             const T beta,
             std::complex<T>* C, int ldc )
{
  char trans = 'N';
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  const int N_ = 2*N;
  const int ldb_ = 2*ldb;
  const int ldc_ = 2*ldc;
  SPBLAS::csrmm( trans, M, N_, K, alpha, "GxxCxx", val, col, row,  row+1, B_, ldb_, beta, C_, ldc_);
}

template<typename T>
inline void product_SpMatM(const int M, const int N, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* rowb,
             const int* rowe,
             const std::complex<T>* B, const int ldb,
             const T beta,
             std::complex<T>* C, int ldc )
{
  char trans = 'N';
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  const int N_ = 2*N;
  const int ldb_ = 2*ldb;
  const int ldc_ = 2*ldc;
  SPBLAS::csrmm( trans, M, N_, K, alpha, "GxxCxx", val, col, rowb,  rowe, B_, ldb_, beta, C_, ldc_);
}

template<typename T>
inline void product_SpMatTV(const int M, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const std::complex<T>* B,
             const T beta,
             std::complex<T>* C )
{
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  product_SpMatTM(M,2,K,alpha,val,col,row,B_,2,beta,C_,2);
}

template<typename T>
inline void product_SpMatV(const int M, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const std::complex<T>* B,
             const T beta,
             std::complex<T>* C )
{
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  product_SpMatM(M,2,K,alpha,val,col,row,B_,2,beta,C_,2);
}

template<typename T>
inline void product_SpMatV(const int M, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* rowb,
             const int* rowe,
             const std::complex<T>* B,
             const T beta,
             std::complex<T>* C )
{
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  product_SpMatM(M,2,K,alpha,val,col,rowb,rowe,B_,2,beta,C_,2);
}

template<typename T>
inline void product_SpMatTM(const int M, const int N, const int K,
             const T alpha,
             const T* val,
             const int* col,
             const int* row,
             const std::complex<T>* B, const int ldb,
             const T beta,
             std::complex<T>* C, int ldc )
{
  char trans = 'T';
  const T* B_ = reinterpret_cast<T*>(const_cast<std::complex<T>*>(B));
  T* C_ = reinterpret_cast<T*>(C);
  const int N_ = 2*N;
  const int ldb_ = 2*ldb;
  const int ldc_ = 2*ldc;
  SPBLAS::csrmm( trans, M, N_, K, alpha, "GxxCxx", val, col, row,  row+1, B_, ldb_, beta, C_, ldc_);
}

// Dot product between 2 sparse vectors
template<class T>
inline T product_SpVSpV(const int n1,const  int* indx1, const T* A1, const int n2, const int* indx2, const T* A2) {
  T res=T(0);
  int i=0, j=0;
  while( i<n1 && j<n2 ) {
    if( *(indx1+i) < *(indx2+j)   )
      ++i;
    else if( *(indx2+j) < *(indx1+i) )
      ++j;
    else {
      res += *(A1+i) * (*(A2+j));
      ++i;++j;
    }
  }
  return res;
}

template<class T>
inline T product_SpVSpVc(const int n1,const  int* indx1, const T* A1, const int n2, const int* indx2, const T* A2) {
  T res = T(0);
  int i=0, j=0;
  if(std::is_same<T,std::complex<RealType>>::value || std::is_same<T,std::complex<SPRealType>>::value) {
    while( i<n1 && j<n2 ) {
      if( *(indx1+i) < *(indx2+j)   )
        ++i;
      else if( *(indx2+j) < *(indx1+i) )
        ++j;
      else {
        res += *(A1+i) * (myconj(*(A2+j)));
        ++i;++j;
      }
    }
  } else {
    while( i<n1 && j<n2 ) {
      if( *(indx1+i) < *(indx2+j)   )
        ++i;
      else if( *(indx2+j) < *(indx1+i) )
        ++j;
      else {
        res += *(A1+i) * (*(A2+j));
        ++i;++j;
      }
    }
  }
  return res;
}

template<class T>
inline void transpose_SpMat(const T& A, T& AT)
{
  AT.clear();
  AT.setDims(A.cols(),A.rows());
  AT.resize_arrays(A.size());
  int n = A.size();
  std::copy(A.values(),A.values()+n,AT.values());  
  std::copy(A.row_data(),A.row_data()+n,AT.column_data());  
  std::copy(A.column_data(),A.column_data()+n,AT.row_data());  
  AT.compress();
}

bool sparseEigenSystem(ComplexSpMat &A, int& m0, RealType *eigval, ComplexType* eigVec, double Emin ); 
bool sparseEigenSystem(RealSpMat &A, int& m0, RealType *eigval, RealType* eigVec, double Emin  ); 

} // namespace SparseMatrixOperators

} // namespace qmcplusplus

#endif

