#ifndef AFQMC_SPARSEMATRIXOPERATORS_H
#define AFQMC_SPARSEMATRIXOPERATORS_H

#include<iterator>
#include<tuple>

#include "AFQMC/config.h"

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
void product_SD(const IndexType K,
             const s2D<ComplexType>* A, const int nterms,
             ComplexType* B, const IndexType LDB,  
             ComplexType* C, IndexType LDC );

template<class T>
void product_SpMatV(int nrows,
             T& A, 
             ComplexType* B, 
             ComplexType* C );

template<class T>
void product_SpMatV(const int M, const int K,
             const ComplexType& alpha,
             const T& A,
             const ComplexType* B,
             const ComplexType& beta, 
             ComplexType* C );

template<class T>
void product_SpMatV(const int M, const int K,
             const RealType& alpha,
             const T& A,
             const RealType* B,
             const RealType& beta,
             RealType* C );

void product_SpMatV(const int M, const int K,
             const RealType& alpha,
             const RealType* val,
             const int* col,
             const int* row,
             const RealType* B,
             const RealType& beta,
             RealType* C );

void product_SpMatV(const int M, const int K,
             const ComplexType& alpha,
             const ComplexType* val,
             const int* col,
             const int* row,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C );

void product_SpMatTV(const int M, const int K,
             const ComplexType& alpha,
             const ComplexType* val,
             const int* col,
             const int* row,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C );

template<class T>
void product_SpMatTV(const int M, const int K,
             const ComplexType& alpha,
             const T& A,
             const ComplexType* B,
             const ComplexType& beta, 
             ComplexType* C );

template<class T>
void product_SpMatM(const int M, const int N, const int K,
             const ComplexType& alpha,
             const T& A,
             const ComplexType* B, const int ldb,
             const ComplexType& beta, 
             ComplexType* C, int ldc );

template<class T>
void product_SpMatM(const int M, const int N, const int K,
             const RealType& alpha,
             const T& A,
             const RealType* B, const int ldb,
             const RealType& beta, 
             RealType* C, int ldc );

template<class T>
void product_SpMatM(const int M, const int N, const int K,
             const float& alpha,
             const T& A,
             const float* B, const int ldb,
             const float& beta,
             float* C, int ldc );



// Performs a product between a sparse matrix stored in format s2D and a dense 
// matrix stored in c format
//   N: number of rows in B/C 
//   M: number of columns in B/C
//   LDB: leading dimension of B
//
//   For eack term in A, aik 
//   C(i,:) += aik * B(k.:) 
void product_SD(const IndexType K,
             const s2D<RealType>* A, const int nterms,
             ComplexType* B, const IndexType LDB,                                            
             ComplexType* C, IndexType LDC );

// Dot product between 2 sparse vectors
template<class T>
T product_SpVSpV(const int n1, const int* indx1, const T* A1, const int n2, const int* indx2, const T* A2);

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

