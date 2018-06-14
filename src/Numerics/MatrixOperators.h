//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_MATRIXOPERATORS_H
#define QMCPLUSPLUS_MATRIXOPERATORS_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "Numerics/OhmmsBlas.h"
#include <simd/simd.hpp>

namespace qmcplusplus
{

template<typename T>
inline T TESTDOT(const T* restrict f, const T* restrict l, const T* restrict b)
{
  T s=0;
  while(f!=l)
    s += (*f++)*(*b++);
  return s;
}

namespace MatrixOperators
{
  /** static function to perform C=AB for real matrices
   *
   * Call dgemm/zgemm/sgemm/cgemm via BLAS::gemm
   */
  template<typename T>
  inline void product(const Matrix<T>& A,
                      const Matrix<T>& B, Matrix<T>& C)
  {
    const char transa = 'N';
    const char transb = 'N';
    const T one(1.0);
    const T zero(0.0);
    BLAS::gemm(transa, transb, B.cols(), A.rows(), B.rows(),
          one, B.data(), B.cols(), A.data(), A.cols(),
          zero, C.data(), C.cols());
  }

  template<typename T>
  inline void product_ABt(const Matrix<T>& A,
      const Matrix<T>& B, Matrix<T>& C)
  {
    const char transa = 't';
    const char transb = 'n';
    const T one(1.0);
    const T zero(0.0);
    BLAS::gemm(transa, transb, B.rows(), A.rows(), B.cols(),
        one, B.data(), B.cols(), A.data(), A.cols(),
        zero, C.data(), C.cols());
  }

  template<typename T>
  inline void product_AtB(const Matrix<T>& A,
      const Matrix<T>& B, Matrix<T>& C)
  {
    const char transa = 'n';
    const char transb = 't';
    const T one(1.0);
    const T zero(0.0);
    BLAS::gemm(transa, transb, B.cols(), A.cols(), B.rows(),
        one, B.data(), B.cols(), A.data(), A.cols(),
        zero, C.data(), C.cols());
  }


  inline void half_outerProduct(const Matrix<double> &M,
                                const Vector<double> &B,
                                int iat,
                                Matrix<double> &C)
  {
    for (int i=0; i<C.rows(); i++)
      for (int j=0; j<C.cols(); j++)
        C(iat,i)+=M(i,j)*B[j];
    //for (int i=0; i<C.rows(); i++)
    //  C(iat,i)+=simd::dot(M[i],B.data(),C.cols());
  }

  inline void other_half_outerProduct(const Matrix<double> &M,
                                      const Vector<double> &B,
                                      int iat,
                                      Matrix<double> &C)
  {
    for (int i=0; i<C.rows(); i++)
      for (int j=0; j<C.cols(); j++)
        C(i,iat)+=M(i,j)*B[j];
    //for (int i=0; i<C.rows(); i++)
    //  C(i,iat)+=simd::dot(M[i],B.data(),C.cols());
  }

  inline double  innerProduct(const Vector<double> &A,
                              const Vector<double> &B)
  {
    double tot=0.0;
    for (int i=0; i<A.size(); i++)
      tot+=A[i]*B[i];
    return tot;
    //return simd::dot(A.data(),B.data(),A.size());
  }



  template<typename T>
    inline void transpose(Matrix<T>& A)
  {
    for (int i=0; i<A.extent(0); i++)
      for (int j=0; j<i; j++)
        std::swap(A(i,j),A(j,i));
  }

  template<typename T>
    inline void transpose(const Matrix<T>& A, Matrix<T>& B)
  {
    simd::transpose(A.data(), A.rows(), A.cols(), B.data(), B.rows(), B.cols());
  }


  /// C = A*diag(B)
  template<typename T1,typename T2,typename T3>
  inline void diag_product(const Matrix<T1>& A,
                           const Vector<T2>& B,
                           Matrix<T3>&       C)
  {
    for (int i=0; i<C.rows(); ++i)
      for (int j=0; j<C.cols(); ++j)
        C(i,j)=A(i,j)*B[j];
    //works?
    //const int ccols = C.cols();
    //const int ijmax = C.size();
    //for (int ij=0; ij<ijmax; ++ij)
    //  C(ij)=A(ij)*B(ij%ccols);
  }


  /// C = diag(A)*B
  template<typename T1,typename T2,typename T3>
  inline void diag_product(const Vector<T1> &A,
                           const Matrix<T2> &B,
                           Matrix<T3> &C)
  {

    for (int i=0; i<C.rows(); ++i)
      for (int j=0; j<C.cols(); ++j)
        C(i,j)=A[i]*B(i,j);


    //const int crows = C.rows();
    //const int ccols = C.cols();
    //for (int i=0,ijb=0; i<crows; ++i,ijb+=ccols)
    //{
    //  const T1 a = A(i);
    //  for (int j=0; j<ccols; ++j)
    //  {
    //    int ij = ijb+j;
    //    C(ij)=a*B(ij);
    //  }
    //}
    //
    //const int crows = C.rows();
    //const int ijmax = C.size();
    //for (int ij=0; ij<ijmax; ++ij)
    //  C(ij)=A(ij%crows)*B(ij);

  }


  /** static function to perform C=AB for complex matrices
   *
   * Call zgemm
   */
  inline void product(const Matrix<double>& A,
                      const Matrix<std::complex<double> >& B, Matrix<double>& C)
  {
    std::cerr << " Undefined C=AB with real A and complex B " << std::endl;
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  template<typename T>
  inline void product(const Matrix<T>& A, const Vector<T>& x, T* restrict yptr)
  {
    const char transa = 'T';
    const T one=1.0;
    const T zero=0.0;
    BLAS::gemv(transa, A.cols(), A.rows(), one, A.data(), A.cols(), x.data(), 1, zero, yptr, 1);
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  template<typename T>
  inline void product(const Matrix<T>& A, const T* restrict xptr, T* restrict yptr)
  {
    const char transa = 'T';
    const T one=1.0;
    const T zero=0.0;
    BLAS::gemv(transa, A.cols(), A.rows(), one, A.data(), A.cols(), xptr, 1, zero, yptr, 1);
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  template<typename T, unsigned D>
    inline void product(const Matrix<T>& A, const TinyVector<T,D>* xvPtr,
                        TinyVector<T,D>* restrict yptr)
  {
    const T one=1.0;
    const T zero=0.0;
    const char transa = 'N';
    const char transb = 'N';
    BLAS::gemm(transa, transb, D, A.rows(), A.cols(),
          one, xvPtr->begin(), D, A.data(), A.cols(),
          zero, yptr->begin(), D);
  }

  template<typename T, unsigned D>
    inline void product(const Matrix<T>& A, const Tensor<T,D>* xvPtr,
                        Tensor<T,D>* restrict yptr)
  {
    const T one=1.0;
    const T zero=0.0;
    const char transa = 'N';
    const char transb = 'N';
    BLAS::gemm(transa, transb, D*D, A.rows(), A.cols(),
          one, xvPtr->begin(), D*D, A.data(), A.cols(),
          zero, yptr->begin(), D*D);
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  template<unsigned D>
    inline void product(const Matrix<double>& A, const Vector<TinyVector<double,D> >& x,
                        TinyVector<double,D>* restrict yptr)
  {
    const double one=1.0;
    const double zero=0.0;
    const char transa = 'N';
    const char transb = 'N';
    dgemm(transa, transb, D, A.rows(), x.size(),
          one, x.data()->begin(), D, A.data(), A.cols(),
          zero, yptr->begin(), D);
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  template<unsigned D>
    inline void product(const Matrix<std::complex<double> >& A,
                        const Vector<TinyVector<std::complex<double>,D> >& x,
                        TinyVector<std::complex<double>,D>* restrict yptr)
  {
    const char transa = 'N';
    const char transb = 'N';
    const std::complex<double> zone(1.0,0.0);
    const std::complex<double> zero(0.0,0.0);
    zgemm(transa, transb, D, A.rows(), x.size(),
          zone, x.data()->begin(), D, A.data(), A.cols(),
          zero, yptr->begin(), D);
  }


  /** static function to perform y=Ax for generic matrix and vector
   */
  inline void product(const Matrix<std::complex<double> >& A,
                      const Vector<std::complex<double> >& x,
                      std::complex<double>* restrict yptr)
  {
    const char transa = 'T';
    const std::complex<double> zone(1.0,0.0);
    const std::complex<double> zero(0.0,0.0);
    zgemv(transa, A.cols(), A.rows(), zone, A.data(), A.cols(), x.data(), 1, zero, yptr, 1);
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  inline void product(const Matrix<std::complex<double> >& A
                      , const Vector<double>& x, std::complex<double>* restrict yptr)
  {
    const int n=A.rows();
    const int m=A.cols();
    const std::complex<double>* restrict aptr=A.data();
    for(int i=0,ij=0; i<n; ++i)
    {
      std::complex<double> t=0.0;
      for(int j=0; j<m; ++j,++ij)
        t += aptr[ij]*x[j];
      yptr[i]=t;
    }
  }

  inline void product(const Matrix<std::complex<double> >& A
                      , const std::complex<double>* restrict x
                      , std::complex<double>* restrict yptr)
  {
    const char transa = 'T';
    const std::complex<double> zone(1.0,0.0);
    const std::complex<double> zero(0.0,0.0);
    zgemv(transa, A.cols(), A.rows(), zone, A.data(), A.cols(), x, 1, zero, yptr, 1);
  }

  /** static function to perform y=Ax for generic matrix and vector
   */
  inline void product(const Matrix<double>& A,
                      const Vector<std::complex<double> >& x, double* restrict yptr)
  {
    std::cerr << " Undefined C=AB with real A and complex x " << std::endl;
  }

  template<typename T>
    inline void product(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, std::vector<int>& offset)
  {
    int nr=C.rows();
    int nb=offset.size()-1;
    for(int i=0; i<nr; i++)
    {
      for(int b=0; b<nb; b++)
      {
        int firstK=offset[b];
        int lastK=offset[b+1];
        const T* restrict firstY=A[i]+firstK;
        const T* restrict lastY=A[i]+lastK;
        for(int k=firstK; k<lastK; k++)
        {
          C(i,k)=TESTDOT(firstY,lastY,B[k]+firstK);
        }
      }
    }
  }

//    template<typename T>
//      inline statis void product(const Matrix<T>& A, const T* restrict x, T* restrict y)
//      {
//        GEMV<T,0>::apply(A.data(),x,y,A.rows(),A.cols());
//      }
}

/** API to handle gemv */
namespace simd
{
  template<typename T>
    inline void gemv(const Matrix<T>& a, const T* restrict v, T* restrict b)
    {
      MatrixOperators::product(a,v,b);
    }

  template<typename T,unsigned D>
    inline void gemv(const Matrix<T>& a, const TinyVector<T,D>* restrict v, TinyVector<T,D>* restrict b)
    {
      MatrixOperators::product(a,v,b);
    }

  template<typename T,unsigned D>
    inline void gemv(const Matrix<T>& a, const Tensor<T,D>* restrict v, Tensor<T,D>* restrict b)
    {
      MatrixOperators::product(a,v,b);
    }

}

}
#endif
