//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_MATRIX_H
#define QMCPLUSPLUS_DIRAC_MATRIX_H

#include "Numerics/Blasf.h"
#include <OhmmsPETE/OhmmsMatrix.h>
#include <type_traits/scalar_traits.h>

namespace qmcplusplus { 

  inline void Xgetrf(int n, int m, float* restrict a, int lda, int* restrict piv)
  {
    int status;
    sgetrf(n,m,a,lda,piv,status);
  }

  inline void Xgetri(int n, float* restrict a, int lda, int* restrict piv, float* restrict work, int& lwork)
  {
    int status;
    sgetri(n,a,lda,piv,work,lwork,status);
  }

  inline void Xgetrf(int n, int m, std::complex<float>* restrict a, int lda, int* restrict piv)
  {
    int status;
    cgetrf(n,m,a,lda,piv,status);
  }

  /** inversion of a float matrix after lu factorization*/
  inline void Xgetri(int n, std::complex<float>* restrict a, int lda, int* restrict piv, std::complex<float>* restrict work, int& lwork)
  {
    int status;
    cgetri(n,a,lda,piv,work,lwork,status);
  }

  inline void Xgetrf(int n, int m, double* restrict a, int lda, int* restrict piv)
  {
    int status;
    dgetrf(n,m,a,lda,piv,status);
  }

  inline void Xgetri(int n, double* restrict a, int lda, int* restrict piv, double* restrict work, int& lwork)
  {
    int status;
    dgetri(n,a,lda,piv,work,lwork,status);
  }

  inline void Xgetrf(int n, int m, std::complex<double>* restrict a, int lda, int* restrict piv)
  {
    int status;
    zgetrf(n,m,a,lda,piv,status);
  }

  /** inversion of a std::complex<double> matrix after lu factorization*/
  inline void Xgetri(int n, std::complex<double>* restrict a, int lda, int* restrict piv, std::complex<double>* restrict work, int& lwork)
  {
    int status;
    zgetri(n,a,lda,piv,work,lwork,status);
  }


  template<typename TIN, typename TOUT>
    inline void TansposeSquare(const TIN* restrict in, TOUT* restrict out, size_t n, size_t lda)
    {
#pragma omp simd 
      for(size_t i=0; i<n; ++i)
        for(size_t j=0; j<n; ++j) 
          out[i*lda+j]=in[i+j*lda];
    }


  template<typename T>
    inline T computeLogDet(const T* restrict X, int n, int lda, const int* restrict pivot, T& phase)
    {
      T logdet(0);
      int sign_det=1;
      for(size_t i=0; i<n; i++)
      {
        const size_t ii=i*lda+i;
        sign_det *= (pivot[i] == i+1)?1:-1;
        sign_det *= (X[ii]>0)?1:-1;
        logdet += std::log(std::abs(X[ii]));
      }
      phase=(sign_det>0)? T(0):M_PI;
      return logdet;
    }

  template<typename T>
    inline T computeLogDet(const std::complex<T>* restrict X, int n, int lda, const int* restrict pivot, T& phase)
    {
      T logdet(0);
      phase=T(0);
      for(size_t i=0; i<n; i++)
      {
        const size_t ii=i*lda+i;
        phase += std::arg(X[ii]);
        if(pivot[i]!=i+1)
          phase += M_PI;
        logdet+=std::log(X[ii].real()*X[ii].real()+X[ii].imag()*X[ii].imag());
        //slightly smaller error with the following
        //        logdet+=2.0*std::log(std::abs(x[ii]);
      }
      CONSTEXPR T one_over_2pi=T(1)/TWOPI;
      phase -= std::floor(phase*one_over_2pi)*TWOPI;
      return 0.5*logdet;
    }

  template<typename T>
    struct DiracMatrix
    {
      typedef typename scalar_traits<T>::real_type real_type;
      aligned_vector<T> m_work;
      aligned_vector<int> m_pivot;
      int Lwork;
      real_type LogDet;
      real_type Phase;

      DiracMatrix():Lwork(0) {}

      inline void invert(Matrix<T>& amat, bool computeDet)
      {
        const int n=amat.rows();
        const int lda=amat.cols();
        if(Lwork<n) reset(amat,n,lda);
        int status;
        Xgetrf(n,n,amat.data(),lda,m_pivot.data());
        if(computeDet)
        {
          LogDet=computeLogDet(amat.data(),n,lda,m_pivot.data(),Phase);
        }
        Xgetri(n,  amat.data(),lda,m_pivot.data(),m_work.data(),Lwork);
      }

      inline void reset(Matrix<T>& amat,int n, int lda)
      {
        //Lwork=n;
        //m_work.resize(n);
        //m_pivot.resize(n);
        m_pivot.resize(lda);
        Lwork=-1;
        T tmp;
        real_type lw;
        Xgetri(n, amat.data(),lda,m_pivot.data(),&tmp,Lwork);
        convert(tmp,lw);
        Lwork=static_cast<int>(lw);
        m_work.resize(Lwork); 
      }

      inline void updateRow(Matrix<T>& a, T* arow, int rowchanged, T c_ratio_in)
      {
        const int m=a.rows();
        const int lda=a.cols();
        CONSTEXPR T cone(1);
        CONSTEXPR T czero(0);
        T temp[lda], rcopy[lda];
        T c_ratio=cone/c_ratio_in;
        BLAS::gemv('T', m, m, c_ratio, a.data(), lda, arow, 1, czero, temp, 1);
        temp[rowchanged]=cone-c_ratio;
        simd::copy_n(a[rowchanged],m,rcopy);
        BLAS::ger(m,m,-cone,rcopy,1,temp,1,a.data(),lda);
      }
    };
}

#endif // OHMMS_PETE_MATRIX_H

