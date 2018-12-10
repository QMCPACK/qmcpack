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
// -*- C++ -*-
/**@file detemrinant.hpp
 *
 * Standalone determinant class for determinant miniapps
 */
#ifndef QMCPLUSPLUS_DETERMINANT_MINIAPPS_H
#define QMCPLUSPLUS_DETERMINANT_MINIAPPS_H
#include <OhmmsPETE/OhmmsMatrix.h>
#include <simd/allocator.hpp>
#include <Numerics/DeterminantOperators.h>

namespace qmcplusplus
{
  /**@{Determinant utilities */
  /** Inversion of a double matrix after LU factorization*/
  inline void getri(int n, double* restrict a, int lda, int* restrict piv, double* restrict work, int& lwork)
  {
    int status;
    dgetri(n,a,lda,piv,work,lwork,status);
  }

  /** Inversion of a float matrix after LU factorization*/
  inline void getri(int n, float* restrict a, int lda, int* restrict piv, float* restrict work, int& lwork)
  {
    int status;
    sgetri(n,a,lda,piv,work,lwork,status);
  }

  /** Inversion of a std::complex<double> matrix after LU factorization*/
  inline void getri(int n, std::complex<double>* restrict a, int lda, int* restrict piv, std::complex<double>* restrict work, int& lwork)
  {
    int status;
    zgetri(n,a,lda,piv,work,lwork,status);
  }

  /** Inversion of a complex<float> matrix after LU factorization*/
  inline void getri(int n, std::complex<float>* restrict a, int lda, int* restrict piv, std::complex<float>* restrict work, int& lwork)
  {
    int status;
    cgetri(n,a,lda,piv,work,lwork,status);
  }

  /** query the size of workspace for Xgetri after LU decompisiton */
  template<class T>
    inline int getGetriWorkspace(T* restrict x, int n, int lda, int* restrict pivot)
    {
      T work;
      int lwork=-1;
      getri(n, x, lda, pivot, &work, lwork);
      lwork=static_cast<int>(work);
      return lwork;
    }

  ///inner product
  template<typename T1, typename T2, typename T3>
    inline T3 inner_product_n(const T1* restrict a, const T2* restrict b, int n, T3 res)
    {
      for(int i=0; i<n; ++i) res += a[i]*b[i];
      return res;
    }

  /** transpose in to out
   *
   * Assume: in[n][lda] and out[n][lda]
   */
  template<typename TIN, typename TOUT>
    inline void transpose(const TIN* restrict in, TOUT* restrict out, int n, int lda)
    {
      for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j) 
          out[i*lda+j]=in[i+j*lda];
    }

  ///used only for debugging or walker move
  template<class T>
    inline T
    InvertWithLog(T* restrict x, int n, int lda,
        T* restrict work, int lwork, int* restrict pivot, T& phase)
    {
      T logdet(0.0);
      LUFactorization(n,n,x,lda,pivot);
      int sign_det=1;
      for(int i=0; i<n; i++)
      {
        sign_det *= (pivot[i] == i+1)?1:-1;
        sign_det *= (x[i*lda+i]>0)?1:-1;
        logdet += std::log(std::abs(x[i*lda+i]));
      }
      getri(n, x, lda, pivot, work, lwork);
      phase=(sign_det>0)?0.0:M_PI;
      return logdet;
    }

  ///recompute inverse, do not evaluate log|det|
  template<class T>
    inline void
    InvertOnly(T* restrict x, int n, int lda, T* restrict work, int* restrict pivot, int lwork)
    {
      LUFactorization(n,n,x,lda,pivot);
      getri(n, x, lda, pivot, work, lwork);
    }

  /** update Row as implemented in the full code */
  template<typename T, typename RT>
  inline void 
  inverseRowUpdate(T* restrict pinv,  const T* restrict tv, int m, int lda, int rowchanged, RT c_ratio_in)
  {
    constexpr T cone(1);
    constexpr T czero(0);
    T temp[m], rcopy[m];
    T c_ratio=cone/c_ratio_in;
    BLAS::gemv('T', m, m, c_ratio, pinv, m, tv, 1, czero, temp, 1);
    temp[rowchanged]=cone-c_ratio;
    copy_n(pinv+m*rowchanged,m,rcopy);
    BLAS::ger(m,m,-cone,rcopy,1,temp,1,pinv,m);
  }
  /**@}*/

template<typename T, typename INVT=double>
struct DiracDet
{
  ///log|det|
  INVT LogValue;
  ///current ratio
  INVT curRatio;
  ///workspace size
  int LWork;
  ///inverse matrix to be update
  Matrix<T> psiMinv; 
  ///a SPO set for the row update
  aligned_vector<T> psiV;
  ///internal storage to perform inversion correctly
  Matrix<INVT> psiM;//matrix to be inverted
  ///random number generator for testing
  RandomGenerator<T> myRandom;

  //temporary workspace for inversion
  aligned_vector<int> pivot;
  aligned_vector<INVT> work;
  Matrix<T> psiMsave;

  explicit DiracDet(int nels)
  {
    psiMinv.resize(nels,nels);
    psiV.resize(nels);
    psiM.resize(nels,nels);

    pivot.resize(nels);
    psiMsave.resize(nels,nels);
  }

  void initialize(RandomGenerator<T> RNG)
  {
    int nels=psiM.rows();
    //get lwork and resize workspace
    LWork=getGetriWorkspace(psiM.data(),nels,nels,pivot.data());
    work.resize(LWork);

    myRandom=RNG;
    constexpr T shift(0.5);
    RNG.generate_uniform(psiMsave.data(),nels*nels);
    psiMsave -= shift;

    INVT phase;
    transpose(psiMsave.data(),psiM.data(),nels,nels);
    LogValue=InvertWithLog(psiM.data(),nels,nels,work.data(),LWork,pivot.data(),phase);
    copy_n(psiM.data(),nels*nels,psiMinv.data());

    if(omp_get_num_threads()==1)
    {
      checkIdentity(psiMsave,psiM,"Psi_0 * psiM(double)");
      checkIdentity(psiMsave,psiMinv,"Psi_0 * psiMinv(T)");
      checkDiff(psiM,psiMinv,"psiM(double)-psiMinv(T)");
    }
  }

  ///recompute the inverse
  inline void recompute()
  {
    const int nels=psiV.size();
    INVT phase;
    transpose(psiMsave.data(),psiM.data(),nels,nels);
    InvertOnly(psiM.data(),nels,nels,work.data(),pivot.data(),LWork);
    copy_n(psiM.data(),nels*nels,psiMinv.data());
  }

  /** return determinant ratio for the row replacement
   * @param iel the row (active particle) index
   */
  inline INVT ratio(int iel)
  {
    const int nels=psiV.size();
    constexpr T shift(0.5);
    constexpr INVT czero(0);
    for(int j=0; j<nels; ++j) psiV[j]=myRandom()-shift;
    curRatio=inner_product_n(psiV.data(),psiMinv[iel],nels,czero);
    return curRatio;
  }

  /** accept the row and update the inverse */
  inline void accept(int iel)
  {
    const int nels=psiV.size();
    inverseRowUpdate(psiMinv.data(),psiV.data(),nels,nels,iel,curRatio);
    copy_n(psiV.data(),nels,psiMsave[iel]);
  }

  void debug()
  {
    const int nels=psiV.size();
    constexpr T shift(0.5);
    constexpr INVT czero(0.0);
    constexpr INVT eps=10.*numeric_limits<float>::epsilon();
    double ratio_error=czero;
    INVT phase;
    {
      for(int i=0; i<nels; ++i)
      {
        auto r=ratio(i);
        double ratio_full;
        if(r>0 && r>myRandom())
        {
          accept(i);
          //update A0
          copy_n(psiV.data(),nels,psiMsave[i]);
          transpose(psiMsave.data(),psiM.data(),nels,nels);
          auto  newlog=InvertWithLog(psiM.data(), nels, nels, work.data(), pivot.data(), phase);

          ratio_full=std::exp(newlog-LogValue);
          LogValue=newlog;
          double err=r/ratio_full-1;
          ratio_error += err;
#pragma omp master
          if(std::abs(err)>eps)
          {
            cout << i << " accepted: curRatio " << r << " error " << err <<endl;
            checkDiff(psiMinv,psiM,"update");
          }
        }
      }
    }
#pragma omp master
    cout << "Cummulative ratio_error " << ratio_error << endl;
  }

};

template<typename MT1, typename MT2>
void checkIdentity(const MT1& a, const MT2& b, const string& tag)
{
  constexpr double czero(0.0);
  constexpr double cone(1.0);
  const int nrows=a.rows();
  const int ncols=a.cols();
  double error=czero;
  for(int i=0; i<nrows; ++i)
  {
    for(int j=0; j<nrows; ++j)
    {
      double e=inner_product_n(a[i],b[j],ncols,czero);
      error += (i==j)? std::abs(e-cone):std::abs(e);
    }
  }
#pragma omp master
  cout << tag << " Identity Error = " << error/nrows/nrows << endl;
}

template<typename MT1, typename MT2>
void checkDiff(const MT1& a, const MT2& b, const string& tag)
{
  const int nrows=a.rows();
  const int ncols=a.cols();
  constexpr double czero(0.0);
  double error=czero;
  for(int i=0; i<nrows; ++i)
    for(int j=0; j<ncols; ++j)
      error += std::abs(static_cast<double>(a(i,j)-b(i,j)));
#pragma omp master
  cout << tag << " diff Error = " << error/nrows/nrows << endl;
}

}
#endif
