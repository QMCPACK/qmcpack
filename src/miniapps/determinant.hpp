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
#ifndef QMCPLUSPLUS_DETERMINANT_MINIAPPS_H
#define QMCPLUSPLUS_DETERMINANT_MINIAPPS_H
#include <OhmmsPETE/OhmmsMatrix.h>
#include <simd/simd.hpp>
#include <simd/allocator.hpp>
#include <simd/algorithm.hpp>
#include <Numerics/DeterminantOperators.h>
namespace qmcplusplus
{
template<typename T, typename INVT=double>
struct DiracDet
{
  ///log|det|
  INVT LogValue;
  ///current ratio
  INVT curRatio;
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
  aligned_vector<T> workV1,workV2;
  Matrix<T> psiMsave;

  explicit DiracDet(int nels)
  {
    psiMinv.resize(nels,nels);
    psiV.resize(nels);
    psiM.resize(nels,nels);

    pivot.resize(nels);
    work.resize(nels);
    workV1.resize(nels);
    workV2.resize(nels);

    psiMsave.resize(nels,nels);
  }

  void initialize(RandomGenerator<T> RNG)
  {
    myRandom=RNG;
    constexpr T shift(0.5);
    int nels=psiM.rows();
    RNG.generate_uniform(psiMsave.data(),nels*nels);
    psiMsave -= shift;

    recompute();

    if(omp_get_num_threads()==1)
    {
      checkIdentity(psiMsave,psiM,"Psi_0 * psiM(double)");
      checkIdentity(psiMsave,psiMinv,"Psi_0 * psiMinv(T)");
      checkDiff(psiM,psiMinv,"psiM(double)-psiMinv(T)");
    }
  }

  void recompute()
  {
    const int nels=psiV.size();
    INVT phase;
    simd::transpose(psiMsave.data(),psiM.data(),nels,nels);
    LogValue=InvertWithLog(psiM.data(),nels,nels,work.data(),pivot.data(),phase);
    simd::copy_n(psiM.data(),nels*nels,psiMinv.data());
  }

  inline INVT ratio(int iel)
  {
    const int nels=psiV.size();
    constexpr T shift(0.5);
    constexpr INVT czero(0);
    for(int j=0; j<nels; ++j) psiV[j]=myRandom()-shift;
    curRatio=simd::inner_product_n(psiV.data(),psiMinv[iel],nels,czero);
    return curRatio;
  }

  inline void accept(int iel)
  {
    InverseUpdateByRow(psiMinv,psiV,workV1,workV2,iel,curRatio);
    simd::copy_n(psiV.data(),psiV.size(),psiMsave[iel]);
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
          simd::copy_n(psiV.data(),nels,psiMsave[i]);
          simd::transpose(psiMsave.data(),psiM.data(),nels,nels);
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
      double e=simd::inner_product_n(a[i],b[j],ncols,czero);
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
