//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantCUDA.h
 * @brief Declaration of DiracDeterminantCUDA with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
#include <typeinfo>
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/Fermion/determinant_update.h"
#include "Numerics/CUDA/cuda_inverse.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{
class DiracDeterminantCUDA: public DiracDeterminantBase
{
public:
  typedef SPOSetBase::IndexVector_t IndexVector_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef SPOSetBase::GradMatrix_t  GradMatrix_t;
  typedef ParticleSet::Walker_t     Walker_t;

  DiracDeterminantCUDA(SPOSetBasePtr const &spos, int first=0);
  DiracDeterminantCUDA(const DiracDeterminantCUDA& s) = delete;

protected:
  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
  int RowStride;
  size_t AOffset, AinvOffset, newRowOffset, AinvDeltaOffset,
         AinvColkOffset, gradLaplOffset, newGradLaplOffset, 
         AWorkOffset, AinvWorkOffset;
  gpu::host_vector<CudaValueType*> UpdateList;
  gpu::device_vector<CudaValueType*> UpdateList_d;
  gpu::host_vector<updateJob> UpdateJobList;
  gpu::device_vector<updateJob> UpdateJobList_d;
  std::vector<CudaValueType*> srcList, destList, AList, AinvList, newRowList,
                              AinvDeltaList, AinvColkList, gradLaplList, newGradLaplList, 
                              AWorkList, AinvWorkList, GLList;
  gpu::device_vector<CudaValueType*> srcList_d, destList_d, AList_d, AinvList_d, newRowList_d, 
                                    AinvDeltaList_d, AinvColkList_d, gradLaplList_d, 
                                    newGradLaplList_d, AWorkList_d, AinvWorkList_d, GLList_d;
  gpu::device_vector<int> PivotArray_d;
  gpu::device_vector<int> infoArray_d;
  gpu::host_vector<int> infoArray_host;
  gpu::device_vector<CudaValueType> ratio_d;
  gpu::host_vector<CudaValueType> ratio_host;
  gpu::device_vector<CudaValueType> gradLapl_d;
  gpu::host_vector<CudaValueType> gradLapl_host;
  gpu::device_vector<int> iatList_d;
  gpu::host_vector<int> iatList;

  // Data members for nonlocal psuedopotential ratio evaluation
  static const int NLrowBufferRows = 4800;

  gpu::device_vector<CudaValueType> NLrowBuffer_d;
  gpu::host_vector<CudaValueType> NLrowBuffer_host;
  gpu::device_vector<CudaValueType*> SplineRowList_d;
  gpu::host_vector<CudaValueType*> SplineRowList_host;
  gpu::device_vector<CudaValueType*> RatioRowList_d;
  gpu::host_vector<CudaValueType*> RatioRowList_host[2];
  gpu::device_vector<CudaRealType> NLposBuffer_d;
  gpu::host_vector<CudaRealType> NLposBuffer_host;
  gpu::device_vector<CudaValueType*> NLAinvList_d;
  gpu::host_vector<CudaValueType*> NLAinvList_host[2];
  gpu::device_vector<int> NLnumRatioList_d;
  gpu::host_vector<int> NLnumRatioList_host[2];
  gpu::device_vector<int> NLelecList_d;
  gpu::host_vector<int> NLelecList_host[2];
  gpu::device_vector<CudaValueType> NLratios_d[2];
  gpu::host_vector<CudaValueType> NLratios_host;
  gpu::device_vector<CudaValueType*> NLratioList_d;
  gpu::host_vector<CudaValueType*> NLratioList_host[2];

  void resizeLists(int numWalkers)
  {
    AList.resize(numWalkers);
    AList_d.resize(numWalkers);
    AinvList.resize(numWalkers);
    AinvList_d.resize(numWalkers);
    newRowList.resize(numWalkers);
    newRowList_d.resize(numWalkers);
    AinvDeltaList.resize(numWalkers);
    AinvDeltaList_d.resize(numWalkers);
    AinvColkList.resize(numWalkers);
    AinvColkList_d.resize(numWalkers);
    ratio_d.resize(5*numWalkers);
    ratio_host.resize(5*numWalkers);
    gradLaplList.resize(numWalkers);
    gradLaplList_d.resize(numWalkers);
    GLList.resize(numWalkers);
    GLList_d.resize(numWalkers);
    newGradLaplList.resize(numWalkers);
    newGradLaplList_d.resize(numWalkers);
    AWorkList.resize(numWalkers);
    AinvWorkList.resize(numWalkers);
    AWorkList_d.resize(numWalkers);
    AinvWorkList_d.resize(numWalkers);
    iatList.resize(numWalkers);
    iatList_d.resize(numWalkers);
    // HACK HACK HACK
    // gradLapl_d.resize   (numWalkers*NumOrbitals*4);
    // gradLapl_host.resize(numWalkers*NumOrbitals*4);
    infoArray_d.resize(numWalkers*2);
    infoArray_host.resize(numWalkers*2);
    PivotArray_d.resize(numWalkers*NumOrbitals);
    gradLapl_d.resize   (numWalkers*RowStride*4);
    gradLapl_host.resize(numWalkers*RowStride*4);
    NLrowBuffer_d.resize(NLrowBufferRows*RowStride);
    NLrowBuffer_host.resize(NLrowBufferRows*RowStride);
    SplineRowList_d.resize(NLrowBufferRows);
    SplineRowList_host.resize(NLrowBufferRows);
    for (int i=0; i<NLrowBufferRows; i++)
      SplineRowList_host[i] = &(NLrowBuffer_d.data()[i*RowStride]);
    SplineRowList_d = SplineRowList_host;
    NLposBuffer_d.resize   (OHMMS_DIM * NLrowBufferRows);
    NLposBuffer_host.resize(OHMMS_DIM * NLrowBufferRows);
    for(int i = 0; i < 2; ++i)
      NLratios_d[i].resize(NLrowBufferRows);
    NLratios_host.resize(NLrowBufferRows);
  }

public:
  ValueType ratio(ParticleSet& P, int iat)
  {
    return DiracDeterminantBase::ratio (P, iat);
  }

  void update (std::vector<Walker_t*> &walkers, int iat);
  void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList);

  void reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool) {
    RowStride = ((NumOrbitals + 31)/32) * 32;
    AOffset           = pool.reserve((size_t)    NumPtcls * RowStride);
    AinvOffset        = pool.reserve((size_t)    NumPtcls * RowStride);
    gradLaplOffset    = pool.reserve((size_t)4 * NumPtcls * RowStride);
    newRowOffset      = pool.reserve((size_t)1            * RowStride);
    AinvDeltaOffset   = pool.reserve((size_t)1            * RowStride);
    AinvColkOffset    = pool.reserve((size_t)1            * RowStride);
    newGradLaplOffset = pool.reserve((size_t)4            * RowStride);
    if (typeid(CudaRealType) == typeid(float))
    {
      AWorkOffset       = pool.reserve((size_t)2 * NumPtcls * RowStride);
      AinvWorkOffset    = pool.reserve((size_t)2 * NumPtcls * RowStride);
    }
    else if (typeid(CudaRealType) == typeid(double))
    {
      AWorkOffset       = pool.reserve((size_t)    NumPtcls * RowStride);
      AinvWorkOffset    = 0;                  // not needed for inversion
    }
    Phi->reserve(pool);
  }

  void recompute (MCWalkerConfiguration &W, bool firstTime);

  void addLog (MCWalkerConfiguration &W, std::vector<RealType> &logPsi);

  void addGradient(MCWalkerConfiguration &W, int iat,
                   std::vector<GradType> &grad);

  void calcGradient(MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);

  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios);

  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad);

  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
              std::vector<ValueType> &lapl);
  void calcRatio (MCWalkerConfiguration &W, int iat,
                  std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
                  std::vector<ValueType> &lapl);
  void addRatio (MCWalkerConfiguration &W, int iat,
                 std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
                 std::vector<ValueType> &lapl);

  void ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
              std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
              std::vector<GradType>  &grad, std::vector<ValueType> &lapl);

  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);

  void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);

  void NLratios_CPU (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                     std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);
};
}
#endif // QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
