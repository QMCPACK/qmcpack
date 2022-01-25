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
 * @brief Declaration of DiracDeterminantCUDA with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
#include <typeinfo>
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/determinant_update.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{
class DiracDeterminantCUDA : public DiracDeterminantBase
{
public:
  using IndexVector = SPOSet::IndexVector;
  using ValueVector = SPOSet::ValueVector;
  using ValueMatrix = SPOSet::ValueMatrix;
  using GradVector  = SPOSet::GradVector;
  using GradMatrix  = SPOSet::GradMatrix;
  using Walker_t    = ParticleSet::Walker_t;

  DiracDeterminantCUDA(std::unique_ptr<SPOSet>&& spos, int first, int last);
  DiracDeterminantCUDA(const DiracDeterminantCUDA& s) = delete;

protected:
  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
  int RowStride;
  size_t AOffset, AinvOffset, LemmaOffset, LemmaLUOffset, LemmaInvOffset, AinvUOffset, newRowOffset, AinvDeltaOffset,
      AinvColkOffset, gradLaplOffset, newGradLaplOffset, AWorkOffset, AinvWorkOffset;
  gpu::host_vector<CTS::ValueType*> UpdateList;
  gpu::device_vector<CTS::ValueType*> UpdateList_d;
  gpu::host_vector<updateJob> UpdateJobList;
  gpu::device_vector<updateJob> UpdateJobList_d;
  std::vector<CTS::ValueType*> srcList, destList, AList, AinvList, newRowList, LemmaList, LemmaLUList, LemmaInvList,
      AinvUList, AinvDeltaList, AinvColkList, gradLaplList, newGradLaplList, AWorkList, AinvWorkList, GLList;
  gpu::device_vector<CTS::ValueType*> srcList_d, destList_d, AList_d, AinvList_d, newRowList_d, LemmaList_d,
      LemmaLUList_d, LemmaInvList_d, AinvUList_d, AinvDeltaList_d, AinvColkList_d, gradLaplList_d, newGradLaplList_d,
      AWorkList_d, AinvWorkList_d, GLList_d;
  gpu::device_vector<int> PivotArray_d;
  gpu::device_vector<int> infoArray_d;
  gpu::host_vector<int> infoArray_host;
  gpu::device_vector<CTS::ValueType> ratio_d;
  gpu::host_vector<CTS::ValueType> ratio_host;
  gpu::device_vector<CTS::ValueType> gradLapl_d;
  gpu::host_vector<CTS::ValueType> gradLapl_host;
  gpu::device_vector<int> iatList_d;
  gpu::host_vector<int> iatList;

  // Data members for nonlocal psuedopotential ratio evaluation
  static const int NLrowBufferRows = 4800;

  gpu::device_vector<CTS::ValueType> NLrowBuffer_d;
  gpu::host_vector<CTS::ValueType> NLrowBuffer_host;
  gpu::device_vector<CTS::ValueType*> SplineRowList_d;
  gpu::host_vector<CTS::ValueType*> SplineRowList_host;
  gpu::device_vector<CTS::ValueType*> RatioRowList_d;
  gpu::host_vector<CTS::ValueType*> RatioRowList_host[2];
  gpu::device_vector<CTS::RealType> NLposBuffer_d;
  gpu::host_vector<CTS::RealType> NLposBuffer_host;
  gpu::device_vector<CTS::ValueType*> NLAinvList_d;
  gpu::host_vector<CTS::ValueType*> NLAinvList_host[2];
  gpu::device_vector<int> NLnumRatioList_d;
  gpu::host_vector<int> NLnumRatioList_host[2];
  gpu::device_vector<int> NLelecList_d;
  gpu::host_vector<int> NLelecList_host[2];
  gpu::device_vector<CTS::ValueType> NLratios_d[2];
  gpu::host_vector<CTS::ValueType> NLratios_host;
  gpu::device_vector<CTS::ValueType*> NLratioList_d;
  gpu::host_vector<CTS::ValueType*> NLratioList_host[2];

  void resizeLists(int numWalkers) { resizeLists(numWalkers, 1); }

  void resizeLists(int numWalkers, int kdelay)
  {
    AList.resize(numWalkers);
    AList_d.resize(numWalkers);
    AinvList.resize(numWalkers * kdelay);
    AinvList_d.resize(numWalkers * kdelay);
    AinvUList.resize(numWalkers);
    AinvUList_d.resize(numWalkers);
    newRowList.resize(numWalkers * kdelay);
    newRowList_d.resize(numWalkers * kdelay);
    AinvDeltaList.resize(numWalkers);
    AinvDeltaList_d.resize(numWalkers);
    AinvColkList.resize(numWalkers);
    AinvColkList_d.resize(numWalkers);
    LemmaList.resize(numWalkers);
    LemmaList_d.resize(numWalkers);
    LemmaLUList.resize(numWalkers);
    LemmaLUList_d.resize(numWalkers);
    LemmaInvList.resize(numWalkers);
    LemmaInvList_d.resize(numWalkers);
    ratio_d.resize(5 * numWalkers);
    ratio_host.resize(5 * numWalkers);
    gradLaplList.resize(numWalkers);
    gradLaplList_d.resize(numWalkers);
    GLList.resize(numWalkers);
    GLList_d.resize(numWalkers);
    newGradLaplList.resize(numWalkers * kdelay);
    newGradLaplList_d.resize(numWalkers * kdelay);
    AWorkList.resize(numWalkers);
    AinvWorkList.resize(numWalkers);
    AWorkList_d.resize(numWalkers);
    AinvWorkList_d.resize(numWalkers);
    iatList.resize(numWalkers);
    iatList_d.resize(numWalkers);
    // HACK HACK HACK
    // gradLapl_d.resize   (numWalkers*NumOrbitals*4);
    // gradLapl_host.resize(numWalkers*NumOrbitals*4);
    infoArray_d.resize(numWalkers * 2);
    infoArray_host.resize(numWalkers * 2);
    PivotArray_d.resize(numWalkers * NumOrbitals);
    gradLapl_d.resize(numWalkers * RowStride * 4);
    gradLapl_host.resize(numWalkers * RowStride * 4);
    NLrowBuffer_d.resize(NLrowBufferRows * RowStride);
    NLrowBuffer_host.resize(NLrowBufferRows * RowStride);
    SplineRowList_d.resize(NLrowBufferRows);
    SplineRowList_host.resize(NLrowBufferRows);
    for (int i = 0; i < NLrowBufferRows; i++)
      SplineRowList_host[i] = &(NLrowBuffer_d.data()[i * RowStride]);
    SplineRowList_d = SplineRowList_host;
    NLposBuffer_d.resize(OHMMS_DIM * NLrowBufferRows);
    NLposBuffer_host.resize(OHMMS_DIM * NLrowBufferRows);
    for (int i = 0; i < 2; ++i)
      NLratios_d[i].resize(NLrowBufferRows);
    NLratios_host.resize(NLrowBufferRows);
  }

public:
  // safe-guard all CPU interfaces
  std::unique_ptr<DiracDeterminantBase> makeCopy(std::unique_ptr<SPOSet>&& spo) const override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::makeCopy is illegal!");
    return nullptr;
  }

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::evaluateLog is illegal!");
    return 0;
  }

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::acceptMove is illegal!");
  }

  void restore(int iat) override { throw std::runtime_error("Calling DiracDeterminantCUDA::restore is illegal!"); }

  PsiValueType ratio(ParticleSet& P, int iat) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::ratio is illegal!");
    return 0;
  }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::ratioGrad is illegal!");
    return 0;
  }

  void registerData(ParticleSet& P, WFBufferType& buf) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::registerData is illegal!");
  }

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::updateBuffer is illegal!");
    return 0;
  }

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override
  {
    throw std::runtime_error("Calling DiracDeterminantCUDA::copyFromBuffer is illegal!");
  }

  void update(MCWalkerConfiguration* W,
              std::vector<Walker_t*>& walkers,
              int iat,
              std::vector<bool>* acc,
              int k) override;
  void update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iatList) override;

  void reserve(PointerPool<gpu::device_vector<CTS::ValueType>>& pool, int kblocksize = 1) override
  {
    RowStride         = ((NumOrbitals + 31) / 32) * 32;
    size_t kblock2    = ((kblocksize * kblocksize + 31) / 32) * 32;
    AOffset           = pool.reserve((size_t)NumPtcls * RowStride);
    AinvOffset        = pool.reserve((size_t)NumPtcls * RowStride);
    LemmaOffset       = pool.reserve((size_t)kblock2);
    LemmaLUOffset     = pool.reserve((size_t)kblock2);
    LemmaInvOffset    = pool.reserve((size_t)kblock2);
    AinvUOffset       = pool.reserve((size_t)1 * RowStride * kblocksize);
    gradLaplOffset    = pool.reserve((size_t)4 * NumPtcls * RowStride);
    newRowOffset      = pool.reserve((size_t)1 * RowStride * kblocksize);
    AinvDeltaOffset   = pool.reserve((size_t)1 * RowStride);
    AinvColkOffset    = pool.reserve((size_t)1 * RowStride);
    newGradLaplOffset = pool.reserve((size_t)4 * RowStride * kblocksize);
    if (typeid(CTS::RealType) == typeid(float))
    {
      AWorkOffset    = pool.reserve((size_t)2 * NumPtcls * RowStride);
      AinvWorkOffset = pool.reserve((size_t)2 * NumPtcls * RowStride);
    }
    else if (typeid(CTS::RealType) == typeid(double))
    {
      AWorkOffset    = pool.reserve((size_t)NumPtcls * RowStride);
      AinvWorkOffset = 0;
    }
    Phi->reserve(pool);
  }

  void recompute(MCWalkerConfiguration& W, bool firstTime) override;

  void addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi) override;

  void addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad) override;

  void calcGradient(MCWalkerConfiguration& W, int iat, int k, std::vector<GradType>& grad) override;

  void ratio(MCWalkerConfiguration& W, int iat, std::vector<ValueType>& psi_ratios) override;

  void ratio(MCWalkerConfiguration& W,
             int iat,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad) override;

  void ratio(MCWalkerConfiguration& W,
             int iat,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad,
             std::vector<ValueType>& lapl) override;
  void calcRatio(MCWalkerConfiguration& W,
                 int iat,
                 std::vector<ValueType>& psi_ratios,
                 std::vector<GradType>& grad,
                 std::vector<ValueType>& lapl) override;
  void addRatio(MCWalkerConfiguration& W,
                int iat,
                int k,
                std::vector<ValueType>& psi_ratios,
                std::vector<GradType>& grad,
                std::vector<ValueType>& lapl) override;

  void ratio(std::vector<Walker_t*>& walkers,
             std::vector<int>& iatList,
             std::vector<PosType>& rNew,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& grad,
             std::vector<ValueType>& lapl) override;

  void gradLapl(MCWalkerConfiguration& W, GradMatrix& grads, ValueMatrix& lapl) override;

  void NLratios(MCWalkerConfiguration& W,
                std::vector<NLjob>& jobList,
                std::vector<PosType>& quadPoints,
                std::vector<ValueType>& psi_ratios) override;

  void NLratios_CPU(MCWalkerConfiguration& W,
                    std::vector<NLjob>& jobList,
                    std::vector<PosType>& quadPoints,
                    std::vector<ValueType>& psi_ratios);

  void det_lookahead(MCWalkerConfiguration& W,
                     std::vector<ValueType>& psi_ratios,
                     std::vector<GradType>& grad,
                     std::vector<ValueType>& lapl,
                     int iat,
                     int k,
                     int kd,
                     int nw) override;
};
} // namespace qmcplusplus
#endif // QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
