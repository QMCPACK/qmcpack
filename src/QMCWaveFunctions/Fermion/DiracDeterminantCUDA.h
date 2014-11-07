//////////////////////////////////////////////////////////////////
// (c) Copyright 2009- by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Ken Esler
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file DiracDeterminantCUDA.h
 * @brief Declaration of DiracDeterminantCUDA with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
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
  DiracDeterminantCUDA(const DiracDeterminantCUDA& s);

protected:
  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
  int RowStride;
  size_t AOffset, AinvOffset, newRowOffset, AinvDeltaOffset,
         AinvColkOffset, gradLaplOffset, newGradLaplOffset, 
         AWorkOffset, AinvWorkOffset;
  gpu::host_vector<CudaRealType*> UpdateList;
  gpu::device_vector<CudaRealType*> UpdateList_d;
  gpu::host_vector<updateJob> UpdateJobList;
  gpu::device_vector<updateJob> UpdateJobList_d;
  vector<CudaRealType*> srcList, destList, AList, AinvList, newRowList,
                        AinvDeltaList, AinvColkList, gradLaplList, newGradLaplList, 
                        AWorkList, AinvWorkList, GLList;
  gpu::device_vector<CudaRealType*> srcList_d, destList_d, AList_d, AinvList_d, newRowList_d, 
                                    AinvDeltaList_d, AinvColkList_d, gradLaplList_d, 
                                    newGradLaplList_d, AWorkList_d, AinvWorkList_d, GLList_d;
  gpu::device_vector<CudaRealType> ratio_d;
  gpu::host_vector<CudaRealType> ratio_host;
  gpu::device_vector<CudaRealType> gradLapl_d;
  gpu::host_vector<CudaRealType> gradLapl_host;
  gpu::device_vector<int> iatList_d;
  gpu::host_vector<int> iatList;

  // Data members for nonlocal psuedopotential ratio evaluation
  static const int NLrowBufferRows = 4800;
  gpu::device_vector<CudaRealType> NLrowBuffer_d;
  gpu::host_vector<CudaRealType> NLrowBuffer_host;

  gpu::device_vector<CudaRealType*> SplineRowList_d;
  gpu::host_vector<CudaRealType*> SplineRowList_host;
  gpu::device_vector<CudaRealType*> RatioRowList_d;
  gpu::host_vector<CudaRealType*> RatioRowList_host[2];
  gpu::device_vector<CudaRealType> NLposBuffer_d;
  gpu::host_vector<CudaRealType> NLposBuffer_host;
  gpu::device_vector<CudaRealType*> NLAinvList_d;
  gpu::host_vector<CudaRealType*> NLAinvList_host[2];
  gpu::device_vector<int> NLnumRatioList_d;
  gpu::host_vector<int> NLnumRatioList_host[2];
  gpu::device_vector<int> NLelecList_d;
  gpu::host_vector<int> NLelecList_host[2];
  gpu::device_vector<CudaRealType> NLratios_d[2];
  gpu::host_vector<CudaRealType> NLratios_host;
  gpu::device_vector<CudaRealType*> NLratioList_d;
  gpu::host_vector<CudaRealType*> NLratioList_host[2];

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

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {
    return DiracDeterminantBase::ratio (P, iat, dG, dL);
  }

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL, int iat)
  {
    DiracDeterminantBase::update (P, dG, dL, iat);
  }





  void update (vector<Walker_t*> &walkers, int iat);
  void update (const vector<Walker_t*> &walkers, const vector<int> &iatList);

  void reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool)
  {
    RowStride = ((NumOrbitals + 31)/32) * 32;
    AOffset           = pool.reserve((size_t)    NumPtcls * RowStride);
    AinvOffset        = pool.reserve((size_t)    NumPtcls * RowStride);
    gradLaplOffset    = pool.reserve((size_t)4 * NumPtcls * RowStride);
    newRowOffset      = pool.reserve((size_t)1            * RowStride);
    AinvDeltaOffset   = pool.reserve((size_t)1            * RowStride);
    AinvColkOffset    = pool.reserve((size_t)1            * RowStride);
    newGradLaplOffset = pool.reserve((size_t)4            * RowStride);
    AWorkOffset       = pool.reserve((size_t)2 * NumPtcls * RowStride);
    AinvWorkOffset    = pool.reserve((size_t)2 * NumPtcls * RowStride);
    Phi->reserve(pool);
  }

  void recompute (MCWalkerConfiguration &W, bool firstTime);

  void addLog (MCWalkerConfiguration &W, vector<RealType> &logPsi);

  void addGradient(MCWalkerConfiguration &W, int iat,
                   vector<GradType> &grad);

  void calcGradient(MCWalkerConfiguration &W, int iat,
                    vector<GradType> &grad);

  void ratio (MCWalkerConfiguration &W, int iat,
              vector<ValueType> &psi_ratios);


  void ratio (MCWalkerConfiguration &W, int iat,
              vector<ValueType> &psi_ratios,	vector<GradType>  &grad);

  void ratio (MCWalkerConfiguration &W, int iat,
              vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
              vector<ValueType> &lapl);
  void calcRatio (MCWalkerConfiguration &W, int iat,
                  vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
                  vector<ValueType> &lapl);
  void addRatio (MCWalkerConfiguration &W, int iat,
                 vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
                 vector<ValueType> &lapl);

  void ratio (vector<Walker_t*> &walkers, vector<int> &iatList,
              vector<PosType> &rNew, vector<ValueType> &psi_ratios,
              vector<GradType>  &grad, vector<ValueType> &lapl);

  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);

  void NLratios (MCWalkerConfiguration &W,  vector<NLjob> &jobList,
                 vector<PosType> &quadPoints, vector<ValueType> &psi_ratios);

  void NLratios_CPU (MCWalkerConfiguration &W,  vector<NLjob> &jobList,
                     vector<PosType> &quadPoints, vector<ValueType> &psi_ratios);
};
}
#endif // QMCPLUSPLUS_DIRAC_DETERMINANT_CUDA_H
