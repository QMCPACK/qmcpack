//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   partially refactored from DiracDeterminantBase.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_EVAL_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_EVAL_H

#include "Configuration.h"
#include "QMCWaveFunctions/BsplineFactory/temp_batch_type.h"
#include "QMCWaveFunctions/SPOSetSingle.h"
#include "QMCWaveFunctions/SPOSetBatched.h"
#include "QMCWaveFunctions/Fermion/determinant_update.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/NLjob.h"
#include "QMCWaveFunctions/WaveFunctionComponentTypeAliases.h"
namespace qmcplusplus
{

class DiracDeterminantEvalDefault
{
public:
  using QMCT = QMCTraits;
  using SSTA = SPOSetTypeAliases;
  using WFCA = WaveFunctionComponentTypeAliases;
  using Walker_t =  WFCA::Walker_t;

  void abortNoSpecialize()
  {
    APP_ABORT("DiracDeterminantEvalDefault methods should not be reached");
  }
  
  virtual void addLog (MCWalkerConfiguration &W, std::vector<QMCT::RealType> &logPsi)
  { abortNoSpecialize(); }
  virtual void addGradient(MCWalkerConfiguration &W, int iat,
			   std::vector<QMCT::GradType> &grad)
  { abortNoSpecialize(); }
  virtual void calcGradient(MCWalkerConfiguration &W, int iat,
			    std::vector<QMCT::GradType> &grad)
  { abortNoSpecialize(); }
  virtual void ratio (MCWalkerConfiguration &W, int iat,
		      std::vector<QMCT::ValueType> &psi_ratios)
  { abortNoSpecialize(); }

  virtual void ratio (MCWalkerConfiguration &W, int iat,
		      std::vector<QMCT::ValueType> &psi_ratios,
		      std::vector<QMCT::GradType>  &grad)
  { abortNoSpecialize(); }
  virtual void ratio (MCWalkerConfiguration &W, int iat,
		      std::vector<QMCT::ValueType> &psi_ratios,
		      std::vector<QMCT::GradType>  &grad,
		      std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }
  virtual void ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
		      std::vector<QMCT::PosType> &rNew, std::vector<QMCT::ValueType> &psi_ratios,
		      std::vector<QMCT::GradType>  &grad, std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }
  virtual  void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
			  std::vector<QMCT::PosType> &quadPoints, std::vector<QMCT::ValueType> &psi_ratios)
  { abortNoSpecialize(); }
  virtual void calcRatio (MCWalkerConfiguration &W, int iat,
			  std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
                  std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }
  virtual void addRatio (MCWalkerConfiguration &W, int iat,
                 std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
                 std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }

  virtual void gradLapl (MCWalkerConfiguration &W, SSTA::GradMatrix_t &grads,
			 SSTA::ValueMatrix_t &lapl)
  { abortNoSpecialize(); }

  virtual QMCT::RealType updateBuffer(ParticleSet& P, WFCA::WFBufferType& buf, bool fromscratch=false)
  { abortNoSpecialize();
    return QMCT::RealType();
  }
  
  virtual void recompute(ParticleSet& P)
  {
    APP_ABORT("Need specialization of DiracDetermiantEval::recompute(ParticleSet& P).\n");
  }
  
  virtual void recompute(MCWalkerConfiguration &W, bool firstTime)
  {
    std::cerr << "Need specialization of DiracDetermiantEval::recompute.\n";
    abort();
  }

  virtual void update (std::vector<Walker_t*> &walkers, int iat)
  {
    APP_ABORT("Need specialization of DiracDetermiantEval::update.\n");
  }

  virtual void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList)
  {
    APP_ABORT("Need specialization of DiracDetermiantEval::update.\n");
  }


  
};

template<Batching batching>
class DiracDeterminantEval : public DiracDeterminantEvalDefault
{

};

template<>
class DiracDeterminantEval<Batching::SINGLE> : public DiracDeterminantEvalDefault
{
};

template<>
class DiracDeterminantEval<Batching::BATCHED> : public DiracDeterminantEvalDefault
{
public:
  DiracDeterminantEval();
  
  using CudaValueType = QMCT::CudaValueType;
  using CudaRealType = QMCT::CudaRealType;
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

  void resizeLists(int numWalkers, int num_orbitals)
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
    // gradLapl_d.resize   (numWalkers*num_orbitals*4);
    // gradLapl_host.resize(numWalkers*num_orbitals*4);
    infoArray_d.resize(numWalkers*2);
    infoArray_host.resize(numWalkers*2);
    PivotArray_d.resize(numWalkers*num_orbitals);
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

  void reserve (PointerPool<gpu::device_vector<QMCT::CudaValueType> > &pool, SPOSetBatched& Phi,
  		int num_particles, int num_orbitals);

  // void update (std::vector<Walker_t*> &walkers, int iat);
  // void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList);
  // void addLog (MCWalkerConfiguration &W, std::vector<QMCT::RealType> &logPsi);

  // void addGradient(MCWalkerConfiguration &W, int iat,
  //                  std::vector<QMCT::GradType> &grad);

  // void calcGradient(MCWalkerConfiguration &W, int iat,
  //                   std::vector<QMCT::GradType> &grad);

  // void ratio (MCWalkerConfiguration &W, int iat,
  //             std::vector<QMCT::ValueType> &psi_ratios);

  // void ratio (MCWalkerConfiguration &W, int iat,
  //             std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad);

  // void ratio (MCWalkerConfiguration &W, int iat,
  //             std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
  //             std::vector<QMCT::ValueType> &lapl);
  // void calcRatio (MCWalkerConfiguration &W, int iat,
  //                 std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
  //                 std::vector<QMCT::ValueType> &lapl);
  // void addRatio (MCWalkerConfiguration &W, int iat,
  //                std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
  //                std::vector<QMCT::ValueType> &lapl);
  // void gradLapl (MCWalkerConfiguration &W, SSTA::GradMatrix_t &grads,
  //                SSTA::ValueMatrix_t &lapl);


};

  
}

#endif
