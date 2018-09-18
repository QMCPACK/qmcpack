//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   initially refactored from SlaterDet.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SLATERDETERMINANT_BATCHED_H 
#define QMCPLUSPLUS_SLATERDETERMINANT_BATCHED_H

#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"

namespace qmcplusplus
{

template<>
class SlaterDet<Batching::BATCHED> : public SlaterDet<Batching::SINGLE>
{
  using Determinant_t = DiracDeterminant<Batching::BATCHED>;
  std::vector<Determinant_t*> Dets;
  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
  GPU_XRAY_TRACE void  recompute(MCWalkerConfiguration &W, bool firstTime)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->recompute(W, firstTime);
  }

  GPU_XRAY_TRACE void  reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->reserve(pool);
  }

  GPU_XRAY_TRACE void  addLog (MCWalkerConfiguration &W, std::vector<RealType> &logPsi)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->addLog(W, logPsi);
  }

  GPU_XRAY_TRACE void ratio (MCWalkerConfiguration &W, int iat
         , std::vector<ValueType> &psi_ratios,std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    Dets[getDetID(iat)]->ratio(W, iat, psi_ratios, grad, lapl);
  }

  GPU_XRAY_TRACE void calcRatio (MCWalkerConfiguration &W, int iat
             , std::vector<ValueType> &psi_ratios,std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    Dets[getDetID(iat)]->calcRatio(W, iat, psi_ratios, grad, lapl);
  }

  GPU_XRAY_TRACE void addRatio (MCWalkerConfiguration &W, int iat
            , std::vector<ValueType> &psi_ratios,std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    Dets[getDetID(iat)]->addRatio(W, iat, psi_ratios, grad, lapl);
  }

  GPU_XRAY_TRACE void  ratio (std::vector<Walker_t*> &walkers,    std::vector<int> &iatList,
              std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
              std::vector<GradType>  &grad, std::vector<ValueType> &lapl);

  GPU_XRAY_TRACE void  calcGradient(MCWalkerConfiguration &W, int iat, std::vector<GradType> &grad)
  {
    Dets[getDetID(iat)]->calcGradient(W, iat, grad);
  }

  GPU_XRAY_TRACE void  addGradient(MCWalkerConfiguration &W, int iat, std::vector<GradType> &grad)
  {
    Dets[getDetID(iat)]->addGradient(W, iat, grad);
  }

  GPU_XRAY_TRACE void  update (std::vector<Walker_t*> &walkers, int iat)
  {
    Dets[getDetID(iat)]->update(walkers, iat);
  }

  GPU_XRAY_TRACE void  update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList);

  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads, ValueMatrix_t &lapl)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->gradLapl(W, grads, lapl);
  }

  GPU_XRAY_TRACE void  NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList
                 , std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->NLratios(W, jobList, quadPoints, psi_ratios);
  }


};
  
}
#endif
