//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   refactored from SPOSet.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSET_EVALUATION_BATCHED_H
#define QMCPLUSPLUS_SPOSET_EVALUATION_BATCHED_H

#include "Configuration.h"
#include "Batching.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/SPOSetEvaluation.h"
namespace qmcplusplus
{

//! SPOSetEvaluation batched walker specialization evaluation
/*!
  SPOSetXXXX inherits this so the correct evaluation function signatures and translations
  for the walker batching strategy are specified.  The actual evaluations
  are implement by SplineAdoptor descendents.
*/
template<>
class SPOSetEvaluation<Batching::BATCHED> : public SPOSetEvaluationDefault
{
public:
  using SSTA = SPOSetTypeAliases;
  using QMCT = QMCTraits;
  virtual void initGPU() {  }

  //////////////////////////////////////////
  // Walker-parallel vectorized functions //
  //////////////////////////////////////////
  virtual void
  reserve (PointerPool<gpu::device_vector<QMCT::CudaValueType> > &pool) { }

  virtual void
  evaluate (std::vector<SSTA::Walker_t*> &walkers, int iat, gpu::device_vector<QMCT::CudaValueType*> &phi);

  virtual void evaluate (std::vector<SSTA::Walker_t*> &walkers, std::vector<QMCT::PosType> &new_pos
                         , gpu::device_vector<QMCT::CudaValueType*> &phi);

  virtual void
  evaluate (std::vector<SSTA::Walker_t*> &walkers,
            std::vector<QMCT::PosType> &new_pos,
            gpu::device_vector<QMCT::CudaValueType*> &phi,
            gpu::device_vector<QMCT::CudaValueType*> &grad_lapl_list,
            int row_stride);

  virtual void
  evaluate (std::vector<QMCT::PosType> &pos, gpu::device_vector<QMCT::CudaRealType*> &phi);
  virtual void
  evaluate (std::vector<QMCT::PosType> &pos, gpu::device_vector<QMCT::CudaComplexType*> &phi);

};




}

#endif
