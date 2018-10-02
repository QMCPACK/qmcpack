//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSET_BATCHED_H
#define QMCPLUSPLUS_SPOSET_BATCHED_H

#include "QMCWaveFunctions/Batching.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/SPOSet.h"
//! SPOSet evaluation interface depends on walker batching strategy
/*!
  SPOSet inherits this so the correct evaluation function signatures and translations
  for the walker batching strategy are specified.  The actual evaluations
  are implement by SplineAdoptor descendents.
*/
namespace qmcplusplus
{

template<>
class SPOSet<Batching::BATCHED>: public SPOSet<Batching::SINGLE>
{
public:
  using SPOSetPtr = SPOSet<Batching::BATCHED>*;

  SPOSet() : WhatAmI("SPOSet<BATCHED>") {}
  
  
  virtual void resetParameters(const opt_variables_type& active)
  { }

  virtual void resetTargetParticleSet(ParticleSet& e)
  { }

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
