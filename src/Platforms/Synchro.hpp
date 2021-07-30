//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 */
#ifndef QMCPLUSPLUS_SYNCHRO_HPP
#define QMCPLUSPLUS_SYNCHRO_HPP

#ifdef ENABLE_CUDA
#include <cuda_runtime_api.h>
#include "CUDA/cudaError.h"
#endif

namespace qmcplusplus
{


/** Abstracts the synchronization context of memory transfer.  This is implicit in OMPTarget
 *  but explicit in CUDA, HIP and likely other ways of writing device code.
 *
 *  This prevents contamination of DualAllocator and API's with device code specific objects.
 */
class Synchro
{
public:
  virtual void sync() {};
};

#ifdef ENABLE_CUDA
class CudaStreamSynchro : public Synchro
{
public:
  CudaStreamSynchro();
  CudaStreamSynchro(cudaStream_t stream);
  void sync() override;
  const cudaStream_t& get_stream() const;
private:
  /// initialize to default stream
  cudaStream_t stream_;
};
#else
class CudaStreamSynchro : public Synchro
{};
#endif

class OpenMPSynchro : public Synchro
{
public:
  void sync() override;
};

} // namespace qmcplusplus

#endif
