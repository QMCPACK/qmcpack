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
#include "config.h"
#include "Synchro.hpp"

namespace qmcplusplus
{

#ifdef ENABLE_CUDA
CudaStreamSynchro::CudaStreamSynchro() : stream_(0) {}
CudaStreamSynchro::CudaStreamSynchro(cudaStream_t stream) : stream_(stream) {}
void CudaStreamSynchro::sync() { cudaErrorCheck(cudaStreamSynchronize(stream_), "cudaStreamSynchronize failed!"); }
const cudaStream_t& CudaStreamSynchro::get_stream() const { return stream_; };
#endif

/// Little concerned about whether the opimizer may ruin this.
void OpenMPSynchro::sync()
{
  PRAGMA_OFFLOAD("omp taskwait")
}

}

