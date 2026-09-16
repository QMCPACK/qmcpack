//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MultiBsplineOffloadMapper.hpp"
#include "OMPTarget/OMPTargetMath.hpp"
#include "OMPTarget/OMPTargetUsage.hpp"

namespace qmcplusplus
{
extern MemoryUsageAccount omptarget_mem_usage;

template<typename T>
MultiBsplineOffloadMapper<T>::MultiBsplineOffloadMapper(const HostBspline& host_bsplines)
    : MultiBsplineOffloadMapperBase<T>(host_bsplines)
{
  mapToDevice();
}

template<typename T>
void MultiBsplineOffloadMapper<T>::mapToDevice()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target enter data map(to: spline_m[:1]) map(alloc: coefs[:spline_m->coefs_size])")
    omptarget_mem_usage.creditUsage(sizeof(decltype(*spline_m)) + spline_m->coefs_size * sizeof(T));
  }
}

template<typename T>
MultiBsplineOffloadMapper<T>::~MultiBsplineOffloadMapper()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target exit data map(delete: spline_m[:1]) map(delete: coefs[:spline_m->coefs_size])")
    omptarget_mem_usage.debitUsage(sizeof(decltype(*spline_m)) + spline_m->coefs_size * sizeof(T));
  }
}

template class MultiBsplineOffloadMapper<float>;
template class MultiBsplineOffloadMapper<double>;
} // namespace qmcplusplus
