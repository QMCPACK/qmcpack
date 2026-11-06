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
// -*- C++ -*-
/** @file MultiBsplineOffloadMapperCCLShared.hpp
 * @brief Defines the MultiBsplineOffloadMapperCCLShared class for GPU offload mapping
 *        using NCCL/CCL shared memory windows.
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINEOFFLOADMAPPER_CCLSHARED_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEOFFLOADMAPPER_CCLSHARED_HPP

#include "MultiBsplineOffloadMapperBase.hpp"
#include "MultiBsplineMPIShared.hpp"
#include "ccl_wrapper.hpp"

namespace qmcplusplus
{
/** @ingroup spline2
 * @brief Offload mapper that manages 3D multi-B-spline coefficients distributed across GPUs
 *        using CCL (NCCL/RCCL/oneCCL) memory windows.
 *
 * Analogous to MultiBsplineMPIShared on the host, MultiBsplineOffloadMapperCCLShared
 * allocates device memory using cclMemAlloc, registers a shared window across ranks via
 * cclCommWindowRegister, and retrieves direct peer device pointers via cclGetPeerDevicePointer.
 *
 * @tparam T The precision of spline coefficients (float or double).
 */
template<typename T>
class MultiBsplineOffloadMapperCCLShared : public MultiBsplineOffloadMapperBase<T>
{
private:
  using Base                 = MultiBsplineOffloadMapperBase<T>;
  using HostBsplineMPIShared = MultiBsplineMPIShared<T>;
  using Base::block_coefs_dev_;
  using Base::host_bsplines_;

  /// Reference to the host distributed and shared splines.
  const HostBsplineMPIShared& host_bsplines_shared_;

  /// CCL communicator created from the host MPI communicator.
  cclComm_t ccl_comm_{nullptr};
  /// CCL memory window handle.
  cclWindow_t win_{nullptr};
  /// Local device memory buffer allocated via cclMemAlloc.
  void* local_dev_ptr_{nullptr};
  /// Whether this rank allocates device memory for a block.
  bool is_allocating_rank_{false};
  /// Block index allocated by this rank (valid if is_allocating_rank_ is true).
  size_t allocated_block_{0};
  /// Whether spline block metadata has been mapped to device.
  bool metadata_mapped_{false};

  /// Map host spline metadata and coefficients to device.
  void mapToDevice();
  /// Release all allocated CCL and device resources.
  void cleanup() noexcept;
  /// Coordinate error status across the communicator; clean up and throw on all ranks if any rank failed.
  void checkAndCoordinateError(bool local_failed, const std::string& local_err_msg);

public:
  MultiBsplineOffloadMapperCCLShared(const MultiBsplineMPIShared<T>& host_bsplines);

  ~MultiBsplineOffloadMapperCCLShared() override;

  /// Update device coefficients for the locally owned block and synchronize.
  void updateToDevice() override;
};

extern template class MultiBsplineOffloadMapperCCLShared<float>;
extern template class MultiBsplineOffloadMapperCCLShared<double>;
} // namespace qmcplusplus

#endif
