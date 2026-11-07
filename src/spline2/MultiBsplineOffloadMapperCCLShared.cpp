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


#include "MultiBsplineOffloadMapperCCLShared.hpp"
#include "OMPTarget/OMPrequires.hpp"
#include "OMPTarget/OMPTargetUsage.hpp"
#include <omp.h>
#include <string>
#include <stdexcept>

namespace qmcplusplus
{
extern MemoryUsageAccount omptarget_mem_usage;

template<typename T>
MultiBsplineOffloadMapperCCLShared<T>::MultiBsplineOffloadMapperCCLShared(
    const MultiBsplineMPIShared<T>& host_bsplines)
    : MultiBsplineOffloadMapperBase<T>(host_bsplines), host_bsplines_shared_(host_bsplines)
{
  mapToDevice();
}

template<typename T>
void MultiBsplineOffloadMapperCCLShared<T>::checkAndCoordinateError(bool local_failed,
                                                                    const std::string& local_err_msg)
{
  int any_failed = local_failed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &any_failed, 1, MPI_INT, MPI_LOR, host_bsplines_shared_.getComm().getMPI());
  if (any_failed)
  {
    cleanup();
    throw std::runtime_error(local_failed ? local_err_msg
                                          : "MultiBsplineOffloadMapperCCLShared setup failed on a peer rank.");
  }
}

template<typename T>
void MultiBsplineOffloadMapperCCLShared<T>::mapToDevice()
{
  const auto& comm                 = host_bsplines_shared_.getComm();
  const int comm_rank              = comm.rank();
  const int comm_size              = comm.size();
  const unsigned distributed_ranks = host_bsplines_shared_.getDistributedRanks();
  const bool contiguous_shared     = host_bsplines_shared_.getContiguousShared();
  const unsigned shared_ranks      = comm_size / distributed_ranks;

  std::string err_msg;
  cclResult_t res = cclSetDevice(omp_get_default_device());
  if (res != cclSuccess)
    err_msg = std::string("cclSetDevice failed: ") + cclGetErrorString(res);

  // 1. Create CCL communicator using the MPI communicator stored in MultiBsplineMPIShared
  cclUniqueId id{};
  if (err_msg.empty() && comm_rank == 0)
  {
    res = cclGetUniqueId(&id);
    if (res != cclSuccess)
      err_msg = std::string("cclGetUniqueId failed: ") + cclGetErrorString(res);
  }
  checkAndCoordinateError(!err_msg.empty(), err_msg);

  MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, comm.getMPI());
  res = cclCommInitRank(&ccl_comm_, comm_size, id, comm_rank);
  if (res != cclSuccess)
  {
    ccl_comm_ = nullptr;
    err_msg   = std::string("cclCommInitRank failed: ") + cclGetErrorString(res);
  }
  checkAndCoordinateError(!err_msg.empty(), err_msg);

  // 2. Identify the allocating ranks and their corresponding blocks
  if (contiguous_shared)
  {
    is_allocating_rank_ = (comm_rank % shared_ranks == 0) && (comm_rank / shared_ranks < distributed_ranks);
    allocated_block_    = comm_rank / shared_ranks;
  }
  else
  {
    is_allocating_rank_ = (comm_rank < distributed_ranks);
    allocated_block_    = comm_rank;
  }

  // 3. Allocate device memory using cclMemAlloc (with 1-element padding if coefs_size is 0)
  const size_t coefs_size  = is_allocating_rank_ ? host_bsplines_shared_.getBlock(allocated_block_).coefs_size : 0;
  const size_t alloc_bytes = (coefs_size > 0 ? coefs_size : 1) * sizeof(T);
  res                      = cclMemAlloc(&local_dev_ptr_, alloc_bytes);
  if (res != cclSuccess)
  {
    local_dev_ptr_ = nullptr;
    err_msg        = std::string("cclMemAlloc failed: ") + cclGetErrorString(res);
  }
  else if (is_allocating_rank_)
    omptarget_mem_usage.creditUsage(alloc_bytes);
  checkAndCoordinateError(!err_msg.empty(), err_msg);

  // 4. Register communication window across ranks
  res = cclCommWindowRegister(ccl_comm_, local_dev_ptr_, alloc_bytes, &win_, CCL_WIN_DEFAULT);
  if (res != cclSuccess)
  {
    win_    = nullptr;
    err_msg = std::string("cclCommWindowRegister failed: ") + cclGetErrorString(res);
  }
  checkAndCoordinateError(!err_msg.empty(), err_msg);

  // 5. Query device pointers for peer access for all distributed blocks
  block_coefs_dev_.resize(distributed_ranks, nullptr);
  for (int i = 0; i < distributed_ranks; i++)
  {
    const int target_rank = contiguous_shared ? i * shared_ranks : i;
    void* peer_ptr        = nullptr;
    res                   = cclGetPeerDevicePointer(win_, 0, target_rank, &peer_ptr);
    if (res != cclSuccess)
    {
      err_msg = std::string("cclGetPeerDevicePointer failed: ") + cclGetErrorString(res);
      break;
    }
    if (!peer_ptr)
    {
      err_msg = "cclGetPeerDevicePointer returned nullptr between rank " + std::to_string(comm_rank) +
          " and target rank " + std::to_string(target_rank) +
          ". Direct GPU peer access may not be supported across these ranks.";
      break;
    }
    block_coefs_dev_[i] = static_cast<const T*>(peer_ptr);
  }
  checkAndCoordinateError(!err_msg.empty(), err_msg);

  // 6. Map spline metadata to device
  for (int ib = 0; ib < host_bsplines_shared_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_shared_.getBlock(ib);
    if (spline_m->coefs_size > 0)
    {
      PRAGMA_OFFLOAD("omp target enter data map(to: spline_m[:1])")
      omptarget_mem_usage.creditUsage(sizeof(decltype(*spline_m)));
    }
  }
  metadata_mapped_ = true;
  comm.barrier();
}

template<typename T>
void MultiBsplineOffloadMapperCCLShared<T>::updateToDevice()
{
  const auto& comm = host_bsplines_shared_.getComm();
  comm.barrier();
  int status = 0;
  if (is_allocating_rank_)
  {
    const auto& spline_m = host_bsplines_shared_.getBlock(allocated_block_);
    if (spline_m.coefs_size > 0)
      status = omp_target_memcpy(local_dev_ptr_, spline_m.coefs, spline_m.coefs_size * sizeof(T), 0, 0,
                                 omp_get_default_device(), omp_get_initial_device());
  }
  MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_LOR, comm.getMPI());
  if (status != 0)
    throw std::runtime_error("omp_target_memcpy failed in MultiBsplineOffloadMapperCCLShared::updateToDevice!");
}

template<typename T>
void MultiBsplineOffloadMapperCCLShared<T>::cleanup() noexcept
{
  if (metadata_mapped_)
  {
    for (int ib = 0; ib < host_bsplines_shared_.getNumBlocks(); ib++)
    {
      auto* spline_m = &host_bsplines_shared_.getBlock(ib);
      if (spline_m->coefs_size > 0)
      {
        PRAGMA_OFFLOAD("omp target exit data map(delete: spline_m[:1])")
        omptarget_mem_usage.debitUsage(sizeof(decltype(*spline_m)));
      }
    }
    metadata_mapped_ = false;
  }

  if (win_)
  {
    cclCommWindowDeregister(ccl_comm_, win_);
    win_ = nullptr;
  }

  if (local_dev_ptr_)
  {
    cclMemFree(local_dev_ptr_);
    local_dev_ptr_ = nullptr;
    if (is_allocating_rank_)
    {
      const size_t coefs_size  = host_bsplines_shared_.getBlock(allocated_block_).coefs_size;
      const size_t alloc_bytes = (coefs_size > 0 ? coefs_size : 1) * sizeof(T);
      omptarget_mem_usage.debitUsage(alloc_bytes);
    }
  }

  if (ccl_comm_)
  {
    cclCommDestroy(ccl_comm_);
    ccl_comm_ = nullptr;
  }
}

template<typename T>
MultiBsplineOffloadMapperCCLShared<T>::~MultiBsplineOffloadMapperCCLShared()
{
  cleanup();
}

template class MultiBsplineOffloadMapperCCLShared<float>;
template class MultiBsplineOffloadMapperCCLShared<double>;
} // namespace qmcplusplus
