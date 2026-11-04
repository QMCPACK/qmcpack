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
#include <nccl_device.h>

namespace qmcplusplus
{
extern MemoryUsageAccount omptarget_mem_usage;

#define NCCLCHECK(cmd)                                                                                \
  do                                                                                                  \
  {                                                                                                   \
    ncclResult_t res = (cmd);                                                                         \
    if (res != ncclSuccess)                                                                           \
      throw std::runtime_error(std::string("NCCL error at ") + __FILE__ + ":" +                      \
                               std::to_string(__LINE__) + " - " + ncclGetErrorString(res));           \
  } while (0)

template<typename T>
MultiBsplineOffloadMapperCCLShared<T>::MultiBsplineOffloadMapperCCLShared(
    const MultiBsplineMPIShared<T>& host_bsplines)
    : MultiBsplineOffloadMapperBase<T>(host_bsplines), host_bsplines_shared_(host_bsplines)
{
  mapToDevice();
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

  // 1. Create NCCL communicator using the MPI communicator stored in MultiBsplineMPIShared
  ncclUniqueId id;
  if (comm_rank == 0)
    NCCLCHECK(ncclGetUniqueId(&id));
  MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, comm.getMPI());
  NCCLCHECK(ncclCommInitRank(&nccl_comm_, comm_size, id, comm_rank));

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

  // 3. Allocate device memory using ncclMemAlloc (with 1-element padding if coefs_size is 0)
  const size_t coefs_size  = is_allocating_rank_ ? host_bsplines_shared_.getBlock(allocated_block_).coefs_size : 0;
  const size_t alloc_bytes = (coefs_size > 0 ? coefs_size : 1) * sizeof(T);
  NCCLCHECK(ncclMemAlloc(&local_dev_ptr_, alloc_bytes));

  if (is_allocating_rank_)
    omptarget_mem_usage.creditUsage(alloc_bytes);

  // 4. Register communication window across ranks
  NCCLCHECK(ncclCommWindowRegister(nccl_comm_, local_dev_ptr_, alloc_bytes, &win_, 0));

  // 5. Query device pointers for peer access for all distributed blocks
  block_coefs_dev_.resize(distributed_ranks);
  for (int i = 0; i < distributed_ranks; i++)
  {
    const int target_rank = contiguous_shared ? i * shared_ranks : i;
    void* peer_ptr        = nullptr;
    NCCLCHECK(ncclGetPeerDevicePointer(win_, 0, target_rank, &peer_ptr));
    block_coefs_dev_[i] = static_cast<const T*>(peer_ptr);
  }

  // 6. Map spline metadata to device
  for (int ib = 0; ib < host_bsplines_shared_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_shared_.getBlock(ib);
    PRAGMA_OFFLOAD("omp target enter data map(to: spline_m[:1])")
    omptarget_mem_usage.creditUsage(sizeof(decltype(*spline_m)));
  }
  comm.barrier();
}

template<typename T>
void MultiBsplineOffloadMapperCCLShared<T>::updateToDevice()
{
  if (is_allocating_rank_)
  {
    const auto& spline_m = host_bsplines_shared_.getBlock(allocated_block_);
    if (spline_m.coefs_size > 0)
    {
      const int status = omp_target_memcpy(local_dev_ptr_, spline_m.coefs, spline_m.coefs_size * sizeof(T), 0, 0,
                                           omp_get_default_device(), omp_get_initial_device());
      if (status != 0)
        throw std::runtime_error("omp_target_memcpy failed in MultiBsplineOffloadMapperCCLShared::updateToDevice!");
    }
  }
  host_bsplines_shared_.getComm().barrier();
}

template<typename T>
MultiBsplineOffloadMapperCCLShared<T>::~MultiBsplineOffloadMapperCCLShared()
{
  for (int ib = 0; ib < host_bsplines_shared_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_shared_.getBlock(ib);
    PRAGMA_OFFLOAD("omp target exit data map(delete: spline_m[:1])")
    omptarget_mem_usage.debitUsage(sizeof(decltype(*spline_m)));
  }

  if (win_)
  {
    ncclCommWindowDeregister(nccl_comm_, win_);
    win_ = nullptr;
  }

  if (local_dev_ptr_)
  {
    ncclMemFree(local_dev_ptr_);
    local_dev_ptr_ = nullptr;
  }

  if (is_allocating_rank_)
  {
    const size_t coefs_size  = host_bsplines_shared_.getBlock(allocated_block_).coefs_size;
    const size_t alloc_bytes = (coefs_size > 0 ? coefs_size : 1) * sizeof(T);
    omptarget_mem_usage.debitUsage(alloc_bytes);
  }

  if (nccl_comm_)
  {
    ncclCommDestroy(nccl_comm_);
    nccl_comm_ = nullptr;
  }
}

template class MultiBsplineOffloadMapperCCLShared<float>;
template class MultiBsplineOffloadMapperCCLShared<double>;
} // namespace qmcplusplus
