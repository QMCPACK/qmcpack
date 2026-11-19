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
/** @file MultiBsplineMPIShared.hpp
 * @brief Defines the MultiBsplineMPIShared class for distributed and shared MPI 3D multi-B-splines.
 *
 * Evaluation functions are defined in MultiBsplineEval.hpp.
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_MPISHARED_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_MPISHARED_HPP

#include <memory>
#include <cassert>
#include "MultiBsplineBase.hpp"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "Message/Communicate.h"
#include "Utilities/FairDivide.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
/** @ingroup spline2
 * @brief Container class holding 3D multi-B-splines distributed and shared via MPI-3 shared memory.
 *
 * MultiBsplineMPIShared partitions a set of 3D multi-B-splines into multiple blocks
 * distributed across MPI ranks within a shared-memory domain (node).
 *
 * Memory allocation and layout:
 * - The splines are partitioned into @c distributed_ranks blocks using FairDivideAligned.
 * - An MPI shared memory window (@c MPI_Win_allocate_shared) is allocated over @c comm_.
 * - Only one designated rank per shared group allocates storage for its assigned block;
 *   other ranks allocate 0 bytes in the window.
 * - Using @c MPI_Win_shared_query, all ranks query and attach to all distributed blocks,
 *   enabling direct zero-copy read access to every spline block during evaluation.
 *
 * Rank layout (@c contiguous_shared):
 * - If @c contiguous_shared is true, ranks sharing each block are contiguous:
 *   shared group @c b spans ranks <tt>[b * shared_ranks, (b + 1) * shared_ranks)</tt>, and
 *   rank <tt>b * shared_ranks</tt> allocates Block @c b.
 * - If @c contiguous_shared is false, ranks allocating distributed blocks are contiguous:
 *   ranks <tt>0, ..., distributed_ranks - 1</tt> allocate Blocks <tt>0, ..., distributed_ranks - 1</tt>,
 *   and ranks sharing each block are interleaved/strided.
 *
 * @tparam T The precision of the spline coefficients (float or double).
 */
template<typename T>
class MultiBsplineMPIShared : public MultiBsplineBase<T>
{
private:
  using Base  = MultiBsplineBase<T>;
  using Alloc = aligned_allocator<T>;
  /// Communicator covering all participating ranks (distributed and shared).
  const std::unique_ptr<Communicate> comm_;

  /// MPI shared-memory window for spline coefficient storage.
  MPI_Win win;
  /// Number of blocks across which splines are distributed.
  const unsigned distributed_ranks_;
  /// Whether ranks sharing memory are contiguous (true) or distributed allocating ranks are contiguous (false).
  const bool contiguous_shared_;

  using Base::offsets_;

public:
  /** @brief Constructor for distributed and shared MPI 3D multi-B-splines.
   * @tparam BCT Boundary condition type.
   * @param grid Uniform grid specification in 3 dimensions.
   * @param bc Boundary conditions.
   * @param num_splines Total number of splines.
   * @param comm_distributed_and_shared Communicator covering all participating ranks (distributed and shared).
   * @param distributed_ranks Number of blocks across which the spline coefficients are distributed.
   * @param contiguous_shared If true, ranks sharing each block are contiguous (allocating rank for block i is
   *                          i * shared_ranks). If false, ranks allocating distributed blocks are contiguous
   *                          (allocating rank for block i is i). Default is false.
   */
  template<typename BCT>
  MultiBsplineMPIShared(const Ugrid grid[3],
                        const BCT& bc,
                        size_t num_splines,
                        std::unique_ptr<Communicate> comm_distributed_and_shared,
                        unsigned distributed_ranks,
                        bool contiguous_shared = false)
      : Base(FairDivideAligned<std::vector<size_t>>(num_splines, getAlignment<T>(), distributed_ranks)),
        comm_(std::move(comm_distributed_and_shared)),
        distributed_ranks_(distributed_ranks),
        contiguous_shared_(contiguous_shared)
  {
    auto& comm           = *comm_;
    const auto comm_rank = comm.rank();

    Base::spline_blocks.resize(distributed_ranks_, nullptr);
    for (int i = 0; i < distributed_ranks_; i++)
    {
      const auto num_splines = offsets_[i + 1] - offsets_[i];
      auto* spline_m         = new typename Base::SplineType;
      Base::spline_blocks[i] = spline_m;
      Base::setMetaData(*spline_m, grid[0], grid[1], grid[2], Base::createBoundaryCondition(bc).data(), num_splines,
                        getAlignedSize<T, Alloc::alignment>(num_splines));
    }

    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "alloc_shared_noncontig", "true");
    MPI_Info_set(info, "mpi_minimum_memory_alignment", std::to_string(Alloc::alignment).c_str());

    assert((comm.size() % distributed_ranks_ == 0) &&
           "MultiBsplineMPIShared::MultiBsplineMPIShared comm.size() must be divisible by distributed_ranks_!");

    const unsigned shared_ranks = comm.size() / distributed_ranks_;

    size_t allocated_block = 0;
    bool is_allocating_rank = false;
    if (contiguous_shared_)
    {
      is_allocating_rank = (comm_rank % shared_ranks == 0) && (comm_rank / shared_ranks < distributed_ranks_);
      allocated_block    = comm_rank / shared_ranks;
    }
    else
    {
      is_allocating_rank = (comm_rank < distributed_ranks_);
      allocated_block    = comm_rank;
    }

    const MPI_Aint allocation_size =
        is_allocating_rank ? Base::spline_blocks[allocated_block]->coefs_size * sizeof(T) + Alloc::alignment : 0;
    void* coefs = nullptr;
    auto err    = MPI_Win_allocate_shared(allocation_size, sizeof(T), info, comm.getMPI(), &coefs, &win);
    MPI_Info_free(&info);
    if (err != MPI_SUCCESS)
      throw UniformCommunicateError("MultiBsplineMPIShared::MultiBsplineMPIShared MPI_Win_allocate_shared failed!");

    for (int i = 0; i < distributed_ranks_; i++)
    {
      MPI_Aint size;
      int disp_unit;
      void* base_ptr;
      const int target_rank = contiguous_shared_ ? i * shared_ranks : i;
      err                   = MPI_Win_shared_query(win, target_rank, &size, &disp_unit, &base_ptr);
      if (err != MPI_SUCCESS)
        throw UniformCommunicateError("MultiBsplineMPIShared::MultiBsplineMPIShared MPI_Win_shared_query failed!");
      // Overallocate above and adjust the alignment here. Required for OpenMPI 4.x
      auto addr_mod  = (uintptr_t)base_ptr % Alloc::alignment;
      auto& spline_m = *Base::spline_blocks[i];
      spline_m.coefs = (T*)((uintptr_t)base_ptr + Alloc::alignment - addr_mod);
      if (((uintptr_t)spline_m.coefs & (Alloc::alignment - 1)) != 0)
        throw UniformCommunicateError(
            "MultiBsplineMPIShared::MultiBsplineMPIShared spline_m.coefs address not aligned!");
    }
  }

  /// Frees the shared memory window and deletes local spline block metadata.
  ~MultiBsplineMPIShared() override;

  const Communicate& getComm() const { return *comm_; }
  unsigned getDistributedRanks() const { return distributed_ranks_; }
  bool getContiguousShared() const { return contiguous_shared_; }
};

extern template class MultiBsplineMPIShared<float>;
extern template class MultiBsplineMPIShared<double>;
} // namespace qmcplusplus

#endif
