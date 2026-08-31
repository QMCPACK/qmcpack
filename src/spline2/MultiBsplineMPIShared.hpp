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
/**@file MultiBspline.hpp
 *
 * define classes MultiBsplineMPIShared
 * The evaluation functions are defined in MultiBsplineEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_MPISHARED_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_MPISHARED_HPP

#include <memory>
#include "MultiBsplineBase.hpp"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "Message/Communicate.h"
#include "Utilities/FairDivide.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
/** container class to hold a 3D multi spline pointer and BsplineAllocator
 * @tparam T the precision of splines
 */
template<typename T>
class MultiBsplineMPIShared : public MultiBsplineBase<T>
{
private:
  using Base  = MultiBsplineBase<T>;
  using Alloc = aligned_allocator<T>;
  /// communicator covering all the ranks included distributed and shared. ranks with distributed splines are contiguous.
  const std::unique_ptr<Communicate> comm_;

  MPI_Win win;
  /// number of ranks among which splines are distributed.
  const unsigned distributed_ranks_;

  using Base::offsets_;

public:
  template<typename BCT>
  MultiBsplineMPIShared(const Ugrid grid[3],
                        const BCT& bc,
                        size_t num_splines,
                        std::unique_ptr<Communicate>&& comm_distributed_and_shared,
                        unsigned distributed_ranks)
      : Base(FairDivideAligned<std::vector<size_t>>(num_splines, getAlignment<T>(), distributed_ranks)),
        comm_(std::move(comm_distributed_and_shared)),
        distributed_ranks_(distributed_ranks)
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

    auto& spline_owned = *Base::spline_blocks[comm_rank % distributed_ranks_];
    // Only the members of the first group of distributed ranks allocate memory.
    // Other ranks are considered references and hence size = 0 is valid.
    const MPI_Aint allocation_size =
        comm_rank < distributed_ranks_ ? spline_owned.coefs_size * sizeof(T) + Alloc::alignment : 0;
    void* coefs = nullptr;
    auto err = MPI_Win_allocate_shared(allocation_size, sizeof(T), info, comm.getMPI(), &coefs, &win);
    MPI_Info_free(&info);
    if (err != MPI_SUCCESS)
      throw UniformCommunicateError("MultiBsplineMPIShared::MultiBsplineMPIShared MPI_Win_allocate_shared failed!");

    spline_owned.coefs = (T*)coefs;

    for (int i = 0; i < distributed_ranks_; i++)
    {
      MPI_Aint size;
      int disp_unit;
      void* base_ptr;
      err = MPI_Win_shared_query(win, i, &size, &disp_unit, &base_ptr);
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

  ~MultiBsplineMPIShared() override;
};

extern template class MultiBsplineMPIShared<float>;
extern template class MultiBsplineMPIShared<double>;
} // namespace qmcplusplus

#endif
