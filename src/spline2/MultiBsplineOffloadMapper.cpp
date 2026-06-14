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
#include "MultiBsplineEval_OMPoffload.hpp"
#include "OMPTarget/OMPTargetMath.hpp"

namespace qmcplusplus
{
template<typename T>
MultiBsplineOffloadMapper<T>::MultiBsplineOffloadMapper(const HostBspline& host_bsplines)
    : host_bsplines_(host_bsplines)
{
  block_coefs_.reserve(host_bsplines_.getNumBlocks());
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* coefs = host_bsplines_.getBlock(ib).coefs;
    block_coefs_.push_back(coefs);
  }
}

template<typename T>
void MultiBsplineOffloadMapper<T>::mapToDevice()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target enter data map(to: spline_m[:1]) map(alloc: coefs[:spline_m->coefs_size])")
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
  }
}

template<typename T>
void MultiBsplineOffloadMapper<T>::updateToDevice()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target update to(coefs[:spline_m->coefs_size])")
  }
}

template<typename T>
void MultiBsplineOffloadMapper<T>::mw_evaluate_v(int num_pos, T* pos_arr, T* spline_v, size_t walker_stride)
{
  const auto block_offsets = host_bsplines_.getBlockOffsets();
  for (size_t ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    const auto& host_block = host_bsplines_.getBlock(ib);
    if (host_block.num_splines == 0)
      continue;

    const size_t num_splines      = host_block.num_splines;
    const size_t ChunkSizePerTeam = 512;
    const int NumTeams            = (num_splines + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    // Ye: need to extract sizes and pointers before entering target region
    const auto* spline_ptr   = &host_block;
    const auto* spline_coefs = block_coefs_[ib];
    const auto block_offset  = block_offsets[ib];

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams * num_pos)")
    for (int iw = 0; iw < num_pos; iw++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const size_t first = ChunkSizePerTeam * team_id;
        const size_t last  = omptarget::min(first + ChunkSizePerTeam, num_splines);

        auto* spline_v_iw = spline_v + walker_stride * iw;
        int ix, iy, iz;
        T a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr, pos_arr[iw * 3], pos_arr[iw * 3 + 1], pos_arr[iw * 3 + 2], ix,
                                              iy, iz, a, b, c);

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_v_impl_v2(spline_ptr, spline_coefs, ix, iy, iz, first + index, a, b, c,
                                             spline_v_iw + block_offset + first + index);
      }
  }
}

template<typename T>
void MultiBsplineOffloadMapper<T>::mw_evaluate_vgh(int num_pos,
                                                   T* pos_arr,
                                                   T* spline_vgh,
                                                   size_t walker_stride,
                                                   size_t field_stride)
{
  const auto block_offsets = host_bsplines_.getBlockOffsets();
  for (size_t ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    const auto& host_block = host_bsplines_.getBlock(ib);
    if (host_block.num_splines == 0)
      continue;

    const size_t num_splines      = host_block.num_splines;
    const size_t ChunkSizePerTeam = 512;
    const int NumTeams            = (num_splines + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    // Ye: need to extract sizes and pointers before entering target region
    const auto* spline_ptr   = &host_block;
    const auto* spline_coefs = block_coefs_[ib];
    const auto block_offset  = block_offsets[ib];

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams * num_pos)")
    for (int iw = 0; iw < num_pos; iw++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const size_t first = ChunkSizePerTeam * team_id;
        const size_t last  = omptarget::min(first + ChunkSizePerTeam, num_splines);

        auto* spline_vgh_iw = spline_vgh + walker_stride * iw;
        int ix, iy, iz;
        T a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
        spline2::computeLocationAndFractional(spline_ptr, pos_arr[iw * 3], pos_arr[iw * 3 + 1], pos_arr[iw * 3 + 2], ix,
                                              iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c);

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_vgh_impl_v2(spline_ptr, spline_coefs, ix, iy, iz, first + index, a, b, c, da, db, dc,
                                               d2a, d2b, d2c, spline_vgh_iw + block_offset + first + index,
                                               field_stride);
      }
  }
}

template class MultiBsplineOffloadMapper<float>;
template class MultiBsplineOffloadMapper<double>;
} // namespace qmcplusplus
