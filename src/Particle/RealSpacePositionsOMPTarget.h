//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file RealSpacePostionsOMPTarget.h
 */
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONS_OMPTARGET_H
#define QMCPLUSPLUS_REALSPACE_POSITIONS_OMPTARGET_H

#include "Particle/DynamicCoordinates.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "OMPTarget/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "ParticleSet.h"

namespace qmcplusplus
{
/** Introduced to handle virtual moves and ratio computations, e.g. for non-local PP evaluations.
   */
class RealSpacePositionsOMPTarget : public DynamicCoordinates
{
public:
  RealSpacePositionsOMPTarget() : DynamicCoordinates(DynamicCoordinateKind::DC_POS_OFFLOAD), is_host_position_changed_(false) {}
  RealSpacePositionsOMPTarget(const RealSpacePositionsOMPTarget& in)
      : DynamicCoordinates(DynamicCoordinateKind::DC_POS_OFFLOAD), RSoA(in.RSoA)
  {
    RSoA_hostview.attachReference(RSoA.size(), RSoA.capacity(), RSoA.data());
    updateH2D();
  }

  std::unique_ptr<DynamicCoordinates> makeClone() override
  {
    return std::make_unique<RealSpacePositionsOMPTarget>(*this);
  }

  void resize(size_t n) override
  {
    if (RSoA.size() != n)
    {
      RSoA.resize(n);
      RSoA_hostview.attachReference(RSoA.size(), RSoA.capacity(), RSoA.data());
    }
  }

  size_t size() override { return RSoA_hostview.size(); }

  void setAllParticlePos(const ParticlePos_t& R) override
  {
    resize(R.size());
    RSoA_hostview.copyIn(R);
    updateH2D();
  }

  void setOneParticlePos(const PosType& pos, size_t iat) override
  {
    RSoA_hostview(iat) = pos;
    is_host_position_changed_ = true;
    /* This was too slow due to overhead.
    RealType x     = pos[0];
    RealType y     = pos[1];
    RealType z     = pos[2];
    RealType* data = RSoA.data();
    size_t offset  = RSoA.capacity();

    PRAGMA_OFFLOAD("omp target map(to : x, y, z, iat)")
    {
      data[iat]              = x;
      data[iat + offset]     = y;
      data[iat + offset * 2] = z;
    }
    */
  }

  void mw_copyActiveOldParticlePos(const RefVector<DynamicCoordinates>& coords_list,
                                   size_t iat,
                                   const std::vector<PosType>& new_positions) override
  {
    const auto nw = coords_list.size();

    nw_new_old_pos.resize(nw);
    nw_active_pos_view.attachReference(nw_new_old_pos.size(), nw_new_old_pos.capacity(), nw_new_old_pos.data(0));
    nw_old_pos_view.attachReference(nw_new_old_pos.size(), nw_new_old_pos.capacity(), nw_new_old_pos.data(QMCTraits::DIM));

    for (int iw = 0; iw < nw; iw++)
    {
      auto& coords = dynamic_cast<const RealSpacePositionsOMPTarget&>(coords_list[iw].get());
      nw_active_pos_view(iw) = new_positions[iw];
      nw_old_pos_view(iw) = coords.RSoA_hostview[iat];
    }

    auto* mw_pos_ptr = nw_new_old_pos.data();
    PRAGMA_OFFLOAD("omp target update to(mw_pos_ptr[:QMCTraits::DIM * 2 * nw_new_old_pos.capacity()])")
  }

  void mw_acceptParticlePos(const RefVector<DynamicCoordinates>& coords_list,
                            size_t iat,
                            const std::vector<PosType>& new_positions,
                            const std::vector<bool>& isAccepted) override
  {
    const size_t nw = coords_list.size();
    nw_accept_index_ptrs.resize((sizeof(int) + sizeof(RealType*)) * nw);
    auto* RSoA_ptr_array = reinterpret_cast<RealType**>(nw_accept_index_ptrs.data());
    auto* id_array = reinterpret_cast<int*>(nw_accept_index_ptrs.data() + sizeof(RealType*) * coords_list.size());

    size_t num_accepted = 0;
    for (int iw = 0; iw < nw; iw++)
      if (isAccepted[iw])
      {
        auto& coords = dynamic_cast<RealSpacePositionsOMPTarget&>(coords_list[iw].get());
        RSoA_ptr_array[num_accepted] = coords.RSoA.device_data();
        id_array[num_accepted] = iw;
        // save new coordinates on host copy
        coords.RSoA_hostview(iat) = nw_active_pos_view[iw];
        num_accepted++;
      }

    //offload to GPU
    auto* restrict w_accept_buffer_ptr = nw_accept_index_ptrs.data();
    auto* restrict mw_pos_ptr = nw_new_old_pos.data();
    const size_t rsoa_stride = RSoA.capacity();
    const size_t mw_pos_stride = nw_new_old_pos.capacity();

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    map(always, to : w_accept_buffer_ptr[:nw_accept_index_ptrs.size()])")
    for (int i = 0; i < num_accepted; i++)
    {
      const int iw = reinterpret_cast<int*>(w_accept_buffer_ptr + sizeof(RealType*) * nw)[i];
      RealType* RSoA_dev_ptr = reinterpret_cast<RealType**>(w_accept_buffer_ptr)[i];
      for (int id = 0; id < QMCTraits::DIM; id++)
        RSoA_dev_ptr[iat + rsoa_stride * id] = mw_pos_ptr[iw + mw_pos_stride * id];
    }
  }

  const PosVectorSoa& getAllParticlePos() const override { return RSoA_hostview; }
  PosType getOneParticlePos(size_t iat) const override { return RSoA_hostview[iat]; }

  void donePbyP() override
  {
    if (is_host_position_changed_)
    {
      updateH2D();
      is_host_position_changed_ = false;
    }
  }

  const RealType* getDevicePtr() const { return RSoA.device_data(); }

private:
  ///particle positions in SoA layout
  VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> RSoA;

  ///one particle new/old positions in SoA layout
  VectorSoaContainer<RealType, QMCTraits::DIM * 2, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> nw_new_old_pos;

  ///nw_new_pos view of nw_new_old_pos
  VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> nw_active_pos_view;

  ///nw_old_pos view of nw_new_old_pos
  VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> nw_old_pos_view;

  /// accept list
  Vector<char, OMPallocator<char, PinnedAlignedAllocator<char>>> nw_accept_index_ptrs;

  ///host view of RSoA
  PosVectorSoa RSoA_hostview;

  ///if true, host position has been changed while device copy not updated.
  bool is_host_position_changed_;

  void updateH2D()
  {
    RealType* data = RSoA.data();
    PRAGMA_OFFLOAD("omp target update to(data[0:RSoA.capacity()*QMCTraits::DIM])")
    is_host_position_changed_ = false;
  }
};
} // namespace qmcplusplus
#endif
