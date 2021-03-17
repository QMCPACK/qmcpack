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
#include "ResourceCollection.h"

namespace qmcplusplus
{
/** Introduced to handle virtual moves and ratio computations, e.g. for non-local PP evaluations.
   */
class RealSpacePositionsOMPTarget : public DynamicCoordinates
{
public:
  RealSpacePositionsOMPTarget()
      : DynamicCoordinates(DynamicCoordinateKind::DC_POS_OFFLOAD), is_host_position_changed_(false)
  {}
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

  size_t size() const override { return RSoA_hostview.size(); }

  void setAllParticlePos(const ParticlePos_t& R) override
  {
    resize(R.size());
    RSoA_hostview.copyIn(R);
    updateH2D();
    is_nw_new_pos_prepared = false;
  }

  void setOneParticlePos(const PosType& pos, size_t iat) override
  {
    RSoA_hostview(iat)        = pos;
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

  void mw_copyActivePos(const RefVectorWithLeader<DynamicCoordinates>& coords_list,
                        size_t iat,
                        const std::vector<PosType>& new_positions) const override
  {
    assert(this == &coords_list.getLeader());
    auto& coords_leader = coords_list.getCastedLeader<RealSpacePositionsOMPTarget>();
    // make this class unit tests friendly without the need of setup resources.
    if (!coords_leader.mw_mem_)
    {
      app_warning()
          << "RealSpacePositionsOMPTarget: This message should not be seen in production (performance bug) runs but "
             "only unit tests (expected)."
          << std::endl;
      coords_leader.mw_mem_ = std::make_unique<MultiWalkerMem>();
    }
    const auto nw = coords_list.size();

    auto& mw_new_pos = coords_leader.mw_mem_->mw_new_pos;
    mw_new_pos.resize(nw);

    for (int iw = 0; iw < nw; iw++)
      mw_new_pos(iw) = new_positions[iw];

    auto* mw_pos_ptr = mw_new_pos.data();
    PRAGMA_OFFLOAD("omp target update to(mw_pos_ptr[:QMCTraits::DIM * mw_new_pos.capacity()])")

    coords_leader.is_nw_new_pos_prepared = true;
  }

  void mw_acceptParticlePos(const RefVectorWithLeader<DynamicCoordinates>& coords_list,
                            size_t iat,
                            const std::vector<PosType>& new_positions,
                            const std::vector<bool>& isAccepted) const override
  {
    assert(this == &coords_list.getLeader());
    auto& coords_leader        = coords_list.getCastedLeader<RealSpacePositionsOMPTarget>();
    auto& mw_new_pos           = coords_leader.mw_mem_->mw_new_pos;
    auto& nw_accept_index_ptrs = coords_leader.mw_mem_->nw_accept_index_ptrs;
    const size_t nw            = coords_list.size();

    if (!is_nw_new_pos_prepared)
    {
      mw_copyActivePos(coords_list, iat, new_positions);
      app_warning() << "This message only appear in unit tests. Report a bug if seen in production code." << std::endl;
    }

    coords_leader.is_nw_new_pos_prepared = false;

    nw_accept_index_ptrs.resize((sizeof(int) + sizeof(RealType*)) * nw);
    auto* RSoA_ptr_array = reinterpret_cast<RealType**>(nw_accept_index_ptrs.data());
    auto* id_array       = reinterpret_cast<int*>(nw_accept_index_ptrs.data() + sizeof(RealType*) * coords_list.size());

    size_t num_accepted = 0;
    for (int iw = 0; iw < nw; iw++)
      if (isAccepted[iw])
      {
        auto& coords                 = coords_list.getCastedElement<RealSpacePositionsOMPTarget>(iw);
        RSoA_ptr_array[num_accepted] = coords.RSoA.device_data();
        id_array[num_accepted]       = iw;
        // save new coordinates on host copy
        coords.RSoA_hostview(iat) = mw_new_pos[iw];
        num_accepted++;
      }

    //offload to GPU
    auto* restrict w_accept_buffer_ptr = nw_accept_index_ptrs.data();
    auto* restrict mw_pos_ptr          = mw_new_pos.data();
    const size_t rsoa_stride           = RSoA.capacity();
    const size_t mw_pos_stride         = mw_new_pos.capacity();

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    map(always, to : w_accept_buffer_ptr[:nw_accept_index_ptrs.size()])")
    for (int i = 0; i < num_accepted; i++)
    {
      const int iw           = reinterpret_cast<int*>(w_accept_buffer_ptr + sizeof(RealType*) * nw)[i];
      RealType* RSoA_dev_ptr = reinterpret_cast<RealType**>(w_accept_buffer_ptr)[i];
      for (int id = 0; id < QMCTraits::DIM; id++)
        RSoA_dev_ptr[iat + rsoa_stride * id] = mw_pos_ptr[iw + mw_pos_stride * id];
    }
  }

  const PosVectorSoa& getAllParticlePos() const override { return RSoA_hostview; }
  PosType getOneParticlePos(size_t iat) const override { return RSoA_hostview[iat]; }

  void donePbyP() override
  {
    is_nw_new_pos_prepared = false;
    if (is_host_position_changed_)
    {
      updateH2D();
      is_host_position_changed_ = false;
    }
  }

  const RealType* getDevicePtr() const { return RSoA.device_data(); }

  const auto& getFusedNewPosBuffer() const { return mw_mem_->mw_new_pos; }

  void createResource(ResourceCollection& collection) const override
  {
    auto resource_index = collection.addResource(std::make_unique<MultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<DynamicCoordinates>& coords_list) const override
  {
    auto res_ptr = dynamic_cast<MultiWalkerMem*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error("RealSpacePositionsOMPTarget::acquireResource dynamic_cast failed");
    coords_list.getCastedLeader<RealSpacePositionsOMPTarget>().mw_mem_.reset(res_ptr);
  }

  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<DynamicCoordinates>& coords_list) const override
  {
    collection.takebackResource(std::move(coords_list.getCastedLeader<RealSpacePositionsOMPTarget>().mw_mem_));
  }

private:
  ///particle positions in SoA layout
  VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> RSoA;

  ///multi walker shared memory buffer
  struct MultiWalkerMem : public Resource
  {
    ///one particle new/old positions in SoA layout
    VectorSoaContainer<RealType, QMCTraits::DIM, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> mw_new_pos;

    /// accept list
    Vector<char, OMPallocator<char, PinnedAlignedAllocator<char>>> nw_accept_index_ptrs;

    MultiWalkerMem() : Resource("MultiWalkerMem") {}

    MultiWalkerMem(const MultiWalkerMem&) : MultiWalkerMem() {}

    Resource* makeClone() const override { return new MultiWalkerMem(*this); }
  };

  std::unique_ptr<MultiWalkerMem> mw_mem_;

  ///host view of RSoA
  PosVectorSoa RSoA_hostview;

  ///if true, host position has been changed while the device copy has not been updated.
  bool is_host_position_changed_;

  ///if true, mw_new_pos has been updated with active positions.
  bool is_nw_new_pos_prepared;

  void updateH2D()
  {
    RealType* data = RSoA.data();
    PRAGMA_OFFLOAD("omp target update to(data[0:RSoA.capacity()*QMCTraits::DIM])")
    is_host_position_changed_ = false;
  }
};
} // namespace qmcplusplus
#endif
