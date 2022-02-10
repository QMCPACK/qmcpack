//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DTDIMPL_AB_OMPTARGET_H
#define QMCPLUSPLUS_DTDIMPL_AB_OMPTARGET_H

#include "Lattice/ParticleBConds3DSoa.h"
#include "DistanceTable.h"
#include "OMPTarget/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Particle/RealSpacePositionsOMPTarget.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a transposed form
 */
template<typename T, unsigned D, int SC>
class SoaDistanceTableABOMPTarget : public DTD_BConds<T, D, SC>, public DistanceTableAB
{
private:
  template<typename DT>
  using OffloadPinnedVector = Vector<DT, OMPallocator<DT, PinnedAlignedAllocator<DT>>>;

  ///accelerator output buffer for r and dr
  OffloadPinnedVector<RealType> r_dr_memorypool_;
  ///accelerator input array for a list of target particle positions, num_targets_ x D
  OffloadPinnedVector<T> target_pos;

  ///multi walker shared memory buffer
  struct DTABMultiWalkerMem : public Resource
  {
    ///accelerator output array for multiple walkers, [1+D][num_targets_][num_padded] (distances, displacements)
    OffloadPinnedVector<T> mw_r_dr;
    ///accelerator input buffer for multiple data set
    OffloadPinnedVector<char> offload_input;

    DTABMultiWalkerMem() : Resource("DTABMultiWalkerMem") {}

    DTABMultiWalkerMem(const DTABMultiWalkerMem&) : DTABMultiWalkerMem() {}

    Resource* makeClone() const override { return new DTABMultiWalkerMem(*this); }
  };

  std::unique_ptr<DTABMultiWalkerMem> mw_mem_;

  void resize()
  {
    if (num_sources_ * num_targets_ == 0)
      return;
    if (distances_.size())
      return;

    // initialize memory containers and views
    const size_t num_padded  = getAlignedSize<T>(num_sources_);
    const size_t stride_size = getPerTargetPctlStrideSize();
    r_dr_memorypool_.resize(stride_size * num_targets_);

    distances_.resize(num_targets_);
    displacements_.resize(num_targets_);
    for (int i = 0; i < num_targets_; ++i)
    {
      distances_[i].attachReference(r_dr_memorypool_.data() + i * stride_size, num_sources_);
      displacements_[i].attachReference(num_sources_, num_padded,
                                        r_dr_memorypool_.data() + i * stride_size + num_padded);
    }
  }

  static void associateResource(const RefVectorWithLeader<DistanceTable>& dt_list)
  {
    auto& dt_leader = dt_list.getCastedLeader<SoaDistanceTableABOMPTarget>();

    // initialize memory containers and views
    size_t count_targets = 0;
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableABOMPTarget>(iw);
      count_targets += dt.targets();
      dt.r_dr_memorypool_.free();
    }

    const size_t num_sources   = dt_leader.num_sources_;
    const size_t num_padded    = getAlignedSize<T>(dt_leader.num_sources_);
    const size_t stride_size   = num_padded * (D + 1);
    const size_t total_targets = count_targets;
    auto& mw_r_dr              = dt_leader.mw_mem_->mw_r_dr;
    mw_r_dr.resize(total_targets * stride_size);

    count_targets = 0;
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableABOMPTarget>(iw);
      assert(num_sources == dt.num_sources_);

      dt.distances_.resize(dt.targets());
      dt.displacements_.resize(dt.targets());

      for (int i = 0; i < dt.targets(); ++i)
      {
        dt.distances_[i].attachReference(mw_r_dr.data() + (i + count_targets) * stride_size, num_sources);
        dt.displacements_[i].attachReference(num_sources, num_padded,
                                             mw_r_dr.data() + (i + count_targets) * stride_size + num_padded);
      }
      count_targets += dt.targets();
    }
  }

public:
  SoaDistanceTableABOMPTarget(const ParticleSet& source, ParticleSet& target)
      : DTD_BConds<T, D, SC>(source.getLattice()),
        DistanceTableAB(source, target, DTModes::NEED_TEMP_DATA_ON_HOST),
        offload_timer_(*timer_manager.createTimer(std::string("DTABOMPTarget::offload_") + name_, timer_level_fine)),
        evaluate_timer_(*timer_manager.createTimer(std::string("DTABOMPTarget::evaluate_") + name_, timer_level_fine)),
        move_timer_(*timer_manager.createTimer(std::string("DTABOMPTarget::move_") + name_, timer_level_fine)),
        update_timer_(*timer_manager.createTimer(std::string("DTABOMPTarget::update_") + name_, timer_level_fine))

  {
    auto* coordinates_soa = dynamic_cast<const RealSpacePositionsOMPTarget*>(&source.getCoordinates());
    if (!coordinates_soa)
      throw std::runtime_error("Source particle set doesn't have OpenMP offload. Contact developers!");
    PRAGMA_OFFLOAD("omp target enter data map(to : this[:1])")

    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy in the update function
    // temp_r_ is padded explicitly while temp_dr_ is padded internally
    const int num_padded = getAlignedSize<T>(num_sources_);
    temp_r_.resize(num_padded);
    temp_dr_.resize(num_sources_);
  }

  SoaDistanceTableABOMPTarget()                                   = delete;
  SoaDistanceTableABOMPTarget(const SoaDistanceTableABOMPTarget&) = delete;

  ~SoaDistanceTableABOMPTarget() { PRAGMA_OFFLOAD("omp target exit data map(delete : this[:1])") }

  void createResource(ResourceCollection& collection) const override
  {
    auto resource_index = collection.addResource(std::make_unique<DTABMultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<DistanceTable>& dt_list) const override
  {
    auto res_ptr = dynamic_cast<DTABMultiWalkerMem*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error("SoaDistanceTableABOMPTarget::acquireResource dynamic_cast failed");
    auto& dt_leader = dt_list.getCastedLeader<SoaDistanceTableABOMPTarget>();
    dt_leader.mw_mem_.reset(res_ptr);
    associateResource(dt_list);
  }

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<DistanceTable>& dt_list) const override
  {
    collection.takebackResource(std::move(dt_list.getCastedLeader<SoaDistanceTableABOMPTarget>().mw_mem_));
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableABOMPTarget>(iw);
      dt.distances_.clear();
      dt.displacements_.clear();
    }
  }

  const T* getMultiWalkerDataPtr() const override
  {
    if (!mw_mem_)
      throw std::runtime_error("SoaDistanceTableABOMPTarget mw_mem_ is nullptr");
    return mw_mem_->mw_r_dr.data();
  }

  size_t getPerTargetPctlStrideSize() const override { return getAlignedSize<T>(num_sources_) * (D + 1); }

  /** evaluate the full table */
  inline void evaluate(ParticleSet& P) override
  {
    resize();

    ScopedTimer local_timer(evaluate_timer_);
    // be aware of the sign of Displacement
    const int num_targets_local = num_targets_;
    const int num_sources_local = num_sources_;
    const int num_padded        = getAlignedSize<T>(num_sources_);

    target_pos.resize(num_targets_ * D);
    for (size_t iat = 0; iat < num_targets_; iat++)
      for (size_t idim = 0; idim < D; idim++)
        target_pos[iat * D + idim] = P.R[iat][idim];

    auto* target_pos_ptr = target_pos.data();
    auto* source_pos_ptr = origin_.getCoordinates().getAllParticlePos().data();
    auto* r_dr_ptr       = distances_[0].data();
    assert(distances_[0].data() + num_padded == displacements_[0].data());

    // To maximize thread usage, the loop over electrons is chunked. Each chunk is sent to an OpenMP offload thread team.
    const int ChunkSizePerTeam = 256;
    const size_t num_teams     = (num_sources_ + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
    const size_t stride_size   = getPerTargetPctlStrideSize();

    {
      ScopedTimer offload(offload_timer_);
      PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(num_targets_*num_teams) \
                        map(to: source_pos_ptr[:num_padded*D]) \
                        map(always, to: target_pos_ptr[:num_targets_*D]) \
                        map(always, from: r_dr_ptr[:num_targets_*stride_size])")
      for (int iat = 0; iat < num_targets_local; ++iat)
        for (int team_id = 0; team_id < num_teams; team_id++)
        {
          const int first = ChunkSizePerTeam * team_id;
          const int last =
              (first + ChunkSizePerTeam) > num_sources_local ? num_sources_local : first + ChunkSizePerTeam;

          T pos[D];
          for (int idim = 0; idim < D; idim++)
            pos[idim] = target_pos_ptr[iat * D + idim];

          auto* r_iat_ptr  = r_dr_ptr + iat * stride_size;
          auto* dr_iat_ptr = r_iat_ptr + num_padded;

          PRAGMA_OFFLOAD("omp parallel for")
          for (int iel = first; iel < last; iel++)
            DTD_BConds<T, D, SC>::computeDistancesOffload(pos, source_pos_ptr, num_padded, r_iat_ptr, dr_iat_ptr,
                                                          num_padded, iel);
        }
    }
  }

  inline void mw_evaluate(const RefVectorWithLeader<DistanceTable>& dt_list,
                          const RefVectorWithLeader<ParticleSet>& p_list) const override
  {
    assert(this == &dt_list.getLeader());
    auto& dt_leader = dt_list.getCastedLeader<SoaDistanceTableABOMPTarget>();
    // multi walker resource must have been acquired
    assert(dt_leader.mw_mem_);

    ScopedTimer local_timer(evaluate_timer_);

    const size_t nw = dt_list.size();
    auto& mw_mem    = *dt_leader.mw_mem_;
    auto& mw_r_dr   = mw_mem.mw_r_dr;

    size_t count_targets = 0;
    for (ParticleSet& p : p_list)
      count_targets += p.getTotalNum();
    const size_t total_targets = count_targets;

    const int num_padded = getAlignedSize<T>(num_sources_);

#ifndef NDEBUG
    const int stride_size = getPerTargetPctlStrideSize();
    count_targets         = 0;
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableABOMPTarget>(iw);

      for (int i = 0; i < dt.targets(); ++i)
      {
        assert(dt.distances_[i].data() == mw_r_dr.data() + (i + count_targets) * stride_size);
        assert(dt.displacements_[i].data() == mw_r_dr.data() + (i + count_targets) * stride_size + num_padded);
      }
      count_targets += dt.targets();
    }
#endif

    // This is horrible optimization putting different data types in a single buffer but allows a single H2D transfer
    const size_t realtype_size = sizeof(RealType);
    const size_t int_size      = sizeof(int);
    const size_t ptr_size      = sizeof(RealType*);
    auto& offload_input        = mw_mem.offload_input;
    offload_input.resize(total_targets * D * realtype_size + total_targets * int_size + nw * ptr_size);
    auto source_ptrs      = reinterpret_cast<RealType**>(offload_input.data());
    auto target_positions = reinterpret_cast<RealType*>(offload_input.data() + ptr_size * nw);
    auto walker_id_ptr =
        reinterpret_cast<int*>(offload_input.data() + ptr_size * nw + total_targets * D * realtype_size);

    count_targets = 0;
    for (size_t iw = 0; iw < nw; iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableABOMPTarget>(iw);
      ParticleSet& pset(p_list[iw]);

      assert(dt.targets() == pset.getTotalNum());
      assert(num_sources_ == dt.num_sources_);

      auto& RSoA_OMPTarget = static_cast<const RealSpacePositionsOMPTarget&>(dt.origin_.getCoordinates());
      source_ptrs[iw]      = const_cast<RealType*>(RSoA_OMPTarget.getDevicePtr());

      for (size_t iat = 0; iat < pset.getTotalNum(); ++iat, ++count_targets)
      {
        walker_id_ptr[count_targets] = iw;
        for (size_t idim = 0; idim < D; idim++)
          target_positions[count_targets * D + idim] = pset.R[iat][idim];
      }
    }

    // To maximize thread usage, the loop over electrons is chunked. Each chunk is sent to an OpenMP offload thread team.
    const int ChunkSizePerTeam = 256;
    const size_t num_teams     = (num_sources_ + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    auto* r_dr_ptr              = mw_r_dr.data();
    auto* input_ptr             = offload_input.data();
    const int num_sources_local = num_sources_;

    {
      ScopedTimer offload(dt_leader.offload_timer_);
      PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(total_targets*num_teams) \
                        map(always, to: input_ptr[:offload_input.size()]) \
                        depend(out:r_dr_ptr[:mw_r_dr.size()]) nowait")
      for (int iat = 0; iat < total_targets; ++iat)
        for (int team_id = 0; team_id < num_teams; team_id++)
        {
          auto* target_pos_ptr = reinterpret_cast<RealType*>(input_ptr + ptr_size * nw);
          const int walker_id =
              reinterpret_cast<int*>(input_ptr + ptr_size * nw + total_targets * D * realtype_size)[iat];
          auto* source_pos_ptr = reinterpret_cast<RealType**>(input_ptr)[walker_id];
          auto* r_iat_ptr      = r_dr_ptr + iat * num_padded * (D + 1);
          auto* dr_iat_ptr     = r_dr_ptr + iat * num_padded * (D + 1) + num_padded;

          const int first = ChunkSizePerTeam * team_id;
          const int last =
              (first + ChunkSizePerTeam) > num_sources_local ? num_sources_local : first + ChunkSizePerTeam;

          T pos[D];
          for (int idim = 0; idim < D; idim++)
            pos[idim] = target_pos_ptr[iat * D + idim];

          PRAGMA_OFFLOAD("omp parallel for")
          for (int iel = first; iel < last; iel++)
            DTD_BConds<T, D, SC>::computeDistancesOffload(pos, source_pos_ptr, num_padded, r_iat_ptr, dr_iat_ptr,
                                                          num_padded, iel);
        }

      if (!(modes_ & DTModes::MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST))
      {
        PRAGMA_OFFLOAD(
            "omp target update from(r_dr_ptr[:mw_r_dr.size()]) depend(inout:r_dr_ptr[:mw_r_dr.size()]) nowait")
      }
      // wait for computing and (optional) transferring back to host.
      // It can potentially be moved to ParticleSet to fuse multiple similar taskwait
      PRAGMA_OFFLOAD("omp taskwait")
    }
  }

  inline void mw_recompute(const RefVectorWithLeader<DistanceTable>& dt_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           const std::vector<bool>& recompute) const override
  {
    mw_evaluate(dt_list, p_list);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);
    DTD_BConds<T, D, SC>::computeDistances(rnew, origin_.getCoordinates().getAllParticlePos(), temp_r_.data(), temp_dr_,
                                           0, num_sources_);
    // If the full table is not ready all the time, overwrite the current value.
    // If this step is missing, DT values can be undefined in case a move is rejected.
    if (!(modes_ & DTModes::NEED_FULL_TABLE_ANYTIME) && prepare_old)
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], origin_.getCoordinates().getAllParticlePos(),
                                             distances_[iat].data(), displacements_[iat], 0, num_sources_);
  }

  ///update the stripe for jat-th particle
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    std::copy_n(temp_r_.data(), num_sources_, distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(temp_dr_.data(idim), num_sources_, displacements_[iat].data(idim));
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override
  {
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < num_sources_; ++jat)
        if (temp_r_[jat] < min_dist)
        {
          min_dist = temp_r_[jat];
          index    = jat;
        }
      if (index >= 0)
      {
        r  = min_dist;
        dr = temp_dr_[index];
      }
    }
    else
    {
      for (int jat = 0; jat < num_sources_; ++jat)
        if (distances_[iat][jat] < min_dist)
        {
          min_dist = distances_[iat][jat];
          index    = jat;
        }
      if (index >= 0)
      {
        r  = min_dist;
        dr = displacements_[iat][index];
      }
    }
    assert(index >= 0 && index < num_sources_);
    return index;
  }

private:
  /// timer for offload portion
  NewTimer& offload_timer_;
  /// timer for evaluate()
  NewTimer& evaluate_timer_;
  /// timer for move()
  NewTimer& move_timer_;
  /// timer for update()
  NewTimer& update_timer_;
};
} // namespace qmcplusplus
#endif
