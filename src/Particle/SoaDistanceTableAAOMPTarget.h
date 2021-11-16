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
#ifndef QMCPLUSPLUS_DTDIMPL_AA_OMPTARGET_H
#define QMCPLUSPLUS_DTDIMPL_AA_OMPTARGET_H

#include "Lattice/ParticleBConds3DSoa.h"
#include "DistanceTable.h"
#include "CPU/SIMD/algorithm.hpp"
#include "OMPTarget/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Particle/RealSpacePositionsOMPTarget.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAAOMPTarget : public DTD_BConds<T, D, SC>, public DistanceTableAA
{
  /// actual memory for dist and displacements_
  aligned_vector<RealType> memory_pool_;

  /// actual memory for temp_r_
  DistRow temp_r_mem_;
  /// actual memory for temp_dr_
  DisplRow temp_dr_mem_;
  /// actual memory for old_r_
  DistRow old_r_mem_;
  /// actual memory for old_dr_
  DisplRow old_dr_mem_;

  ///multi walker shared memory buffer
  struct DTAAMultiWalkerMem : public Resource
  {
    ///dist displ for temporary and old pairs
    Vector<RealType, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> nw_new_old_dist_displ;

    Vector<const RealType*, OMPallocator<const RealType*, PinnedAlignedAllocator<const RealType*>>> rsoa_dev_list;

    DTAAMultiWalkerMem() : Resource("DTAAMultiWalkerMem") {}

    DTAAMultiWalkerMem(const DTAAMultiWalkerMem&) : DTAAMultiWalkerMem() {}

    Resource* makeClone() const override { return new DTAAMultiWalkerMem(*this); }
  };

  std::unique_ptr<DTAAMultiWalkerMem> mw_mem_;

  SoaDistanceTableAAOMPTarget(ParticleSet& target)
      : DTD_BConds<T, D, SC>(target.Lattice),
        DistanceTableAA(target, DTModes::ALL_OFF),
        num_targets_padded_(getAlignedSize<T>(num_targets_)),
#if !defined(NDEBUG)
        old_prepared_elec_id_(-1),
#endif
        offload_timer_(
            *timer_manager.createTimer(std::string("SoaDistanceTableAAOMPTarget::offload_") + name_, timer_level_fine)),
        evaluate_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAAOMPTarget::evaluate_") + name_,
                                                   timer_level_fine)),
        move_timer_(
            *timer_manager.createTimer(std::string("SoaDistanceTableAAOMPTarget::move_") + name_, timer_level_fine)),
        update_timer_(
            *timer_manager.createTimer(std::string("SoaDistanceTableAAOMPTarget::update_") + name_, timer_level_fine))

  {
    auto* coordinates_soa = dynamic_cast<const RealSpacePositionsOMPTarget*>(&target.getCoordinates());
    if (!coordinates_soa)
      throw std::runtime_error("Source particle set doesn't have OpenMP offload. Contact developers!");
    resize();
    PRAGMA_OFFLOAD("omp target enter data map(to : this[:1])")
  }

  SoaDistanceTableAAOMPTarget()                                   = delete;
  SoaDistanceTableAAOMPTarget(const SoaDistanceTableAAOMPTarget&) = delete;
  ~SoaDistanceTableAAOMPTarget(){PRAGMA_OFFLOAD("omp target exit data map(delete : this[:1])")}

  size_t compute_size(int N) const
  {
    const size_t num_padded = getAlignedSize<T>(N);
    const size_t Alignment  = getAlignment<T>();
    return (num_padded * (2 * N - num_padded + 1) + (Alignment - 1) * num_padded) / 2;
  }

  void resize()
  {
    // initialize memory containers and views
    const size_t total_size = compute_size(num_targets_);
    memory_pool_.resize(total_size * (1 + D));
    distances_.resize(num_targets_);
    displacements_.resize(num_targets_);
    for (int i = 0; i < num_targets_; ++i)
    {
      distances_[i].attachReference(memory_pool_.data() + compute_size(i), i);
      displacements_[i].attachReference(i, total_size, memory_pool_.data() + total_size + compute_size(i));
    }

    old_r_mem_.resize(num_targets_);
    old_dr_mem_.resize(num_targets_);
    temp_r_mem_.resize(num_targets_);
    temp_dr_mem_.resize(num_targets_);
  }

  const RealType* getMultiWalkerTempDataPtr() const override
  {
    if (!mw_mem_)
      throw std::runtime_error("SoaDistanceTableAAOMPTarget mw_mem_ is nullptr");
    return mw_mem_->nw_new_old_dist_displ.data();
  }

  void createResource(ResourceCollection& collection) const override
  {
    auto resource_index = collection.addResource(std::make_unique<DTAAMultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<DistanceTable>& dt_list) const override
  {
    auto res_ptr = dynamic_cast<DTAAMultiWalkerMem*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error("SoaDistanceTableAAOMPTarget::acquireResource dynamic_cast failed");
    assert(this == &dt_list.getLeader());
    auto& dt_leader = dt_list.getCastedLeader<SoaDistanceTableAAOMPTarget>();
    dt_leader.mw_mem_.reset(res_ptr);
    auto& mw_mem             = *dt_leader.mw_mem_;
    const size_t nw          = dt_list.size();
    const size_t stride_size = num_targets_padded_ * (D + 1);

    for (int iw = 0; iw < nw; iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableAAOMPTarget>(iw);
      dt.temp_r_.free();
      dt.temp_dr_.free();
      dt.old_r_.free();
      dt.old_dr_.free();
    }

    auto& nw_new_old_dist_displ = mw_mem.nw_new_old_dist_displ;
    nw_new_old_dist_displ.resize(nw * 2 * stride_size);
    for (int iw = 0; iw < nw; iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableAAOMPTarget>(iw);
      dt.temp_r_.attachReference(nw_new_old_dist_displ.data() + stride_size * iw, num_targets_padded_);
      dt.temp_dr_.attachReference(num_targets_, num_targets_padded_,
                                  nw_new_old_dist_displ.data() + stride_size * iw + num_targets_padded_);
      dt.old_r_.attachReference(nw_new_old_dist_displ.data() + stride_size * (iw + nw), num_targets_padded_);
      dt.old_dr_.attachReference(num_targets_, num_targets_padded_,
                                 nw_new_old_dist_displ.data() + stride_size * (iw + nw) + num_targets_padded_);
    }
  }

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<DistanceTable>& dt_list) const override
  {
    collection.takebackResource(std::move(dt_list.getCastedLeader<SoaDistanceTableAAOMPTarget>().mw_mem_));
    const size_t nw = dt_list.size();
    for (int iw = 0; iw < nw; iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableAAOMPTarget>(iw);
      dt.temp_r_.free();
      dt.temp_dr_.free();
      dt.old_r_.free();
      dt.old_dr_.free();
    }
  }

  inline void evaluate(ParticleSet& P) override
  {
    ScopedTimer local_timer(evaluate_timer_);

    constexpr T BigR = std::numeric_limits<T>::max();
    for (int iat = 1; iat < num_targets_; ++iat)
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), distances_[iat].data(),
                                             displacements_[iat], 0, iat, iat);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);

#if !defined(NDEBUG)
    old_prepared_elec_id_ = prepare_old ? iat : -1;
#endif
    temp_r_.attachReference(temp_r_mem_.data(), temp_r_mem_.size());
    temp_dr_.attachReference(temp_dr_mem_.size(), temp_dr_mem_.capacity(), temp_dr_mem_.data());

    assert((prepare_old && iat >= 0 && iat < num_targets_) || !prepare_old);
    DTD_BConds<T, D, SC>::computeDistances(rnew, P.getCoordinates().getAllParticlePos(), temp_r_.data(), temp_dr_, 0,
                                           num_targets_, iat);
    // set up old_r_ and old_dr_ for moves may get accepted.
    if (prepare_old)
    {
      old_r_.attachReference(old_r_mem_.data(), old_r_mem_.size());
      old_dr_.attachReference(old_dr_mem_.size(), old_dr_mem_.capacity(), old_dr_mem_.data());
      //recompute from scratch
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), old_r_.data(), old_dr_,
                                             0, num_targets_, iat);
      old_r_[iat] = std::numeric_limits<T>::max(); //assign a big number
    }
  }

  /** evaluate the temporary pair relations when a move is proposed
   * this implementation is asynchronous and the synchronization is managed at ParticleSet.
   * Transfering results to host depends on DTModes::NEED_TEMP_DATA_ON_HOST.
   * If the temporary pair distance are consumed on the device directly, the device to host data transfer can be
   * skipped as an optimization.
   */
  void mw_move(const RefVectorWithLeader<DistanceTable>& dt_list,
               const RefVectorWithLeader<ParticleSet>& p_list,
               const std::vector<PosType>& rnew_list,
               const IndexType iat = 0,
               bool prepare_old    = true) const override
  {
    assert(this == &dt_list.getLeader());
    auto& dt_leader = dt_list.getCastedLeader<SoaDistanceTableAAOMPTarget>();
    // multi walker resource must have been acquired;
    assert(dt_leader.mw_mem_);
    auto& mw_mem      = *dt_leader.mw_mem_;
    auto& pset_leader = p_list.getLeader();

    ScopedTimer local_timer(move_timer_);
    const size_t nw          = dt_list.size();
    const size_t stride_size = num_targets_padded_ * (D + 1);

    auto& nw_new_old_dist_displ = mw_mem.nw_new_old_dist_displ;
    auto& rsoa_dev_list         = mw_mem.rsoa_dev_list;
    rsoa_dev_list.resize(nw);

    for (int iw = 0; iw < nw; iw++)
    {
      auto& dt = dt_list.getCastedElement<SoaDistanceTableAAOMPTarget>(iw);
#if !defined(NDEBUG)
      dt.old_prepared_elec_id_ = prepare_old ? iat : -1;
#endif
      auto& coordinates_soa = static_cast<const RealSpacePositionsOMPTarget&>(p_list[iw].getCoordinates());
      rsoa_dev_list[iw]     = coordinates_soa.getDevicePtr();
    }

    const int ChunkSizePerTeam = 256;
    const size_t num_teams     = (num_targets_ + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    auto& coordinates_leader = static_cast<const RealSpacePositionsOMPTarget&>(pset_leader.getCoordinates());

    const auto activePtcl_local  = pset_leader.activePtcl;
    const auto num_sources_local = num_targets_;
    const auto num_padded        = num_targets_padded_;
    auto* rsoa_dev_list_ptr      = rsoa_dev_list.data();
    auto* r_dr_ptr               = nw_new_old_dist_displ.data();
    auto* new_pos_ptr            = coordinates_leader.getFusedNewPosBuffer().data();
    const size_t new_pos_stride  = coordinates_leader.getFusedNewPosBuffer().capacity();

    {
      ScopedTimer offload(offload_timer_);
      PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(nw * num_teams) \
                        map(always, to: rsoa_dev_list_ptr[:rsoa_dev_list.size()]) \
                        nowait depend(out: r_dr_ptr[:nw_new_old_dist_displ.size()])")
      for (int iw = 0; iw < nw; ++iw)
        for (int team_id = 0; team_id < num_teams; team_id++)
        {
          auto* source_pos_ptr = rsoa_dev_list_ptr[iw];
          const int first      = ChunkSizePerTeam * team_id;
          const int last =
              (first + ChunkSizePerTeam) > num_sources_local ? num_sources_local : first + ChunkSizePerTeam;

          { // temp
            auto* r_iw_ptr  = r_dr_ptr + iw * stride_size;
            auto* dr_iw_ptr = r_dr_ptr + iw * stride_size + num_padded;

            T pos[D];
            for (int idim = 0; idim < D; idim++)
              pos[idim] = new_pos_ptr[idim * new_pos_stride + iw];

            PRAGMA_OFFLOAD("omp parallel for")
            for (int iel = first; iel < last; iel++)
              DTD_BConds<T, D, SC>::computeDistancesOffload(pos, source_pos_ptr, num_padded, r_iw_ptr, dr_iw_ptr,
                                                            num_padded, iel, activePtcl_local);
          }

          if (prepare_old)
          { // old
            auto* r_iw_ptr  = r_dr_ptr + (iw + nw) * stride_size;
            auto* dr_iw_ptr = r_dr_ptr + (iw + nw) * stride_size + num_padded;

            T pos[D];
            for (int idim = 0; idim < D; idim++)
              pos[idim] = source_pos_ptr[idim * num_padded + iat];

            PRAGMA_OFFLOAD("omp parallel for")
            for (int iel = first; iel < last; iel++)
              DTD_BConds<T, D, SC>::computeDistancesOffload(pos, source_pos_ptr, num_padded, r_iw_ptr, dr_iw_ptr,
                                                            num_padded, iel, iat);
            r_iw_ptr[iat] = std::numeric_limits<T>::max(); //assign a big number
          }
        }
    }

    if (modes_ & DTModes::NEED_TEMP_DATA_ON_HOST)
    {
      PRAGMA_OFFLOAD("omp target update nowait depend(inout: r_dr_ptr[:nw_new_old_dist_displ.size()]) \
                      from(r_dr_ptr[:nw_new_old_dist_displ.size()])")
    }
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override
  {
    //ensure there are neighbors
    assert(num_targets_ > 1);
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < num_targets_; ++jat)
        if (temp_r_[jat] < min_dist && jat != iat)
        {
          min_dist = temp_r_[jat];
          index    = jat;
        }
      assert(index >= 0);
      dr = temp_dr_[index];
    }
    else
    {
      for (int jat = 0; jat < iat; ++jat)
        if (distances_[iat][jat] < min_dist)
        {
          min_dist = distances_[iat][jat];
          index    = jat;
        }
      for (int jat = iat + 1; jat < num_targets_; ++jat)
        if (distances_[jat][iat] < min_dist)
        {
          min_dist = distances_[jat][iat];
          index    = jat;
        }
      assert(index != iat && index >= 0);
      if (index < iat)
        dr = displacements_[iat][index];
      else
        dr = displacements_[index][iat];
    }
    r = min_dist;
    return index;
  }

  /** After accepting the iat-th particle, update the iat-th row of distances_ and displacements_.
   * Upper triangle is not needed in the later computation and thus not updated
   */
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    //update [0, iat) columns
    const int nupdate = iat;
    //copy row
    assert(nupdate <= temp_r_.size());
    std::copy_n(temp_r_.data(), nupdate, distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(temp_dr_.data(idim), nupdate, displacements_[iat].data(idim));
    //copy column
    for (size_t i = iat + 1; i < num_targets_; ++i)
    {
      distances_[i][iat]     = temp_r_[i];
      displacements_[i](iat) = -temp_dr_[i];
    }
  }

  void updatePartial(IndexType jat, bool from_temp) override
  {
    ScopedTimer local_timer(update_timer_);

    //update [0, jat)
    const int nupdate = jat;
    if (from_temp)
    {
      //copy row
      assert(nupdate <= temp_r_.size());
      std::copy_n(temp_r_.data(), nupdate, distances_[jat].data());
      for (int idim = 0; idim < D; ++idim)
        std::copy_n(temp_dr_.data(idim), nupdate, displacements_[jat].data(idim));
    }
    else
    {
      assert(old_prepared_elec_id_ == jat);
      //copy row
      assert(nupdate <= old_r_.size());
      std::copy_n(old_r_.data(), nupdate, distances_[jat].data());
      for (int idim = 0; idim < D; ++idim)
        std::copy_n(old_dr_.data(idim), nupdate, displacements_[jat].data(idim));
    }
  }

  void mw_updatePartial(const RefVectorWithLeader<DistanceTable>& dt_list,
                        IndexType jat,
                        const std::vector<bool>& from_temp) override
  {
    // if temp data on host is not updated by mw_move during p-by-p moves, there is no need to update distance table
    if (!(modes_ & DTModes::NEED_TEMP_DATA_ON_HOST))
      return;

#pragma omp parallel for
    for (int iw = 0; iw < dt_list.size(); iw++)
      dt_list[iw].updatePartial(jat, from_temp[iw]);
  }

  void mw_finalizePbyP(const RefVectorWithLeader<DistanceTable>& dt_list,
                       const RefVectorWithLeader<ParticleSet>& p_list) const override
  {
    // if the distance table is not updated by mw_move during p-by-p, needs to recompute the whole table
    // before being used by Hamiltonian.
    if (!(modes_ & DTModes::NEED_TEMP_DATA_ON_HOST))
      mw_evaluate(dt_list, p_list);
  }

private:
  ///number of targets with padding
  const size_t num_targets_padded_;
#if !defined(NDEBUG)
  /** set to particle id after move() with prepare_old = true. -1 means not prepared.
   * It is intended only for safety checks, not for codepath selection.
   */
  int old_prepared_elec_id_;
#endif
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
