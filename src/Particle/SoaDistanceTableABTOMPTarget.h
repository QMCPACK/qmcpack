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
#ifndef QMCPLUSPLUS_DTDIMPL_ABT_OMPTARGET_H
#define QMCPLUSPLUS_DTDIMPL_ABT_OMPTARGET_H

#include "DistanceTableT.h"
#include "Lattice/ParticleBConds3DSoa.h"
#include "OMPTarget/OMPTargetMath.hpp"
#include "OMPTarget/OMPallocator.hpp"
#include "Particle/RealSpacePositionsTOMPTarget.h"
#include "Platforms/PinnedAllocator.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a
 * transposed form
 */
template<typename T, unsigned D, int SC>
class SoaDistanceTableABTOMPTarget : public DTD_BConds<typename ParticleSetTraits<T>::RealType, D, SC>,
                                     public DistanceTableABT<T>
{
private:
  template<typename DT>
  using OffloadPinnedVector = Vector<DT, OMPallocator<DT, PinnedAlignedAllocator<DT>>>;

  using RealType  = typename DistanceTableABT<T>::RealType;
  using PosType   = typename DistanceTableABT<T>::PosType;
  using IndexType = typename DistanceTableABT<T>::IndexType;

  /// accelerator output buffer for r and dr
  OffloadPinnedVector<RealType> r_dr_memorypool_;
  /// accelerator input array for a list of target particle positions,
  /// num_targets_ x D
  OffloadPinnedVector<RealType> target_pos;

  /// multi walker shared memory buffer
  struct DTABMultiWalkerMem : public Resource
  {
    /// accelerator output array for multiple walkers,
    /// [1+D][num_targets_][num_padded] (distances, displacements)
    OffloadPinnedVector<RealType> mw_r_dr;
    /// accelerator input buffer for multiple data set
    OffloadPinnedVector<char> offload_input;

    DTABMultiWalkerMem() : Resource("DTABMultiWalkerMem") {}

    DTABMultiWalkerMem(const DTABMultiWalkerMem&) : DTABMultiWalkerMem() {}

    std::unique_ptr<Resource> makeClone() const override { return std::make_unique<DTABMultiWalkerMem>(*this); }
  };

  ResourceHandle<DTABMultiWalkerMem> mw_mem_handle_;

  void resize()
  {
    if (this->num_sources_ * this->num_targets_ == 0)
      return;
    if (this->distances_.size())
      return;

    // initialize memory containers and views
    const size_t num_padded  = getAlignedSize<RealType>(this->num_sources_);
    const size_t stride_size = getPerTargetPctlStrideSize();
    r_dr_memorypool_.resize(stride_size * this->num_targets_);

    this->distances_.resize(this->num_targets_);
    this->displacements_.resize(this->num_targets_);
    for (int i = 0; i < this->num_targets_; ++i)
    {
      this->distances_[i].attachReference(r_dr_memorypool_.data() + i * stride_size, this->num_sources_);
      this->displacements_[i].attachReference(this->num_sources_, num_padded,
                                              r_dr_memorypool_.data() + i * stride_size + num_padded);
    }
  }

  static void associateResource(const RefVectorWithLeader<DistanceTableT<T>>& dt_list)
  {
    auto& dt_leader = dt_list.template getCastedLeader<SoaDistanceTableABTOMPTarget>();

    // initialize memory containers and views
    size_t count_targets = 0;
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.template getCastedElement<SoaDistanceTableABTOMPTarget>(iw);
      count_targets += dt.targets();
      dt.r_dr_memorypool_.free();
    }

    const size_t num_sources   = dt_leader.num_sources_;
    const size_t num_padded    = getAlignedSize<RealType>(dt_leader.num_sources_);
    const size_t stride_size   = num_padded * (D + 1);
    const size_t total_targets = count_targets;
    auto& mw_r_dr              = dt_leader.mw_mem_handle_.getResource().mw_r_dr;
    mw_r_dr.resize(total_targets * stride_size);

    count_targets = 0;
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.template getCastedElement<SoaDistanceTableABTOMPTarget>(iw);
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
  SoaDistanceTableABTOMPTarget(const ParticleSetT<T>& source, ParticleSetT<T>& target)
      : DTD_BConds<RealType, D, SC>(source.getLattice()),
        DistanceTableABT<T>(source, target, DTModes::ALL_OFF),
        offload_timer_(createGlobalTimer(std::string("DTABOMPTarget::offload_") + this->name_, timer_level_fine)),
        evaluate_timer_(createGlobalTimer(std::string("DTABOMPTarget::evaluate_") + this->name_, timer_level_fine)),
        move_timer_(createGlobalTimer(std::string("DTABOMPTarget::move_") + this->name_, timer_level_fine)),
        update_timer_(createGlobalTimer(std::string("DTABOMPTarget::update_") + this->name_, timer_level_fine))

  {
    auto* coordinates_soa = dynamic_cast<const RealSpacePositionsTOMPTarget<T>*>(&source.getCoordinates());
    if (!coordinates_soa)
      throw std::runtime_error("Source particle set doesn't have OpenMP "
                               "offload. Contact developers!");
    PRAGMA_OFFLOAD("omp target enter data map(to : this[:1])")

    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy
    // in the update function temp_r_ is padded explicitly while temp_dr_ is
    // padded internally
    const int num_padded = getAlignedSize<RealType>(this->num_sources_);
    this->temp_r_.resize(num_padded);
    this->temp_dr_.resize(this->num_sources_);
  }

  SoaDistanceTableABTOMPTarget()                                    = delete;
  SoaDistanceTableABTOMPTarget(const SoaDistanceTableABTOMPTarget&) = delete;

  ~SoaDistanceTableABTOMPTarget() { PRAGMA_OFFLOAD("omp target exit data map(delete : this[:1])") }

  void createResource(ResourceCollection& collection) const override
  {
    auto resource_index = collection.addResource(std::make_unique<DTABMultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<DistanceTableT<T>>& dt_list) const override
  {
    auto& dt_leader          = dt_list.template getCastedLeader<SoaDistanceTableABTOMPTarget>();
    dt_leader.mw_mem_handle_ = collection.lendResource<DTABMultiWalkerMem>();
    associateResource(dt_list);
  }

  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<DistanceTableT<T>>& dt_list) const override
  {
    collection.takebackResource(dt_list.template getCastedLeader<SoaDistanceTableABTOMPTarget>().mw_mem_handle_);
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.template getCastedElement<SoaDistanceTableABTOMPTarget>(iw);
      dt.distances_.clear();
      dt.displacements_.clear();
    }
  }

  const RealType* getMultiWalkerDataPtr() const override { return mw_mem_handle_.getResource().mw_r_dr.data(); }

  size_t getPerTargetPctlStrideSize() const override { return getAlignedSize<RealType>(this->num_sources_) * (D + 1); }

  /** evaluate the full table */
  inline void evaluate(ParticleSetT<T>& P) override
  {
    resize();

    ScopedTimer local_timer(evaluate_timer_);
    // be aware of the sign of Displacement
    const int num_targets_local = this->num_targets_;
    const int num_sources_local = this->num_sources_;
    const int num_padded        = getAlignedSize<RealType>(this->num_sources_);

    target_pos.resize(this->num_targets_ * D);
    for (size_t iat = 0; iat < this->num_targets_; iat++)
      for (size_t idim = 0; idim < D; idim++)
        target_pos[iat * D + idim] = P.R[iat][idim];

    auto* target_pos_ptr = target_pos.data();
    auto* source_pos_ptr = this->origin_.getCoordinates().getAllParticlePos().data();
    auto* r_dr_ptr       = this->distances_[0].data();
    assert(this->distances_[0].data() + num_padded == this->displacements_[0].data());

    // To maximize thread usage, the loop over electrons is chunked. Each
    // chunk is sent to an OpenMP offload thread team.
    const int ChunkSizePerTeam = 512;
    const size_t num_teams     = (this->num_sources_ + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
    const size_t stride_size   = getPerTargetPctlStrideSize();

    {
      ScopedTimer offload(offload_timer_);
      PRAGMA_OFFLOAD("omp target teams distribute collapse(2) \
                num_teams(this->num_targets_*num_teams) \
                map(to: source_pos_ptr[:num_padded*D]) \
                map(always, to: target_pos_ptr[:this->num_targets_*D]) \
                map(always, from: r_dr_ptr[:this->num_targets_*stride_size])")
      for (int iat = 0; iat < num_targets_local; ++iat)
        for (int team_id = 0; team_id < num_teams; team_id++)
        {
          const int first = ChunkSizePerTeam * team_id;
          const int last  = omptarget::min(first + ChunkSizePerTeam, num_sources_local);

          RealType pos[D];
          for (int idim = 0; idim < D; idim++)
            pos[idim] = target_pos_ptr[iat * D + idim];

          auto* r_iat_ptr  = r_dr_ptr + iat * stride_size;
          auto* dr_iat_ptr = r_iat_ptr + num_padded;

          PRAGMA_OFFLOAD("omp parallel for")
          for (int iel = first; iel < last; iel++)
            DTD_BConds<RealType, D, SC>::computeDistancesOffload(pos, source_pos_ptr, num_padded, r_iat_ptr, dr_iat_ptr,
                                                                 num_padded, iel);
        }
    }
  }

  inline void mw_evaluate(const RefVectorWithLeader<DistanceTableT<T>>& dt_list,
                          const RefVectorWithLeader<ParticleSetT<T>>& p_list) const override
  {
    assert(this == &dt_list.getLeader());
    auto& dt_leader = dt_list.template getCastedLeader<SoaDistanceTableABTOMPTarget>();

    ScopedTimer local_timer(evaluate_timer_);

    const size_t nw            = dt_list.size();
    DTABMultiWalkerMem& mw_mem = dt_leader.mw_mem_handle_;
    auto& mw_r_dr              = mw_mem.mw_r_dr;

    size_t count_targets = 0;
    for (ParticleSetT<T>& p : p_list)
      count_targets += p.getTotalNum();
    const size_t total_targets = count_targets;

    const int num_padded = getAlignedSize<RealType>(this->num_sources_);

#ifndef NDEBUG
    const int stride_size = getPerTargetPctlStrideSize();
    count_targets         = 0;
    for (size_t iw = 0; iw < dt_list.size(); iw++)
    {
      auto& dt = dt_list.template getCastedElement<SoaDistanceTableABTOMPTarget>(iw);

      for (int i = 0; i < dt.targets(); ++i)
      {
        assert(dt.distances_[i].data() == mw_r_dr.data() + (i + count_targets) * stride_size);
        assert(dt.displacements_[i].data() == mw_r_dr.data() + (i + count_targets) * stride_size + num_padded);
      }
      count_targets += dt.targets();
    }
#endif

    // This is horrible optimization putting different data types in a
    // single buffer but allows a single H2D transfer
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
      auto& dt = dt_list.template getCastedElement<SoaDistanceTableABTOMPTarget>(iw);
      ParticleSetT<T>& pset(p_list[iw]);

      assert(dt.targets() == pset.getTotalNum());
      assert(this->num_sources_ == dt.num_sources_);

      auto& RSoA_OMPTarget = static_cast<const RealSpacePositionsTOMPTarget<T>&>(dt.origin_.getCoordinates());
      source_ptrs[iw]      = const_cast<RealType*>(RSoA_OMPTarget.getDevicePtr());

      for (size_t iat = 0; iat < pset.getTotalNum(); ++iat, ++count_targets)
      {
        walker_id_ptr[count_targets] = iw;
        for (size_t idim = 0; idim < D; idim++)
          target_positions[count_targets * D + idim] = pset.R[iat][idim];
      }
    }

    // To maximize thread usage, the loop over electrons is chunked. Each
    // chunk is sent to an OpenMP offload thread team.
    const int ChunkSizePerTeam = 512;
    const size_t num_teams     = (this->num_sources_ + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    auto* r_dr_ptr              = mw_r_dr.data();
    auto* input_ptr             = offload_input.data();
    const int num_sources_local = this->num_sources_;

    {
      ScopedTimer offload(dt_leader.offload_timer_);
      PRAGMA_OFFLOAD("omp target teams distribute collapse(2) \
                num_teams(total_targets*num_teams) \
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
          const int last  = omptarget::min(first + ChunkSizePerTeam, num_sources_local);

          RealType pos[D];
          for (int idim = 0; idim < D; idim++)
            pos[idim] = target_pos_ptr[iat * D + idim];

          PRAGMA_OFFLOAD("omp parallel for")
          for (int iel = first; iel < last; iel++)
            DTD_BConds<RealType, D, SC>::computeDistancesOffload(pos, source_pos_ptr, num_padded, r_iat_ptr, dr_iat_ptr,
                                                                 num_padded, iel);
        }

      if (!(this->modes_ & DTModes::MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST))
      {
        PRAGMA_OFFLOAD("omp target update from(r_dr_ptr[:mw_r_dr.size()]) \
                    depend(inout:r_dr_ptr[:mw_r_dr.size()]) nowait")
      }
      // wait for computing and (optional) transferring back to host.
      // It can potentially be moved to ParticleSet to fuse multiple
      // similar taskwait
      PRAGMA_OFFLOAD("omp taskwait")
    }
  }

  inline void mw_recompute(const RefVectorWithLeader<DistanceTableT<T>>& dt_list,
                           const RefVectorWithLeader<ParticleSetT<T>>& p_list,
                           const std::vector<bool>& recompute) const override
  {
    mw_evaluate(dt_list, p_list);
  }

  /// evaluate the temporary pair relations
  inline void move(const ParticleSetT<T>& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);
    DTD_BConds<RealType, D, SC>::computeDistances(rnew, this->origin_.getCoordinates().getAllParticlePos(),
                                                  this->temp_r_.data(), this->temp_dr_, 0, this->num_sources_);
    // If the full table is not ready all the time, overwrite the current
    // value. If this step is missing, DT values can be undefined in case a
    // move is rejected.
    if (!(this->modes_ & DTModes::NEED_FULL_TABLE_ANYTIME) && prepare_old)
      DTD_BConds<RealType, D, SC>::computeDistances(P.R[iat], this->origin_.getCoordinates().getAllParticlePos(),
                                                    this->distances_[iat].data(), this->displacements_[iat], 0,
                                                    this->num_sources_);
  }

  /// update the stripe for jat-th particle
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    std::copy_n(this->temp_r_.data(), this->num_sources_, this->distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(this->temp_dr_.data(idim), this->num_sources_, this->displacements_[iat].data(idim));
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override
  {
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < this->num_sources_; ++jat)
        if (this->temp_r_[jat] < min_dist)
        {
          min_dist = this->temp_r_[jat];
          index    = jat;
        }
      if (index >= 0)
      {
        r  = min_dist;
        dr = this->temp_dr_[index];
      }
    }
    else
    {
      for (int jat = 0; jat < this->num_sources_; ++jat)
        if (this->distances_[iat][jat] < min_dist)
        {
          min_dist = this->distances_[iat][jat];
          index    = jat;
        }
      if (index >= 0)
      {
        r  = min_dist;
        dr = this->displacements_[iat][index];
      }
    }
    assert(index >= 0 && index < this->num_sources_);
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
