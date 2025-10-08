//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DTDIMPL_AA_H
#define QMCPLUSPLUS_DTDIMPL_AA_H

#include "Lattice/ParticleBConds3DSoa.h"
#include "DistanceTable.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAA : public DTD_BConds<T, D, SC>, public DistanceTableAA
{
  /// actual memory for dist and displacements_
  aligned_vector<RealType> memory_pool_;

  SoaDistanceTableAA(const ParticleSet& target)
      : DTD_BConds<T, D, SC>(target.getLattice()),
        DistanceTableAA(target, DTModes::ALL_OFF),
#if !defined(NDEBUG)
        old_prepared_elec_id_(-1),
#endif
        evaluate_timer_(createGlobalTimer("DTAA::evaluate_" + name_, timer_level_fine)),
        move_timer_(createGlobalTimer("DTAA::move_" + name_, timer_level_fine)),
        update_timer_(createGlobalTimer("DTAA::update_" + name_, timer_level_fine))
  {}

  SoaDistanceTableAA()                          = delete;
  SoaDistanceTableAA(const SoaDistanceTableAA&) = delete;
  ~SoaDistanceTableAA() override {}

  inline void evaluate(const DynamicCoordinates& coords) override
  {
    ScopedTimer local_timer(evaluate_timer_);
    if (num_targets_ != coords.size())
      resize(coords.size());

    constexpr T BigR = std::numeric_limits<T>::max();
    for (int iat = 1; iat < num_targets_; ++iat)
      DTD_BConds<T, D, SC>::computeDistances(coords.getOneParticlePos(iat), coords.getAllParticlePos(),
                                             distances_[iat].data(), displacements_[iat], 0, iat, iat);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);

#if !defined(NDEBUG)
    old_prepared_elec_id_ = prepare_old ? iat : -1;
#endif
    DTD_BConds<T, D, SC>::computeDistances(rnew, P.getCoordinates().getAllParticlePos(), temp_r_.data(), temp_dr_, 0,
                                           num_targets_, iat);
    // set up old_r_ and old_dr_ for moves may get accepted.
    if (prepare_old)
    {
      //recompute from scratch
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), old_r_.data(), old_dr_,
                                             0, num_targets_, iat);
      old_r_[iat] = std::numeric_limits<T>::max(); //assign a big number
    }
  }

  /** After accepting the iat-th particle, update the iat-th row of distances_ and displacements_.
   * Upper triangle is not needed in the later computation and thus not updated
   */
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    //update [0, iat)
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

  /// compute internal storage size
  size_t compute_size(int N) const
  {
    const size_t num_padded = getAlignedSize<T>(N);
    const size_t Alignment  = getAlignment<T>();
    return (num_padded * (2 * N - num_padded + 1) + (Alignment - 1) * num_padded) / 2;
  }

private:
  void resize(const size_t num_targets) override
  {
    // initialize memory containers and views
    num_targets_            = num_targets;
    const size_t total_size = compute_size(num_targets_);
    memory_pool_.resize(total_size * (1 + D));
    distances_.resize(num_targets_);
    displacements_.resize(num_targets_);
    for (int i = 0; i < num_targets_; ++i)
    {
      distances_[i].attachReference(memory_pool_.data() + compute_size(i), i);
      displacements_[i].attachReference(i, total_size, memory_pool_.data() + total_size + compute_size(i));
    }

    old_r_.resize(num_targets_);
    old_dr_.resize(num_targets_);
    temp_r_.resize(num_targets_);
    temp_dr_.resize(num_targets_);
  }

#if !defined(NDEBUG)
  /** set to particle id after move() with prepare_old = true. -1 means not prepared.
   * It is intended only for safety checks, not for codepath selection.
   */
  int old_prepared_elec_id_;
#endif
  /// timer for evaluate()
  NewTimer& evaluate_timer_;
  /// timer for move()
  NewTimer& move_timer_;
  /// timer for update()
  NewTimer& update_timer_;
};
} // namespace qmcplusplus
#endif
