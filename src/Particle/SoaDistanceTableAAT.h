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
#ifndef QMCPLUSPLUS_DTDIMPL_AAT_H
#define QMCPLUSPLUS_DTDIMPL_AAT_H

#include "CPU/SIMD/algorithm.hpp"
#include "Lattice/ParticleBConds3DSoa.h"
#include "Particle/DistanceTableT.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAAT : public DTD_BConds<typename ParticleSetTraits<T>::RealType, D, SC>,
                             public DistanceTableAAT<T>
{
  using RealType  = typename DistanceTableAAT<T>::RealType;
  using PosType   = typename DistanceTableAAT<T>::PosType;
  using IndexType = typename DistanceTableAAT<T>::IndexType;

  /// actual memory for dist and displacements_
  aligned_vector<RealType> memory_pool_;

  SoaDistanceTableAAT(ParticleSetT<T>& target)
      : DTD_BConds<RealType, D, SC>(target.getLattice()),
        DistanceTableAAT<T>(target, DTModes::ALL_OFF),
        num_targets_padded_(getAlignedSize<RealType>(this->num_targets_)),
#if !defined(NDEBUG)
        old_prepared_elec_id_(-1),
#endif
        evaluate_timer_(createGlobalTimer(std::string("DTAA::evaluate_") + target.getName() + "_" + target.getName(),
                                          timer_level_fine)),
        move_timer_(createGlobalTimer(std::string("DTAA::move_") + target.getName() + "_" + target.getName(),
                                      timer_level_fine)),
        update_timer_(createGlobalTimer(std::string("DTAA::update_") + target.getName() + "_" + target.getName(),
                                        timer_level_fine))
  {
    resize();
  }

  SoaDistanceTableAAT()                           = delete;
  SoaDistanceTableAAT(const SoaDistanceTableAAT&) = delete;
  ~SoaDistanceTableAAT() override {}

  size_t compute_size(int N) const
  {
    const size_t num_padded = getAlignedSize<RealType>(N);
    const size_t Alignment  = getAlignment<RealType>();
    return (num_padded * (2 * N - num_padded + 1) + (Alignment - 1) * num_padded) / 2;
  }

  void resize()
  {
    // initialize memory containers and views
    const size_t total_size = compute_size(this->num_targets_);
    memory_pool_.resize(total_size * (1 + D));
    this->distances_.resize(this->num_targets_);
    this->displacements_.resize(this->num_targets_);
    for (int i = 0; i < this->num_targets_; ++i)
    {
      this->distances_[i].attachReference(memory_pool_.data() + compute_size(i), i);
      this->displacements_[i].attachReference(i, total_size, memory_pool_.data() + total_size + compute_size(i));
    }

    this->old_r_.resize(this->num_targets_);
    this->old_dr_.resize(this->num_targets_);
    this->temp_r_.resize(this->num_targets_);
    this->temp_dr_.resize(this->num_targets_);
  }

  inline void evaluate(ParticleSetT<T>& P) override
  {
    ScopedTimer local_timer(evaluate_timer_);
    constexpr RealType BigR = std::numeric_limits<RealType>::max();
    for (int iat = 1; iat < this->num_targets_; ++iat)
      DTD_BConds<RealType, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(),
                                                    this->distances_[iat].data(), this->displacements_[iat], 0, iat,
                                                    iat);
  }

  /// evaluate the temporary pair relations
  inline void move(const ParticleSetT<T>& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);

#if !defined(NDEBUG)
    old_prepared_elec_id_ = prepare_old ? iat : -1;
#endif
    DTD_BConds<RealType, D, SC>::computeDistances(rnew, P.getCoordinates().getAllParticlePos(), this->temp_r_.data(),
                                                  this->temp_dr_, 0, this->num_targets_, iat);
    // set up old_r_ and old_dr_ for moves may get accepted.
    if (prepare_old)
    {
      // recompute from scratch
      DTD_BConds<RealType, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(),
                                                    this->old_r_.data(), this->old_dr_, 0, this->num_targets_, iat);
      this->old_r_[iat] = std::numeric_limits<RealType>::max(); // assign a big number
    }
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override
  {
    // ensure there are neighbors
    assert(this->num_targets_ > 1);
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < this->num_targets_; ++jat)
        if (this->temp_r_[jat] < min_dist && jat != iat)
        {
          min_dist = this->temp_r_[jat];
          index    = jat;
        }
      assert(index >= 0);
      dr = this->temp_dr_[index];
    }
    else
    {
      for (int jat = 0; jat < iat; ++jat)
        if (this->distances_[iat][jat] < min_dist)
        {
          min_dist = this->distances_[iat][jat];
          index    = jat;
        }
      for (int jat = iat + 1; jat < this->num_targets_; ++jat)
        if (this->distances_[jat][iat] < min_dist)
        {
          min_dist = this->distances_[jat][iat];
          index    = jat;
        }
      assert(index != iat && index >= 0);
      if (index < iat)
        dr = this->displacements_[iat][index];
      else
        dr = this->displacements_[index][iat];
    }
    r = min_dist;
    return index;
  }

  /** After accepting the iat-th particle, update the iat-th row of distances_
     * and displacements_. Upper triangle is not needed in the later computation
     * and thus not updated
     */
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    // update [0, iat)
    const int nupdate = iat;
    // copy row
    assert(nupdate <= this->temp_r_.size());
    std::copy_n(this->temp_r_.data(), nupdate, this->distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(this->temp_dr_.data(idim), nupdate, this->displacements_[iat].data(idim));
    // copy column
    for (size_t i = iat + 1; i < this->num_targets_; ++i)
    {
      this->distances_[i][iat]     = this->temp_r_[i];
      this->displacements_[i](iat) = -this->temp_dr_[i];
    }
  }

  void updatePartial(IndexType jat, bool from_temp) override
  {
    ScopedTimer local_timer(update_timer_);
    // update [0, jat)
    const int nupdate = jat;
    if (from_temp)
    {
      // copy row
      assert(nupdate <= this->temp_r_.size());
      std::copy_n(this->temp_r_.data(), nupdate, this->distances_[jat].data());
      for (int idim = 0; idim < D; ++idim)
        std::copy_n(this->temp_dr_.data(idim), nupdate, this->displacements_[jat].data(idim));
    }
    else
    {
      assert(old_prepared_elec_id_ == jat);
      // copy row
      assert(nupdate <= this->old_r_.size());
      std::copy_n(this->old_r_.data(), nupdate, this->distances_[jat].data());
      for (int idim = 0; idim < D; ++idim)
        std::copy_n(this->old_dr_.data(idim), nupdate, this->displacements_[jat].data(idim));
    }
  }

private:
  /// number of targets with padding
  const size_t num_targets_padded_;
#if !defined(NDEBUG)
  /** set to particle id after move() with prepare_old = true. -1 means not
     * prepared. It is intended only for safety checks, not for codepath
     * selection.
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
