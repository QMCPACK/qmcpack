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
#include "DistanceTableData.h"
#include "CPU/SIMD/algorithm.hpp"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAA : public DTD_BConds<T, D, SC>, public DistanceTableData
{
  ///number of targets with padding
  int Ntargets_padded;

  ///actual memory for dist and displacements_
  aligned_vector<RealType> memory_pool_;

  /// old distances
  DistRow old_r_;

  /// old displacements
  DisplRow old_dr_;

  SoaDistanceTableAA(ParticleSet& target)
      : DTD_BConds<T, D, SC>(target.Lattice),
        DistanceTableData(target, target, DTModes::NEED_TEMP_DATA_ON_HOST),
        evaluate_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAA::evaluate_") + target.getName() +
                                                       "_" + target.getName(),
                                                   timer_level_fine)),
        move_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAA::move_") + target.getName() + "_" +
                                                   target.getName(),
                                               timer_level_fine)),
        update_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAA::update_") + target.getName() + "_" +
                                                     target.getName(),
                                                 timer_level_fine))
  {
    resize();
  }

  SoaDistanceTableAA()                          = delete;
  SoaDistanceTableAA(const SoaDistanceTableAA&) = delete;
  ~SoaDistanceTableAA() override {}

  size_t compute_size(int N) const
  {
    const size_t N_padded  = getAlignedSize<T>(N);
    const size_t Alignment = getAlignment<T>();
    return (N_padded * (2 * N - N_padded + 1) + (Alignment - 1) * N_padded) / 2;
  }

  void resize()
  {
    // initialize memory containers and views
    Ntargets_padded         = getAlignedSize<T>(N_targets);
    const size_t total_size = compute_size(N_targets);
    memory_pool_.resize(total_size * (1 + D));
    distances_.resize(N_targets);
    displacements_.resize(N_targets);
    for (int i = 0; i < N_targets; ++i)
    {
      distances_[i].attachReference(memory_pool_.data() + compute_size(i), i);
      displacements_[i].attachReference(i, total_size, memory_pool_.data() + total_size + compute_size(i));
    }

    old_r_.resize(N_targets);
    old_dr_.resize(N_targets);
    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy in the update function
    // temp_r_ is padded explicitly while temp_dr_ is padded internally
    temp_r_.resize(Ntargets_padded);
    temp_dr_.resize(N_targets);
  }

  const DistRow& getOldDists() const override { return old_r_; }
  const DisplRow& getOldDispls() const override { return old_dr_; }

  inline void evaluate(ParticleSet& P) override
  {
    ScopedTimer local_timer(evaluate_timer_);
    constexpr T BigR = std::numeric_limits<T>::max();
    for (int iat = 1; iat < N_targets; ++iat)
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), distances_[iat].data(),
                                             displacements_[iat], 0, iat, iat);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);

    old_prepared_elec_id = prepare_old ? iat : -1;

    DTD_BConds<T, D, SC>::computeDistances(rnew, P.getCoordinates().getAllParticlePos(), temp_r_.data(), temp_dr_, 0,
                                           N_targets, iat);
    // set up old_r_ and old_dr_ for moves may get accepted.
    if (prepare_old)
    {
      //recompute from scratch
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), old_r_.data(), old_dr_,
                                             0, N_targets, iat);
      old_r_[iat] = std::numeric_limits<T>::max(); //assign a big number
    }
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override
  {
    //ensure there are neighbors
    assert(N_targets > 1);
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < N_targets; ++jat)
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
      for (int jat = iat + 1; jat < N_targets; ++jat)
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
   * Since the upper triangle is not needed in the later computation,
   * only the [0,iat-1) columns need to save the new values.
   * The memory copy goes up to the padded size only for better performance.
   */
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    //update by a cache line
    const int nupdate = getAlignedSize<T>(iat);
    //copy row
    std::copy_n(temp_r_.data(), nupdate, distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(temp_dr_.data(idim), nupdate, displacements_[iat].data(idim));
    //copy column
    for (size_t i = iat + 1; i < N_targets; ++i)
    {
      distances_[i][iat]     = temp_r_[i];
      displacements_[i](iat) = -temp_dr_[i];
    }
  }

  void updatePartial(IndexType jat, bool from_temp) override
  {
    ScopedTimer local_timer(update_timer_);
    //update by a cache line
    const int nupdate = getAlignedSize<T>(jat);
    if (from_temp)
    {
      //copy row
      std::copy_n(temp_r_.data(), nupdate, distances_[jat].data());
      for (int idim = 0; idim < D; ++idim)
        std::copy_n(temp_dr_.data(idim), nupdate, displacements_[jat].data(idim));
    }
    else
    {
      assert(old_prepared_elec_id == jat);
      //copy row
      std::copy_n(old_r_.data(), nupdate, distances_[jat].data());
      for (int idim = 0; idim < D; ++idim)
        std::copy_n(old_dr_.data(idim), nupdate, displacements_[jat].data(idim));
    }
  }

private:
  /// timer for evaluate()
  NewTimer& evaluate_timer_;
  /// timer for move()
  NewTimer& move_timer_;
  /// timer for update()
  NewTimer& update_timer_;
};
} // namespace qmcplusplus
#endif
