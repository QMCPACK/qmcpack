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
#ifndef QMCPLUSPLUS_DTDIMPL_AB_H
#define QMCPLUSPLUS_DTDIMPL_AB_H

#include "Lattice/ParticleBConds3DSoa.h"
#include "Utilities/FairDivide.h"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a transposed form
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAB : public DTD_BConds<T, D, SC>, public DistanceTableAB
{
  SoaDistanceTableAB(const ParticleSet& source, ParticleSet& target)
      : DTD_BConds<T, D, SC>(source.getLattice()),
        DistanceTableAB(source, target, DTModes::NEED_TEMP_DATA_ON_HOST),
        evaluate_timer_(
            *timer_manager.createTimer(std::string("DTAB::evaluate_") + target.getName() + "_" + source.getName(),
                                       timer_level_fine)),
        move_timer_(*timer_manager.createTimer(std::string("DTAB::move_") + target.getName() + "_" + source.getName(),
                                               timer_level_fine)),
        update_timer_(
            *timer_manager.createTimer(std::string("DTAB::update_") + target.getName() + "_" + source.getName(),
                                       timer_level_fine))
  {
    resize();
  }

  void resize()
  {
    if (num_sources_ * num_targets_ == 0)
      return;

    // initialize memory containers and views
    const int num_sources_padded = getAlignedSize<T>(num_sources_);
    distances_.resize(num_targets_);
    displacements_.resize(num_targets_);
    for (int i = 0; i < num_targets_; ++i)
    {
      distances_[i].resize(num_sources_padded);
      displacements_[i].resize(num_sources_padded);
    }

    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy in the update function
    // temp_r_ is padded explicitly while temp_dr_ is padded internally
    temp_r_.resize(num_sources_padded);
    temp_dr_.resize(num_sources_);
  }

  SoaDistanceTableAB()                          = delete;
  SoaDistanceTableAB(const SoaDistanceTableAB&) = delete;

  /** evaluate the full table */
  inline void evaluate(ParticleSet& P) override
  {
    ScopedTimer local_timer(evaluate_timer_);
#pragma omp parallel
    {
      int first, last;
      FairDivideAligned(num_sources_, getAlignment<T>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

      //be aware of the sign of Displacement
      for (int iat = 0; iat < num_targets_; ++iat)
        DTD_BConds<T, D, SC>::computeDistances(P.R[iat], origin_.getCoordinates().getAllParticlePos(),
                                               distances_[iat].data(), displacements_[iat], first, last);
    }
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
  /// timer for evaluate()
  NewTimer& evaluate_timer_;
  /// timer for move()
  NewTimer& move_timer_;
  /// timer for update()
  NewTimer& update_timer_;
};
} // namespace qmcplusplus
#endif
