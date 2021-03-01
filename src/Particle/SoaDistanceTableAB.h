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

#include "Utilities/FairDivide.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a transposed form
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAB : public DTD_BConds<T, D, SC>, public DistanceTableData
{
  SoaDistanceTableAB(const ParticleSet& source, ParticleSet& target)
      : DTD_BConds<T, D, SC>(source.Lattice),
        DistanceTableData(source, target),
        evaluate_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAB::evaluate_") + target.getName() +
                                                       "_" + source.getName(),
                                                   timer_level_fine)),
        move_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAB::move_") + target.getName() + "_" +
                                                   source.getName(),
                                               timer_level_fine)),
        update_timer_(*timer_manager.createTimer(std::string("SoaDistanceTableAB::update_") + target.getName() + "_" +
                                                     source.getName(),
                                                 timer_level_fine))
  {
    resize(source.getTotalNum(), target.getTotalNum());
  }

  void resize(int ns, int nt)
  {
    N_sources = ns;
    N_targets = nt;
    if (N_sources * N_targets == 0)
      return;

    // initialize memory containers and views
    const int Nsources_padded = getAlignedSize<T>(N_sources);
    distances_.resize(N_targets);
    displacements_.resize(N_targets);
    for (int i = 0; i < N_targets; ++i)
    {
      distances_[i].resize(Nsources_padded);
      displacements_[i].resize(Nsources_padded);
    }

    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy in the update function
    // temp_r_ is padded explicitly while temp_dr_ is padded internally
    temp_r_.resize(Nsources_padded);
    temp_dr_.resize(N_sources);
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
      FairDivideAligned(N_sources, getAlignment<T>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

      //be aware of the sign of Displacement
      for (int iat = 0; iat < N_targets; ++iat)
        DTD_BConds<T, D, SC>::computeDistances(P.R[iat], Origin->getCoordinates().getAllParticlePos(),
                                               distances_[iat].data(), displacements_[iat], first, last);
    }
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) override
  {
    ScopedTimer local_timer(move_timer_);
    DTD_BConds<T, D, SC>::computeDistances(rnew, Origin->getCoordinates().getAllParticlePos(), temp_r_.data(), temp_dr_,
                                           0, N_sources);
    // If the full table is not ready all the time, overwrite the current value.
    // If this step is missing, DT values can be undefined in case a move is rejected.
    if (!need_full_table_ && prepare_old)
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], Origin->getCoordinates().getAllParticlePos(),
                                             distances_[iat].data(), displacements_[iat], 0, N_sources);
  }

  ///update the stripe for jat-th particle
  inline void update(IndexType iat) override
  {
    ScopedTimer local_timer(update_timer_);
    std::copy_n(temp_r_.data(), N_sources, distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(temp_dr_.data(idim), N_sources, displacements_[iat].data(idim));
  }

  size_t get_neighbors(int iat,
                       RealType rcut,
                       int* restrict jid,
                       RealType* restrict dist,
                       PosType* restrict displ) const override
  {
    constexpr T cminus(-1);
    size_t nn = 0;
    for (int jat = 0; jat < N_targets; ++jat)
    {
      const RealType rij = distances_[jat][iat];
      if (rij < rcut)
      { //make the compact list
        jid[nn]   = jat;
        dist[nn]  = rij;
        displ[nn] = cminus * displacements_[jat][iat];
        nn++;
      }
    }
    return nn;
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const override
  {
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < N_sources; ++jat)
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
      for (int jat = 0; jat < N_sources; ++jat)
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
    return index;
  }

  size_t get_neighbors(int iat, RealType rcut, RealType* restrict dist) const
  {
    size_t nn = 0;
    for (int jat = 0; jat < N_targets; ++jat)
    {
      const RealType rij = distances_[jat][iat];
      if (rij < rcut)
      { //make the compact list
        dist[nn] = rij;
        nn++;
      }
    }
    return nn;
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
