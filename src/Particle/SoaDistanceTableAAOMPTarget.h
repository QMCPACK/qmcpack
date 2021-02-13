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

#include "CPU/SIMD/algorithm.hpp"
#include "OMPTarget/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Particle/RealSpacePositionsOMPTarget.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAAOMPTarget : public DTD_BConds<T, D, SC>, public DistanceTableData
{
  ///number of targets with padding
  int Ntargets_padded;

  ///actual memory for displacements_
  aligned_vector<RealType> memory_pool_displs_;

  /// old distances
  DistRow old_r_mem_;
  DistRow old_r_;

  /// old displacements
  DisplRow old_dr_mem_;
  DisplRow old_dr_;

  DistRow temp_r_mem_;
  DisplRow temp_dr_mem_;

  SoaDistanceTableAAOMPTarget(ParticleSet& target)
      : DTD_BConds<T, D, SC>(target.Lattice), DistanceTableData(target, target)
  {
    auto* coordinates_soa = dynamic_cast<const RealSpacePositionsOMPTarget*>(&target.getCoordinates());
    if (!coordinates_soa)
      throw std::runtime_error("Source particle set doesn't have OpenMP offload. Contact developers!");
    resize(target.getTotalNum());
    PRAGMA_OFFLOAD("omp target enter data map(to : this[:1])")
  }

  SoaDistanceTableAAOMPTarget()                                   = delete;
  SoaDistanceTableAAOMPTarget(const SoaDistanceTableAAOMPTarget&) = delete;
  ~SoaDistanceTableAAOMPTarget(){PRAGMA_OFFLOAD("omp target exit data map(delete : this[:1])")}

  size_t compute_size(int N) const
  {
    const size_t N_padded  = getAlignedSize<T>(N);
    const size_t Alignment = getAlignment<T>();
    return (N_padded * (2 * N - N_padded + 1) + (Alignment - 1) * N_padded) / 2;
  }

  void resize(int n)
  {
    N_sources = N_targets = n;

    // initialize memory containers and views
    Ntargets_padded         = getAlignedSize<T>(n);
    const size_t total_size = compute_size(N_targets);
    memory_pool_displs_.resize(total_size * D);
    distances_.resize(N_targets);
    displacements_.resize(N_targets);
    for (int i = 0; i < N_targets; ++i)
    {
      distances_[i].resize(Ntargets_padded);
      displacements_[i].attachReference(i, total_size, memory_pool_displs_.data() + compute_size(i));
    }

    old_r_mem_.resize(N_targets);
    old_dr_mem_.resize(N_targets);
    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy in the update function
    // temp_r_ is padded explicitly while temp_dr_ is padded internally
    temp_r_mem_.resize(Ntargets_padded);
    temp_dr_mem_.resize(N_targets);
  }

  const DistRow& getOldDists() const { return old_r_; }
  const DisplRow& getOldDispls() const { return old_dr_; }

  inline void evaluate(ParticleSet& P)
  {
    constexpr T BigR = std::numeric_limits<T>::max();
    for (int iat = 0; iat < N_targets; ++iat)
    {
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), distances_[iat].data(),
                                             displacements_[iat], 0, N_targets, iat);
      distances_[iat][iat] = BigR; //assign big distance
    }
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old)
  {
    temp_r_.attachReference(temp_r_mem_.data(), temp_r_mem_.size());
    temp_dr_.attachReference(temp_dr_mem_.size(), temp_dr_mem_.capacity(), temp_dr_mem_.data());

    DTD_BConds<T, D, SC>::computeDistances(rnew, P.getCoordinates().getAllParticlePos(), temp_r_.data(), temp_dr_, 0,
                                           N_targets, P.activePtcl);
    // set up old_r_ and old_dr_ for moves may get accepted.
    if (prepare_old)
    {
      old_r_.attachReference(old_r_mem_.data(), old_r_mem_.size());
      old_dr_.attachReference(old_dr_mem_.size(), old_dr_mem_.capacity(), old_dr_mem_.data());
      //recompute from scratch
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.getCoordinates().getAllParticlePos(), old_r_.data(), old_dr_,
                                             0, N_targets, iat);
      old_r_[iat] = std::numeric_limits<T>::max(); //assign a big number

      // If the full table is not ready all the time, overwrite the current value.
      // If this step is missing, DT values can be undefined in case a move is rejected.
      if (!need_full_table_)
      {
        //copy row
        std::copy_n(old_r_.data(), iat, distances_[iat].data());
        for (int idim = 0; idim < D; ++idim)
          std::copy_n(old_dr_.data(idim), iat, displacements_[iat].data(idim));
      }
    }
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
  {
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
      if (index >= 0)
        dr = temp_dr_[index];
    }
    else
    {
      for (int jat = 0; jat < N_targets; ++jat)
        if (distances_[iat][jat] < min_dist && jat != iat)
        {
          min_dist = distances_[iat][jat];
          index    = jat;
        }
      if (index >= 0)
        dr = displacements_[iat][index];
    }
    r = min_dist;
    return index;
  }

  /** After accepting the iat-th particle, update the iat-th row of distances_ and displacements_.
   * Since the upper triangle is not needed in the later computation,
   * only the [0,iat-1) columns need to save the new values.
   * The memory copy goes up to the padded size only for better performance.
   */
  inline void update(IndexType iat, bool partial_update)
  {
    //update by a cache line
    const int nupdate = getAlignedSize<T>(iat);
    //copy row
    std::copy_n(temp_r_.data(), nupdate, distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(temp_dr_.data(idim), nupdate, displacements_[iat].data(idim));
    // This is an optimization to reduce update >iat rows during p-by-p forward move when no consumer needs full table.
    if (need_full_table_ || !partial_update)
    {
      //copy column
      for (size_t i = iat + 1; i < N_targets; ++i)
      {
        distances_[i][iat]     = temp_r_[i];
        displacements_[i](iat) = -temp_dr_[i];
      }
    }
  }
};
} // namespace qmcplusplus
#endif
