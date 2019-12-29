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
#ifndef QMCPLUSPLUS_DTDIMPL_BA_H
#define QMCPLUSPLUS_DTDIMPL_BA_H

#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a transposed form
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableBA : public DTD_BConds<T, D, SC>, public DistanceTableData
{
  SoaDistanceTableBA(const ParticleSet& source, ParticleSet& target)
      : DTD_BConds<T, D, SC>(source.Lattice), DistanceTableData(source, target)
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
    Distances.resize(N_targets);
    Displacements.resize(N_targets);
    for (int i = 0; i < N_targets; ++i)
    {
      Distances[i].resize(Nsources_padded);
      Displacements[i].resize(Nsources_padded);
    }

    // The padding of Temp_r and Temp_dr is necessary for the memory copy in the update function
    // Temp_r is padded explicitly while Temp_dr is padded internally
    Temp_r.resize(Nsources_padded);
    Temp_dr.resize(N_sources);
  }

  SoaDistanceTableBA()                          = delete;
  SoaDistanceTableBA(const SoaDistanceTableBA&) = delete;
  ~SoaDistanceTableBA() {}

  /** evaluate the full table */
  inline void evaluate(ParticleSet& P)
  {
#pragma omp parallel
    {
      int first, last;
      FairDivideAligned(N_sources, getAlignment<T>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

      //be aware of the sign of Displacement
      for (int iat = 0; iat < N_targets; ++iat)
        DTD_BConds<T, D, SC>::computeDistances(P.R[iat], Origin->RSoA, Distances[iat].data(), Displacements[iat], first, last);
    }
  }

  /** evaluate the iat-row with the current position
   *
   * Fill Temp_r and Temp_dr and copy them Distances & Displacements
   */
  inline void evaluate(ParticleSet& P, IndexType iat)
  {
    //DTD_BConds<T, D, SC>::computeDistances(P.R[iat], Origin->RSoA, Distances[iat].data(), Displacements[iat], 0, N_sources);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old)
  {
    DTD_BConds<T, D, SC>::computeDistances(rnew, Origin->RSoA, Temp_r.data(), Temp_dr, 0, N_sources);
  }

  ///update the stripe for jat-th particle
  inline void update(IndexType iat)
  {
    std::copy_n(Temp_r.data(), N_sources, Distances[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(Temp_dr.data(idim), N_sources, Displacements[iat].data(idim));
  }

  size_t get_neighbors(int iat,
                       RealType rcut,
                       int* restrict jid,
                       RealType* restrict dist,
                       PosType* restrict displ) const
  {
    constexpr T cminus(-1);
    size_t nn = 0;
    for (int jat = 0; jat < N_targets; ++jat)
    {
      const RealType rij = Distances[jat][iat];
      if (rij < rcut)
      { //make the compact list
        jid[nn]   = jat;
        dist[nn]  = rij;
        displ[nn] = cminus * Displacements[jat][iat];
        nn++;
      }
    }
    return nn;
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
  {
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < N_sources; ++jat)
        if (Temp_r[jat] < min_dist)
        {
          min_dist = Temp_r[jat];
          index    = jat;
        }
      if (index >= 0)
      {
        r  = min_dist;
        dr = Temp_dr[index];
      }
    }
    else
    {
      for (int jat = 0; jat < N_sources; ++jat)
        if (Distances[iat][jat] < min_dist)
        {
          min_dist = Distances[iat][jat];
          index    = jat;
        }
      if (index >= 0)
      {
        r  = min_dist;
        dr = Displacements[iat][index];
      }
    }
    return index;
  }

  size_t get_neighbors(int iat, RealType rcut, RealType* restrict dist) const
  {
    size_t nn = 0;
    for (int jat = 0; jat < N_targets; ++jat)
    {
      const RealType rij = Distances[jat][iat];
      if (rij < rcut)
      { //make the compact list
        dist[nn] = rij;
        nn++;
      }
    }
    return nn;
  }
};
} // namespace qmcplusplus
#endif
