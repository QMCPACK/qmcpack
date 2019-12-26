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
#include "simd/algorithm.hpp"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAA : public DTD_BConds<T, D, SC>, public DistanceTableData
{
  int Ntargets_padded;

  ///actual memory for Displacements
  aligned_vector<RealType> memoryPool_displs_;

  SoaDistanceTableAA(ParticleSet& target) : DTD_BConds<T, D, SC>(target.Lattice), DistanceTableData(target, target)
  {
    resize(target.getTotalNum());
  }

  SoaDistanceTableAA()                          = delete;
  SoaDistanceTableAA(const SoaDistanceTableAA&) = delete;
  ~SoaDistanceTableAA() {}

  size_t compute_size(int N)
  {
    const size_t N_padded  = getAlignedSize<T>(N);
    const size_t Alignment = getAlignment<T>();
    return (N_padded * (2 * N - N_padded + 1) + (Alignment - 1) * N_padded) / 2;
  }

  void resize(int n)
  {
    N_sources = N_targets = n;

    // initialize memory containers and views
    Ntargets_padded                             = getAlignedSize<T>(n);
    const size_t total_size = compute_size(N_targets);
    memoryPool_displs_.resize(total_size * D);
    Distances.resize(N_targets);
    Displacements.resize(N_targets);
    for (int i = 0; i < N_targets; ++i)
    {
      Distances[i].resize(Ntargets_padded);
      Displacements[i].attachReference(i, total_size, memoryPool_displs_.data() + compute_size(i));
    }

    // The padding of Temp_r and Temp_dr is necessary for the memory copy in the update function
    // Temp_r is padded explicitly while Temp_dr is padded internally
    Temp_r.resize(Ntargets_padded);
    Temp_dr.resize(N_targets);
  }

  inline void evaluate(ParticleSet& P)
  {
    constexpr T BigR = std::numeric_limits<T>::max();
    //P.RSoA.copyIn(P.R);
    for (int iat = 0; iat < N_targets; ++iat)
    {
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], P.RSoA, Distances[iat].data(), Displacements[iat], 0, N_targets, iat);
      Distances[iat][iat] = BigR; //assign big distance
    }
  }

  inline void evaluate(ParticleSet& P, IndexType jat)
  {
    DTD_BConds<T, D, SC>::computeDistances(P.R[jat], P.RSoA, Distances[jat].data(), Displacements[jat], 0, N_targets, jat);
    Distances[jat][jat] = std::numeric_limits<T>::max(); //assign a big number
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew)
  {
    DTD_BConds<T, D, SC>::computeDistances(rnew, P.RSoA, Temp_r.data(), Temp_dr, 0, N_targets, P.activePtcl);
  }

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
  {
    RealType min_dist = std::numeric_limits<RealType>::max();
    int index         = -1;
    if (newpos)
    {
      for (int jat = 0; jat < N_targets; ++jat)
        if (Temp_r[jat] < min_dist && jat != iat)
        {
          min_dist = Temp_r[jat];
          index    = jat;
        }
      if (index >= 0)
        dr = Temp_dr[index];
    }
    else
    {
      for (int jat = 0; jat < N_targets; ++jat)
        if (Distances[iat][jat] < min_dist && jat != iat)
        {
          min_dist = Distances[iat][jat];
          index    = jat;
        }
      if (index >= 0)
        dr = Displacements[iat][index];
    }
    r = min_dist;
    return index;
  }

  /** After accepting the iat-th particle, update the iat-th row of Distances and Displacements.
   * Since the upper triangle is not needed in the later computation,
   * only the [0,iat-1) columns need to save the new values.
   * The memory copy goes up to the padded size only for better performance.
   */
  inline void update(IndexType iat)
  {
    if (iat == 0)
      return;
    //update by a cache line
    const int nupdate = getAlignedSize<T>(iat);
    std::copy_n(Temp_r.data(), nupdate, Distances[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(Temp_dr.data(idim), nupdate, Displacements[iat].data(idim));
  }

  inline void storeCurrent(IndexType iat)
  {
  }
};
} // namespace qmcplusplus
#endif
