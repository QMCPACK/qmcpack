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
#ifndef QMCPLUSPLUS_DTDIMPL_AB_OMP_H
#define QMCPLUSPLUS_DTDIMPL_AB_OMP_H

#include "OpenMP/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"

namespace qmcplusplus
{
/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for AB using a transposed form
 */
template<typename T, unsigned D, int SC>
class SoaDistanceTableABOMP : public DTD_BConds<T, D, SC>, public DistanceTableData
{
private:
  Vector<RealType, OMPallocator<RealType, PinnedAlignedAllocator<RealType>>> memoryPool;

public:
  SoaDistanceTableABOMP(const ParticleSet& source, ParticleSet& target)
      : DTD_BConds<T, D, SC>(source.Lattice), DistanceTableData(source, target)
  {
    resize(source.getTotalNum(), target.getTotalNum());
    #pragma omp target enter data map(to:this[:1])
  }

  void resize(int ns, int nt)
  {
    N_sources = ns;
    N_targets = nt;
    if (N_sources * N_targets == 0)
      return;

    // initialize memory containers and views
    const int N_sources_padded = getAlignedSize<T>(N_sources);
    // for both distances and displacements
    memoryPool.resize(N_targets * N_sources_padded * (D + 1));
    // first part is for distances
    const size_t head_offset = N_targets * N_sources_padded;

    distances_.resize(N_targets);
    displacements_.resize(N_targets);
    for (int i = 0; i < N_targets; ++i)
    {
      distances_[i].attachReference(memoryPool.data() + i * N_sources_padded, N_sources);
      displacements_[i].attachReference(N_sources, N_sources_padded,
                                        memoryPool.data() + head_offset + i * N_sources_padded * D);
    }

    // The padding of temp_r_ and temp_dr_ is necessary for the memory copy in the update function
    // temp_r_ is padded explicitly while temp_dr_ is padded internally
    temp_r_.resize(N_sources_padded);
    temp_dr_.resize(N_sources);
  }

  SoaDistanceTableABOMP()                             = delete;
  SoaDistanceTableABOMP(const SoaDistanceTableABOMP&) = delete;

  ~SoaDistanceTableABOMP()
  {
    #pragma omp target exit data map(delete:this[:1])
  }

  /** evaluate the full table */
  inline void evaluate(ParticleSet& P)
  {
    // be aware of the sign of Displacement
    int N_targets_local  = N_targets;
    int N_targets_padded = getAlignedSize<T>(N_targets);
    int N_sources_local  = N_sources;
    int N_sources_padded = getAlignedSize<T>(N_sources);

    auto* target_pos_ptr       = P.RSoA->getAllParticlePos().data();
    const auto* source_pos_ptr = Origin->RSoA->getAllParticlePos().data();
    auto* r_dr_ptr             = memoryPool.data();

    const int ChunkSizePerTeam = 256;
    const int NumTeams         = (N_sources + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    #pragma omp target teams distribute collapse(2) num_teams(N_targets*NumTeams) \
      map(to: source_pos_ptr[:N_sources_padded*D], target_pos_ptr[:N_targets_padded*D]) \
      map(always, from: r_dr_ptr[:memoryPool.size()])
    for (int iat = 0; iat < N_targets; ++iat)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > N_sources_local ? N_sources_local : first + ChunkSizePerTeam;

        T pos[D];
        for (int idim = 0; idim < D; idim++)
          pos[idim] = target_pos_ptr[idim * N_targets_padded + iat];

        auto* r_iat_ptr  = r_dr_ptr + N_sources_padded * iat;
        auto* dr_iat_ptr = r_dr_ptr + N_sources_padded * N_targets + N_sources_padded * D * iat;

        DTD_BConds<T, D, SC>::computeDistancesOffload(pos, source_pos_ptr, r_iat_ptr, dr_iat_ptr, N_sources_padded,
                                                      first, last);
      }

    /*
    #pragma omp critical
    {
      std::cout << "tid = " << omp_get_thread_num() << " source[0] " << Origin->R[0] << " SoA " << Origin->RSoA[0] << std::endl;
      for (int iat = 0; iat < N_targets; ++iat)
      {
        std::cout << iat << " : " << P.R[iat];
        for (int isource = 0; isource < N_sources; ++isource)
          std::cout << " " << Distances[iat][isource];
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    abort();
    */
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old)
  {
    DTD_BConds<T, D, SC>::computeDistances(rnew, Origin->RSoA->getAllParticlePos(), temp_r_.data(), temp_dr_, 0,
                                           N_sources);
    // If the full table is not ready all the time, overwrite the current value.
    // If this step is missing, DT values can be undefined in case a move is rejected.
    if (!need_full_table_)
      DTD_BConds<T, D, SC>::computeDistances(P.R[iat], Origin->RSoA->getAllParticlePos(), distances_[iat].data(),
                                             displacements_[iat], 0, N_sources);
  }

  ///update the stripe for jat-th particle
  inline void update(IndexType iat, bool partial_update)
  {
    std::copy_n(temp_r_.data(), N_sources, distances_[iat].data());
    for (int idim = 0; idim < D; ++idim)
      std::copy_n(temp_dr_.data(idim), N_sources, displacements_[iat].data(idim));
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

  int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
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
};
} // namespace qmcplusplus
#endif
