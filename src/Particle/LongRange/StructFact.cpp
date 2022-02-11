//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "StructFact.h"
#include "CPU/math.hpp"
#include "CPU/e2iphi.h"
#include "CPU/SIMD/vmath.hpp"
#include "CPU/BLAS.hpp"
#include "Utilities/qmc_common.h"

namespace qmcplusplus
{
//Constructor - pass arguments to k_lists_' constructor
StructFact::StructFact(int ns, int nptcls, const ParticleLayout& lattice, const KContainer& k_lists)
    : SuperCellEnum(SUPERCELL_BULK),
      k_lists_(k_lists),
      num_ptcls(nptcls),
      num_species(ns),
      StorePerParticle(false),
      update_all_timer_(*timer_manager.createTimer("StructFact::update_all_part", timer_level_fine))
{
  if (qmc_common.use_ewald && lattice.SuperCellEnum == SUPERCELL_SLAB)
  {
    app_log() << "  Setting StructFact::SuperCellEnum=SUPERCELL_SLAB " << std::endl;
    SuperCellEnum = SUPERCELL_SLAB;
  }
}

//Destructor
StructFact::~StructFact() = default;

void StructFact::resize(int nkpts)
{
  rhok_r.resize(num_species, nkpts);
  rhok_i.resize(num_species, nkpts);
  if (StorePerParticle)
  {
    eikr_r.resize(num_ptcls, nkpts);
    eikr_i.resize(num_ptcls, nkpts);
  }
}


void StructFact::updateAllPart(const ParticleSet& P)
{
  ScopedTimer local(update_all_timer_);
  computeRhok(P);
}

void StructFact::mw_updateAllPart(const RefVectorWithLeader<StructFact>& sk_list,
                                  const RefVectorWithLeader<ParticleSet>& p_list)
{
  auto& sk_leader = sk_list.getLeader();
  ScopedTimer local(sk_leader.update_all_timer_);
  for (int iw = 0; iw < sk_list.size(); iw++)
    sk_list[iw].computeRhok(p_list[iw]);
}


/** evaluate rok per species, eikr  per particle
 */
void StructFact::computeRhok(const ParticleSet& P)
{
  resize(k_lists_.numk);

  rhok_r = 0.0;
  rhok_i = 0.0;
  const int npart = P.getTotalNum();
  const int nk = k_lists_.numk;
  if (StorePerParticle)
  {
    // save per particle and species value
    for (int i = 0; i < npart; ++i)
    {
      const auto& pos           = P.R[i];
      auto* restrict eikr_r_ptr = eikr_r[i];
      auto* restrict eikr_i_ptr = eikr_i[i];
      auto* restrict rhok_r_ptr = rhok_r[P.getGroupID(i)];
      auto* restrict rhok_i_ptr = rhok_i[P.getGroupID(i)];
#pragma omp simd
      for (int ki = 0; ki < nk; ki++)
      {
        qmcplusplus::sincos(dot(k_lists_.kpts_cart[ki], pos), &eikr_i_ptr[ki], &eikr_r_ptr[ki]);
        rhok_r_ptr[ki] += eikr_r_ptr[ki];
        rhok_i_ptr[ki] += eikr_i_ptr[ki];
      }
    }
  }
  else
  {
    // save per species value
    for (int i = 0; i < npart; ++i)
    {
      const auto& pos           = P.R[i];
      auto* restrict rhok_r_ptr = rhok_r[P.getGroupID(i)];
      auto* restrict rhok_i_ptr = rhok_i[P.getGroupID(i)];
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#pragma omp simd
      for (int ki = 0; ki < nk; ki++)
      {
        RealType s, c;
        qmcplusplus::sincos(dot(k_lists_.kpts_cart[ki], pos), &s, &c);
        rhok_r_ptr[ki] += c;
        rhok_i_ptr[ki] += s;
      }
#else
      // make the compute over nk by blocks
      constexpr size_t kblock_size = 512;
      const size_t num_kblocks = (nk + kblock_size) / kblock_size;
      RealType phiV[kblock_size], eikr_r_temp[kblock_size], eikr_i_temp[kblock_size];

      for(int ib = 0; ib < num_kblocks; ib++)
      {
        const size_t offset = ib * kblock_size;
        const size_t this_block_size = std::min(kblock_size, nk - offset);
        for (int ki = 0; ki < this_block_size; ki++)
          phiV[ki] = dot(k_lists_.kpts_cart[ki + offset], pos);
        eval_e2iphi(this_block_size, phiV, eikr_r_temp, eikr_i_temp);
        for (int ki = 0; ki < this_block_size; ki++)
        {
          rhok_r_ptr[ki + offset] += eikr_r_temp[ki];
          rhok_i_ptr[ki + offset] += eikr_i_temp[ki];
        }
      }
#endif
    }
  }
}

void StructFact::turnOnStorePerParticle(const ParticleSet& P)
{
  if (!StorePerParticle)
  {
    StorePerParticle = true;
    computeRhok(P);
  }
}

} // namespace qmcplusplus
