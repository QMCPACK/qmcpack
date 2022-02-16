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
#include "OMPTarget/OMPTargetMath.hpp"
#include "RealSpacePositionsOMPTarget.h"

namespace qmcplusplus
{
//Constructor - pass arguments to k_lists_' constructor
StructFact::StructFact(const ParticleLayout& lattice, const KContainer& k_lists)
    : SuperCellEnum(SUPERCELL_BULK),
      k_lists_(k_lists),
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

void StructFact::resize(int nkpts, int num_species, int num_ptcls)
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
                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                  SKMultiWalkerMem& mw_mem)
{
  auto& sk_leader = sk_list.getLeader();
  auto& p_leader  = p_list.getLeader();
  ScopedTimer local(sk_leader.update_all_timer_);
  if (p_leader.getCoordinates().getKind() != DynamicCoordinateKind::DC_POS_OFFLOAD || sk_leader.StorePerParticle)
    for (int iw = 0; iw < sk_list.size(); iw++)
      sk_list[iw].computeRhok(p_list[iw]);
  else
  {
    const size_t nw          = p_list.size();
    const size_t num_species = p_leader.groups();
    const auto& kpts_cart    = sk_leader.k_lists_.get_kpts_cart_soa();
    const size_t nk          = sk_leader.k_lists_.numk;
    const size_t nk_padded   = kpts_cart.capacity();

    auto& coordinates_leader = static_cast<const RealSpacePositionsOMPTarget&>(p_leader.getCoordinates());
    auto& mw_rsoa_dev_ptrs   = coordinates_leader.getMultiWalkerRSoADevicePtrs();
    const size_t np_padded   = p_leader.getCoordinates().getAllParticlePos().capacity();

    constexpr size_t cplx_stride = 2;
    mw_mem.nw_rhok.resize(nw * num_species * cplx_stride, nk_padded);

    // make the compute over nk by blocks
    constexpr size_t kblock_size = 512;
    const size_t num_kblocks     = (nk + kblock_size) / kblock_size;

    auto* mw_rsoa_ptr   = mw_rsoa_dev_ptrs.data();
    auto* kpts_cart_ptr = kpts_cart.data();
    auto* mw_rhok_ptr   = mw_mem.nw_rhok.data();
    auto* group_offsets = p_leader.get_group_offsets().data();

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) map(always, from : mw_rhok_ptr[:mw_mem.nw_rhok.size()])")
    for (int iw = 0; iw < nw; iw++)
      for (int ib = 0; ib < num_kblocks; ib++)
      {
        const size_t offset          = ib * kblock_size;
        const size_t this_block_size = omptarget::min(kblock_size, nk - offset);
        const auto* rsoa_ptr         = mw_rsoa_ptr[iw];

        PRAGMA_OFFLOAD("omp parallel for")
        for (int ik = 0; ik < this_block_size; ik++)
          for (int is = 0; is < num_species; is++)
          {
            RealType rhok_r(0), rhok_i(0);

            for (int ip = group_offsets[is]; ip < group_offsets[is + 1]; ip++)
            {
              RealType s, c, phase(0);
              for (int idim = 0; idim < DIM; idim++)
                phase += kpts_cart_ptr[ik + offset + nk_padded * idim] * rsoa_ptr[ip + idim * np_padded];
              omptarget::sincos(phase, &s, &c);
              rhok_r += c;
              rhok_i += s;
            }

            mw_rhok_ptr[(iw * num_species + is) * cplx_stride * nk_padded + offset + ik]             = rhok_r;
            mw_rhok_ptr[(iw * num_species + is) * cplx_stride * nk_padded + nk_padded + offset + ik] = rhok_i;
          }
      }

    for (int iw = 0; iw < nw; iw++)
      for (int is = 0; is < num_species; is++)
      {
        std::copy_n(mw_mem.nw_rhok[(iw * num_species + is) * cplx_stride], nk, sk_list[iw].rhok_r[is]);
        std::copy_n(mw_mem.nw_rhok[(iw * num_species + is) * cplx_stride + 1], nk, sk_list[iw].rhok_i[is]);
      }
  }
}


/** evaluate rok per species, eikr  per particle
 */
void StructFact::computeRhok(const ParticleSet& P)
{
  const size_t num_ptcls   = P.getTotalNum();
  const size_t num_species = P.groups();
  const size_t nk          = k_lists_.numk;
  resize(nk, num_species, num_ptcls);

  rhok_r = 0.0;
  rhok_i = 0.0;
  if (StorePerParticle)
  {
    // save per particle and species value
    for (int i = 0; i < num_ptcls; ++i)
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
    for (int i = 0; i < num_ptcls; ++i)
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
      const size_t num_kblocks     = (nk + kblock_size) / kblock_size;
      RealType phiV[kblock_size], eikr_r_temp[kblock_size], eikr_i_temp[kblock_size];

      for (int ib = 0; ib < num_kblocks; ib++)
      {
        const size_t offset          = ib * kblock_size;
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
