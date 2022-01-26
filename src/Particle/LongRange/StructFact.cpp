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
StructFact::StructFact(int nptcls, int ns, const ParticleLayout& lattice, const KContainer& k_lists)
    : DoUpdate(false),
      SuperCellEnum(SUPERCELL_BULK),
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

  resize(k_lists_.numk);
}

//Destructor
StructFact::~StructFact() = default;

void StructFact::resize(int nkpts)
{
  phiV.resize(nkpts);
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r.resize(num_species, nkpts);
  rhok_i.resize(num_species, nkpts);
  if (StorePerParticle)
  {
    eikr_r.resize(num_ptcls, nkpts);
    eikr_i.resize(num_ptcls, nkpts);
  }
  eikr_r_temp.resize(nkpts);
  eikr_i_temp.resize(nkpts);
#else
  rhok.resize(num_species, nkpts);
  eikr.resize(num_ptcls, nkpts);
  eikr_temp.resize(nkpts);
#endif
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
  int npart = P.getTotalNum();
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r = 0.0;
  rhok_i = 0.0;
  //algorithmA
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
      for (int ki = 0; ki < nk; ki++)
        phiV[ki] = dot(k_lists_.kpts_cart[ki], pos);
      eval_e2iphi(nk, phiV.data(), eikr_r_temp.data(), eikr_i_temp.data());
      for (int ki = 0; ki < nk; ki++)
      {
        rhok_r_ptr[ki] += eikr_r_temp[ki];
        rhok_i_ptr[ki] += eikr_i_temp[ki];
      }
#endif
    }
  }
#else
  rhok = 0.0;
  for (int i = 0; i < npart; i++)
  {
    PosType pos(P.R[i]);
    RealType s, c; //get sin and cos
    ComplexType* restrict eikr_ref = eikr[i];
    ComplexType* restrict rhok_ref = rhok[P.getGroupID(i)];
    for (int ki = 0; ki < k_lists_.numk; ki++)
    {
      qmcplusplus::sincos(dot(k_lists_.kpts_cart[ki], pos), &s, &c);
      eikr_ref[ki] = ComplexType(c, s);
      rhok_ref[ki] += eikr_ref[ki];
    }
  }
#endif
}


void StructFact::makeMove(int active, const PosType& pos)
{
#if defined(USE_REAL_STRUCT_FACTOR)
#pragma omp simd
  for (int ki = 0; ki < k_lists_.numk; ki++)
    qmcplusplus::sincos(dot(k_lists_.kpts_cart[ki], pos), &eikr_i_temp[ki], &eikr_r_temp[ki]);
#else
  RealType s, c; //get sin and cos
  for (int ki = 0; ki < k_lists_.numk; ++ki)
  {
    qmcplusplus::sincos(dot(k_lists_.kpts_cart[ki], pos), &s, &c);
    eikr_temp[ki] = ComplexType(c, s);
  }
#endif
}

void StructFact::acceptMove(int active, int gid, const PosType& rold)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  if (StorePerParticle)
  {
    RealType* restrict eikr_ptr_r = eikr_r[active];
    RealType* restrict eikr_ptr_i = eikr_i[active];
    RealType* restrict rhok_ptr_r(rhok_r[gid]);
    RealType* restrict rhok_ptr_i(rhok_i[gid]);
    for (int ki = 0; ki < k_lists_.numk; ++ki)
    {
      rhok_ptr_r[ki] += (eikr_r_temp[ki] - eikr_ptr_r[ki]);
      rhok_ptr_i[ki] += (eikr_i_temp[ki] - eikr_ptr_i[ki]);
      eikr_ptr_r[ki] = eikr_r_temp[ki];
      eikr_ptr_i[ki] = eikr_i_temp[ki];
    }
  }
  else
  {
    RealType* restrict rhok_ptr_r(rhok_r[gid]);
    RealType* restrict rhok_ptr_i(rhok_i[gid]);

// add the new value and subtract the old value
#pragma omp simd
    for (int ki = 0; ki < k_lists_.numk; ++ki)
    {
      RealType s, c;
      qmcplusplus::sincos(dot(k_lists_.kpts_cart[ki], rold), &s, &c);
      rhok_ptr_r[ki] += eikr_r_temp[ki] - c;
      rhok_ptr_i[ki] += eikr_i_temp[ki] - s;
    }
  }
#else
  ComplexType* restrict eikr_ptr = eikr[active];
  ComplexType* restrict rhok_ptr(rhok[gid]);
  for (int ki = 0; ki < k_lists_.numk; ++ki)
  {
    rhok_ptr[ki] += (eikr_temp[ki] - eikr_ptr[ki]);
    eikr_ptr[ki] = eikr_temp[ki];
  }
#endif
}

void StructFact::rejectMove(int active, int gid)
{
  //do nothing
}

void StructFact::turnOnStorePerParticle(const ParticleSet& P)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  if (!StorePerParticle)
  {
    StorePerParticle = true;
    const int nptcl  = P.getTotalNum();
    eikr_r.resize(nptcl, k_lists_.numk);
    eikr_i.resize(nptcl, k_lists_.numk);
    computeRhok(P);
  }
#endif
}

} // namespace qmcplusplus
