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


#include <LongRange/StructFact.h>
#include <config/stdlib/math.hpp>
#include <Numerics/e2iphi.h>
#include <simd/vmath.hpp>
#include <Numerics/OhmmsBlas.h>
#include <qmc_common.h>

namespace qmcplusplus
{
//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& P, RealType kc)
    : DoUpdate(false), StorePerParticle(false), SuperCellEnum(SUPERCELL_BULK)
{
  if (qmc_common.use_ewald && P.LRBox.SuperCellEnum == SUPERCELL_SLAB)
  {
    app_log() << "  Setting StructFact::SuperCellEnum=SUPERCELL_SLAB " << std::endl;
    SuperCellEnum = SUPERCELL_SLAB;
  }

  UpdateNewCell(P, kc);
}

//Destructor
StructFact::~StructFact() {}

void StructFact::UpdateNewCell(ParticleSet& P, RealType kc)
{
  //Generate the lists of k-vectors
  KLists.UpdateKLists(P.LRBox, kc);
  //resize any array
  resize(P.getSpeciesSet().size(), P.getTotalNum(), KLists.numk);
  //Compute the entire Rhok
  FillRhok(P);
}

void StructFact::resize(int ns, int nptcl, int nkpts)
{
  phiV.resize(nkpts);
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r.resize(ns, nkpts);
  rhok_i.resize(ns, nkpts);
  if (StorePerParticle)
  {
    eikr_r.resize(nptcl, nkpts);
    eikr_i.resize(nptcl, nkpts);
  }
  eikr_r_temp.resize(nkpts);
  eikr_i_temp.resize(nkpts);
#else
  rhok.resize(ns, nkpts);
  eikr.resize(nptcl, nkpts);
  eikr_temp.resize(nkpts);
#endif
}


void StructFact::UpdateAllPart(ParticleSet& P) { FillRhok(P); }


/** evaluate rok per species, eikr  per particle
 */
void StructFact::FillRhok(ParticleSet& P)
{
  int npart = P.getTotalNum();
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r = 0.0;
  rhok_i = 0.0;
  //algorithmA
  const int nk = KLists.numk;
  if (StorePerParticle)
  {
    // save per particle and species value
    for (int i = 0; i < npart; ++i)
    {
      const auto& pos           = P.R[i];
      auto* restrict eikr_r_ptr = eikr_r[i];
      auto* restrict eikr_i_ptr = eikr_i[i];
      auto* restrict rhok_r_ptr = rhok_r[P.GroupID[i]];
      auto* restrict rhok_i_ptr = rhok_i[P.GroupID[i]];
#pragma omp simd
      for (int ki = 0; ki < nk; ki++)
      {
        sincos(dot(KLists.kpts_cart[ki], pos), &eikr_i_ptr[ki], &eikr_r_ptr[ki]);
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
      auto* restrict rhok_r_ptr = rhok_r[P.GroupID[i]];
      auto* restrict rhok_i_ptr = rhok_i[P.GroupID[i]];
#ifdef __INTEL_COMPILER
#pragma omp simd
      for (int ki = 0; ki < nk; ki++)
      {
        RealType s, c;
        sincos(dot(KLists.kpts_cart[ki], pos), &s, &c);
        rhok_r_ptr[ki] += c;
        rhok_i_ptr[ki] += s;
      }
#else
      for (int ki = 0; ki < nk; ki++)
        phiV[ki] = dot(KLists.kpts_cart[ki], pos);
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
    ComplexType* restrict rhok_ref = rhok[P.GroupID[i]];
    for (int ki = 0; ki < KLists.numk; ki++)
    {
      sincos(dot(KLists.kpts_cart[ki], pos), &s, &c);
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
  for (int ki = 0; ki < KLists.numk; ki++)
    sincos(dot(KLists.kpts_cart[ki], pos), &eikr_i_temp[ki], &eikr_r_temp[ki]);
#else
  RealType s, c; //get sin and cos
  for (int ki = 0; ki < KLists.numk; ++ki)
  {
    sincos(dot(KLists.kpts_cart[ki], pos), &s, &c);
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
    for (int ki = 0; ki < KLists.numk; ++ki)
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
    for (int ki = 0; ki < KLists.numk; ++ki)
    {
      RealType s, c;
      sincos(dot(KLists.kpts_cart[ki], rold), &s, &c);
      rhok_ptr_r[ki] += eikr_r_temp[ki] - c;
      rhok_ptr_i[ki] += eikr_i_temp[ki] - s;
    }
  }
#else
  ComplexType* restrict eikr_ptr = eikr[active];
  ComplexType* restrict rhok_ptr(rhok[gid]);
  for (int ki = 0; ki < KLists.numk; ++ki)
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

void StructFact::turnOnStorePerParticle(ParticleSet& P)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  if (!StorePerParticle)
  {
    StorePerParticle = true;
    const int nptcl  = P.getTotalNum();
    eikr_r.resize(nptcl, KLists.numk);
    eikr_i.resize(nptcl, KLists.numk);
    FillRhok(P);
  }
#endif
}

} // namespace qmcplusplus
