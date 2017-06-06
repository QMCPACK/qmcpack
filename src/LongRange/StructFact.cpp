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
#include <config/stdlib/math.h>
#include <Numerics/e2iphi.h>
#include <simd/vmath.hpp>
#include <Numerics/OhmmsBlas.h>
#include <qmc_common.h>

namespace qmcplusplus
{

//Constructor - pass arguments to KLists' constructor
StructFact::StructFact(ParticleSet& P, RealType kc):
  DoUpdate(false), StorePerParticle(false), SuperCellEnum(SUPERCELL_BULK)
{
  if(qmc_common.use_ewald && P.LRBox.SuperCellEnum == SUPERCELL_SLAB)
  {
    app_log() << "  Setting StructFact::SuperCellEnum=SUPERCELL_SLAB " << std::endl;
    SuperCellEnum=SUPERCELL_SLAB;
  }

  UpdateNewCell(P,kc);
}

//Destructor
StructFact::~StructFact() { }

void
StructFact::UpdateNewCell(ParticleSet& P, RealType kc)
{
  //Generate the lists of k-vectors
  KLists.UpdateKLists(P.LRBox,kc);
  //resize any arrary
  resize(P.getSpeciesSet().size(),P.getTotalNum(),KLists.numk);
  //Compute the entire Rhok
  FillRhok(P);
}

void StructFact::resize(int ns, int nptcl, int nkpts)
{
  phiV.resize(nkpts);
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r.resize(ns,nkpts);
  rhok_i.resize(ns,nkpts);
  if(StorePerParticle)
  {
    eikr_r.resize(nptcl,nkpts);
    eikr_i.resize(nptcl,nkpts);
  }
  eikr_r_temp.resize(nkpts);
  eikr_i_temp.resize(nkpts);
#else
  rhok.resize(ns,nkpts);
  eikr.resize(nptcl,nkpts);
  eikr_temp.resize(nkpts);
#endif
  int maxdim=KLists.mmax[DIM];
  C.resize(DIM,2*maxdim+1);
}



void
StructFact::UpdateAllPart(ParticleSet& P)
{
  FillRhok(P);
}


/** evaluate rok per species, eikr  per particle
 */
void
StructFact::FillRhok(ParticleSet& P)
{
  int npart = P.getTotalNum();
#if defined(QMC_SK_USE_RECURSIVE)
  rhok=0.0;
  for(int i=0; i<npart; i++)
  {
    //operate with a reduced positon
    PosType tau_red=P.LRBox.toUnit(P.R[i]);
    for(int idim=0; idim<DIM; idim++)
    {
      RealType phi=TWOPI*tau_red[idim];
      ComplexType ctemp(std::cos(phi),std::sin(phi));
      C(idim,KLists.mmax[idim])=1.0;
      for(int n=1; n<=KLists.mmax[idim]; n++)
      {
        C(idim,KLists.mmax[idim]+n) = ctemp*C(idim,KLists.mmax[idim]+n-1);
        C(idim,KLists.mmax[idim]-n) = conj(C(idim,KLists.mmax[idim]+n));
      }
    }
    ComplexType* restrict eikr_ref=eikr[i];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      eikr_ref[ki]=C(0,KLists.kpts[ki][0]+KLists.mmax[0]);
      for(idim=1; idim<DIM; id++)
        eikr_ref[ki] *= C(idim,KLists.kpts[ki][idim]+KLists.mmax[idim]);
    }
    accumulate_elements(eikr_ref,eikr_ref+KLists.numk,rhok[P.GroupID[i]]);
  } //End particle loop
#else
#if defined(USE_REAL_STRUCT_FACTOR)
  rhok_r=0.0;
  rhok_i=0.0;
  //algorithmA
  const int nk=KLists.numk;
  if(StorePerParticle)
  {
    // save per particle and species value
    for(int i=0; i<npart; ++i)
    {
      //defined in this file
      simd::get_phase(nk,&(KLists.kpts_cart[0][0]), P.R[i].data(), phiV.data());
      //get_phase simply encapsulate this
      eval_e2iphi(nk, phiV.data(), eikr_r[i], eikr_i[i]);
      simd::add(nk,eikr_r[i],rhok_r[P.GroupID[i]]);
      simd::add(nk,eikr_i[i],rhok_i[P.GroupID[i]]);
    }
  }
  else
  {
    // save per species value
    for(int i=0; i<npart; ++i)
    {
      //defined in this file
      simd::get_phase(nk,&(KLists.kpts_cart[0][0]), P.R[i].data(), phiV.data());
      //get_phase simply encapsulate this
      eval_e2iphi(nk, phiV.data(), eikr_r_temp.data(), eikr_i_temp.data());
      simd::add(nk,eikr_r_temp.data(),rhok_r[P.GroupID[i]]);
      simd::add(nk,eikr_i_temp.data(),rhok_i[P.GroupID[i]]);
    }
  }
  //use dgemm: vtune shows algorithmA is better
#else
  rhok=0.0;
  for(int i=0; i<npart; i++)
  {
    PosType pos(P.R[i]);
    RealType s,c;//get sin and cos
    ComplexType* restrict eikr_ref=eikr[i];
    ComplexType* restrict rhok_ref=rhok[P.GroupID[i]];
    for(int ki=0; ki<KLists.numk; ki++)
    {
      sincos(dot(KLists.kpts_cart[ki],pos),&s,&c);
      eikr_ref[ki]=ComplexType(c,s);
      rhok_ref[ki]+= eikr_ref[ki];
    }
  }
#endif
#endif
}


void
StructFact::UpdateRhok(const PosType& rold,const PosType& rnew,int iat,int GroupID)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  APP_ABORT("WHO IS USING UpdateRhok");
#else
#endif
}

void StructFact::makeMove(int active, const PosType& pos)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  for(int ki=0; ki<KLists.numk; ki++)
    phiV[ki]=dot(KLists.kpts_cart[ki],pos);
  eval_e2iphi(KLists.numk, phiV.data(), eikr_r_temp.data(), eikr_i_temp.data());
#else
  RealType s,c;//get sin and cos
  for(int ki=0; ki<KLists.numk; ++ki)
  {
    sincos(dot(KLists.kpts_cart[ki],pos),&s,&c);
    eikr_temp[ki]=ComplexType(c,s);
  }
#endif
}

void StructFact::acceptMove(int active, int gid, const PosType& rold)
{
#if defined(USE_REAL_STRUCT_FACTOR)
  if(StorePerParticle)
  {
    RealType* restrict eikr_ptr_r=eikr_r[active];
    RealType* restrict eikr_ptr_i=eikr_i[active];
    RealType* restrict rhok_ptr_r(rhok_r[gid]);
    RealType* restrict rhok_ptr_i(rhok_i[gid]);
    for(int ki=0; ki<KLists.numk; ++ki)
    {
      rhok_ptr_r[ki] += (eikr_r_temp[ki]-eikr_ptr_r[ki]);
      rhok_ptr_i[ki] += (eikr_i_temp[ki]-eikr_ptr_i[ki]);
      eikr_ptr_r[ki]=eikr_r_temp[ki];
      eikr_ptr_i[ki]=eikr_i_temp[ki];
    }
  }
  else
  {
    RealType* restrict rhok_ptr_r(rhok_r[gid]);
    RealType* restrict rhok_ptr_i(rhok_i[gid]);

    // add the new value
    for(int ki=0; ki<KLists.numk; ++ki)
    {
      rhok_ptr_r[ki] += eikr_r_temp[ki];
      rhok_ptr_i[ki] += eikr_i_temp[ki];
    }

    for(int ki=0; ki<KLists.numk; ki++)
      phiV[ki]=dot(KLists.kpts_cart[ki],rold);
    eval_e2iphi(KLists.numk, phiV.data(), eikr_r_temp.data(), eikr_i_temp.data());

    // remove the old value
    for(int ki=0; ki<KLists.numk; ++ki)
    {
      rhok_ptr_r[ki] -= eikr_r_temp[ki];
      rhok_ptr_i[ki] -= eikr_i_temp[ki];
    }
  }
#else
  ComplexType* restrict eikr_ptr=eikr[active];
  ComplexType* restrict rhok_ptr(rhok[gid]);
  for(int ki=0; ki<KLists.numk; ++ki)
  {
    rhok_ptr[ki] += (eikr_temp[ki]-eikr_ptr[ki]);
    eikr_ptr[ki]=eikr_temp[ki];
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
  if(!StorePerParticle)
  {
    StorePerParticle=true;
    const int nptcl=P.getTotalNum();
    eikr_r.resize(nptcl,KLists.numk);
    eikr_i.resize(nptcl,KLists.numk);
    FillRhok(P);
  }
#endif
}

}
