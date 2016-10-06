//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//                     Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LRCOULOMBAA_H
#define QMCPLUSPLUS_LRCOULOMBAA_H

#include "LongRange/LRHandler.h"
#include "Particle/ParticleSet.h"
#include "LongRange/StructFact.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

template<class BreakupBasis>
class LRCoulombAA: public LRHandler<BreakupBasis>
{
private:
  //Typedef for the lattice-type.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;

  //Import types from the base-class
  typedef LRHandler<BreakupBasis> base_type;
  typedef typename base_type::RealType RealType;
  //Import members from the base-class
  using base_type::Basis;
  using base_type::coefs;
  using base_type::Fk;

  //Private members
  ParticleSet* PtclRef;
  StructFact& PtclRhoK;
  DistanceTableData* d_aa;
  SpeciesSet& tspecies;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NParticles;

  std::vector<RealType> Zat,Zspec;
  std::vector<int> NofSpecies;

public:
  RealType evalLR();
  RealType evalSR();
  RealType evalConsts();

  void initBreakup();
  void resetTargetParticleSet(ParticleSet& ref);

  //Constructor
  LRCoulombAA(ParticleSet& ref):
    //non-default base-class construct. We use 1 function for each species.
    LRHandler<BreakupBasis>(ref.Lattice),
    PtclRhoK(*ref.SK),
    PtclRef(&ref),
    tspecies(ref.getSpeciesSet())
  {
    //Set up the internal distance-table for the particle pair.
    //This is used in evalSR.
    d_aa = DistanceTable::add(*PtclRef);
    //Safety Check: ensure SK was created.
    if(!ref.SK)
    {
      LOGMSG("Structure factor hasn't been created for LRCoulombAA");
      OHMMS::Controller->abort();
    }
    //Things that don't change with lattice are done here instead of InitBreakup()
    ChargeAttribIndx = tspecies.addAttribute("charge");
    MemberAttribIndx = tspecies.addAttribute("membersize");
    NParticles = PtclRef->getTotalNum();
    NumSpecies = tspecies.TotalNum;
    Zat.resize(NParticles);
    Zspec.resize(NumSpecies);
    NofSpecies.resize(NumSpecies);
    for(int spec=0; spec<NumSpecies; spec++)
    {
      Zspec[spec] = tspecies(ChargeAttribIndx,spec);
      NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx,spec));
    }
    for(int iat=0; iat<NParticles; iat++)
      Zat[iat] = Zspec[PtclRef->GroupID[iat]];
    //Initialise the breakup. Can be recalled later if lattice changes.
    //Otherwise all subsequent potential evaluations can reuse existing data.
    //This calls and stores result of evalConsts at end. Ensure that
    //required data is set up already!
    initBreakup();
  }

private:
  RealType evalFk(RealType k,int FunctionIndex);
  RealType evalXk(RealType k,int FunctionIndex);

};

template<class BreakupBasis>
void
LRCoulombAA<BreakupBasis>::initBreakup()
{
  //In this case, all functionality of the base-class can be reused.
  //We only breakup 1 function: The bare interaction for q1=q2=1.
  //We put the charges in manually when the potentials are evaluated.
  LRHandler<BreakupBasis>::InitBreakup(PtclRef->Lattice,1);
  //Additional work:
  //Fill Fk with the FT of V_l(r). This is used in evalLR.
  //This fills Fk for all species. Needs to be repeated if lattice
  //changes - InitBreakup will be called again in that case.
  //GNU C++ needs this-> to access template base-class functions.
  LRHandler<BreakupBasis>::fillFk(PtclRhoK.KLists);
}

//This function is called when the particleset is swapped with another one
//in the Hamiltonian.
template<class BreakupBasis>
void
LRCoulombAA<BreakupBasis>::resetTargetParticleSet(ParticleSet& newP)
{
  PtclRef = &newP;
  //Update the distancetable pointer to use the new Particleset.
  d_aa = DistanceTable::add(*PtclRef);
}


//Evaluate F(k) using basis+coefs.
//This is for the Bare Coulomb Interaction with q1=q2=1
template<class BreakupBasis>
inline typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalFk(RealType k,int FunctionIndex)
{
  RealType FatK;
  FatK = 4.0*M_PI/(Basis.get_CellVolume()*k*k)*
         std::cos(k*Basis.get_rc());
  for(int n=0; n<Basis.NumBasisElem(); n++)
    FatK += coefs[0][n]*Basis.c(n,k);
  return(FatK);
}

//Evaluate x(k) using analytic FT of -V(r) from rc->infinity.
//This is for the Bare Coulomb Interaction with q1=q2=1
template<class BreakupBasis>
inline typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalXk(RealType k,int FunctionIndex)
{
  RealType FatK;
  FatK = -4.0*M_PI/(Basis.get_CellVolume()*k*k)*
         std::cos(k*Basis.get_rc());
  return (FatK);
}


template<class BreakupBasis>
typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalLR()
{
  RealType LR=0.0,temp;
  //Evaluate LR
  //Species loop
  for(int spec1=0; spec1<NumSpecies; spec1++)
  {
    RealType Z1 = Zspec[spec1];
    for(int spec2=0; spec2<NumSpecies; spec2++)
    {
      RealType Z2 = Zspec[spec2];
      temp = 0.0;
      //For each k, find -k and add term to LR.
      for(int ki=0; ki<PtclRhoK.KLists.kpts_cart.size(); ki++)
      {
        int kj=PtclRhoK.KLists.minusk[ki];
        temp += (PtclRhoK.rhok(spec1,ki)*PtclRhoK.rhok(spec2,kj)).real()*Fk[0][ki];
      } //ki
      LR += Z1*Z2*temp;
    } //spec2
  }//spec1
  LR*=0.5;
  return LR;
}


template<class BreakupBasis>
typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalSR()
{
  RealType SR=0.0;
  for(int ipart=0; ipart<NParticles; ipart++)
  {
    RealType esum = 0.0;
    for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++)
    {
      RealType sep=d_aa->r(nn);
      //If r>r_c then skip this pair.
      if(sep >= Basis.get_rc())
        continue;
      if(sep < 1.e-12)
      {
        LOGMSG("WARNING, excluding AA pair with distance < 1.e-12");
        continue;
      }
      //Within cutoff.
      //Now add the short-range pair-wise part.
      double vspair = 0.0;
      //First set equal to the bare-Coulombic interaction with q1=q2=1
      vspair = d_aa->rinv(nn);
      //The subtract off the long-range term (with q1=q2=1):
      //for(int n=0; n<coefs.size(); n++)
      for(int n=0; n<Basis.NumBasisElem(); n++)
        vspair -= coefs[0][n]*Basis.h(n,sep);
      esum += Zat[jpart]*vspair;
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart]*esum;
  }
  return SR;
}


template<class BreakupBasis>
typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalConsts()
{
  RealType Consts=0.0, V0=0.0;
  for(int n=0; n<Basis.NumBasisElem(); n++)
    V0 += coefs[0][n]*Basis.h(n,0.0); //For charge q1=q2=1
  for(int spec=0; spec<NumSpecies; spec++)
  {
    RealType z = Zspec[spec];
    RealType n = NofSpecies[spec];
    Consts += -V0*0.5*z*z*n;
  }
  V0 = Basis.get_rc()*Basis.get_rc()*0.5;
  for(int n=0; n<Basis.NumBasisElem(); n++)
    V0 -= coefs[0][n]*Basis.hintr2(n);
  V0 *= 2.0*TWOPI/Basis.get_CellVolume(); //For charge q1=q2=1
  for(int spec=0; spec<NumSpecies; spec++)
  {
    RealType z = Zspec[spec];
    int n = NofSpecies[spec];
    Consts += -V0*z*z*0.5*n*n;
  }
  //If we have more than one species in this particleset then there is also a
  //single AB term that should be added to the last constant...
  //=-Na*Nb*V0*Za*Zb
  //This accounts for the partitioning of the neutralizing background...
  for(int speca=0; speca<NumSpecies; speca++)
  {
    RealType za = Zspec[speca];
    int na = NofSpecies[speca];
    for(int specb=speca+1; specb<NumSpecies; specb++)
    {
      RealType zb = Zspec[specb];
      int nb = NofSpecies[specb];
      Consts += -V0*za*zb*na*nb;
    }
  }
  return Consts;
}
}
#endif
