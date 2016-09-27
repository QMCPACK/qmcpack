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


#ifndef QMCPLUSPLUS_LRCOULOMBAB_H
#define QMCPLUSPLUS_LRCOULOMBAB_H

#include "LongRange/LRHandler.h"
#include "Particle/ParticleSet.h"
#include "LongRange/StructFact.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

template<class BreakupBasis>
class LRCoulombAB: public LRHandler<BreakupBasis>
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
  ParticleSet* Ions;
  ParticleSet* Elns;
  StructFact& IonsRhoK;
  StructFact& ElnsRhoK;
  DistanceTableData* d_table;
  SpeciesSet& tspeciesIons;
  SpeciesSet& tspeciesElns;
  int NumSpeciesIons;
  int NumSpeciesElns;
  int ChargeAttribIndxIons;
  int ChargeAttribIndxElns;
  int MemberAttribIndxIons;
  int MemberAttribIndxElns;
  int NIons;
  int NElns;

  std::vector<RealType> Zat,Zspec;
  std::vector<int> NofSpeciesIons;
  std::vector<int> NofSpeciesElns;
  //This is set to true if the K_c of structure-factors are different
  bool kcdifferent;
  RealType minkc;

public:
  RealType evalLR();
  RealType evalSR();
  RealType evalConsts();

  void initBreakup();
  void resetTargetParticleSet(ParticleSet& ref);

  //Constructor
  LRCoulombAB(ParticleSet& ions,
              ParticleSet& elns):
    //non-default base-class construct.
    LRHandler<BreakupBasis>(ions.Lattice),
    IonsRhoK(*ions.SK),
    Ions(&ions),
    ElnsRhoK(*elns.SK),
    Elns(&elns),
    tspeciesIons(ions.getSpeciesSet()),
    tspeciesElns(elns.getSpeciesSet())
  {
    //Set up the internal distance-table for the particle pair.
    //This is used in evalSR.
    d_table = DistanceTable::add(*Ions,*Elns);
    //Safety check: ensure SK was created.
    if(!ions.SK || !elns.SK)
    {
      LOGMSG("Structure factor hasn't been created for LRCoulombAB");
      OHMMS::Controller->abort();
    }
    //Things that don't change with lattice are done here instead of InitBreakup()
    ChargeAttribIndxIons = tspeciesIons.addAttribute("charge");
    MemberAttribIndxIons = tspeciesIons.addAttribute("membersize");
    ChargeAttribIndxElns = tspeciesElns.addAttribute("charge");
    MemberAttribIndxElns = tspeciesElns.addAttribute("membersize");
    //Store total number of particles in each set.
    NIons = Ions->getTotalNum();
    NElns = Elns->getTotalNum();
    if(NElns != NIons)
    {
      LOGMSG("PBCs not yet finished for non-neutral cells");
      OHMMS::Controller->abort();
    }
    //Store the number of unique species in each set.
    NumSpeciesIons = tspeciesIons.TotalNum;
    NumSpeciesElns = tspeciesElns.TotalNum;
    //Store information about charges and number of each species
    Zat.resize(NIons);
    Zspec.resize(NIons);
    NofSpeciesIons.resize(NIons);
    NofSpeciesElns.resize(NElns);
    for(int spec=0; spec<NumSpeciesIons; spec++)
    {
      Zspec[spec] = tspeciesIons(ChargeAttribIndxIons,spec);
      NofSpeciesIons[spec] = static_cast<int>(tspeciesIons(MemberAttribIndxIons,spec));
    }
    for(int spec=0; spec<NumSpeciesElns; spec++)
      NofSpeciesElns[spec] = static_cast<int>(tspeciesElns(MemberAttribIndxElns,spec));
    for(int iat=0; iat<NIons; iat++)
      Zat[iat] = Zspec[Ions->GroupID[iat]];
    //Test if the box sizes are same (=> kcut same for fixed dimcut)
    kcdifferent = false;
    if(std::abs(Ions->Lattice.LR_kc - Elns->Lattice.LR_kc) > 1.e-6)
    {
      kcdifferent = true;
      minkc = std::min(Ions->Lattice.LR_kc,Elns->Lattice.LR_kc);
    }
    //Initialise the breakup. Can be re-called later if lattice changes.
    //Otherwise all subsequent potential evaluations can reuse existing data.
    initBreakup();
  }

private:
  RealType evalFk(RealType k,int FunctionIndex);
  RealType evalXk(RealType k,int FunctionIndex);

};

template<class BreakupBasis>
void
LRCoulombAB<BreakupBasis>::initBreakup()
{
  //In this case, all functionality of the base-class can be reused.
  //We only breakup 1 function: The bare interaction for q1=q2=1.
  //We put the charges in manually when the potentials are evaluated.
  //Breakup the potential using the lattice with smallest volume.
  //Allows extra periodicity of ions to be exploited.
  StructFact *BreakupRhoK;
  ParticleSet *BreakupPtclSet;
  if(Ions->Lattice.Volume < Elns->Lattice.Volume)
  {
    BreakupPtclSet = Ions;
    BreakupRhoK = &IonsRhoK;
  }
  else
  {
    BreakupPtclSet = Elns;
    BreakupRhoK = &ElnsRhoK;
  }
  LRHandler<BreakupBasis>::InitBreakup(BreakupPtclSet->Lattice,1);
  //Additional work:
  //Fill Fk with the FT of V_l(r). This is used in evalLR.
  //This fills Fk for all species.
  LRHandler<BreakupBasis>::fillFk(BreakupRhoK->KLists);
}

//This function is called when the particleset is swapped with another one
//in the Hamiltonian.
template<class BreakupBasis>
void
LRCoulombAB<BreakupBasis>::resetTargetParticleSet(ParticleSet& newP)
{
  //"Target" particleset is always electrons
  Elns = &newP;
  //Update the distancetable pointer to use the new Particleset.
  d_table = DistanceTable::add(*Ions,*Elns);
}


//Evaluate F(k) using basis+coefs.
//This is for the Bare Coulomb Interaction with q1=q2=1
template<class BreakupBasis>
typename LRCoulombAB<BreakupBasis>::RealType
LRCoulombAB<BreakupBasis>::evalFk(RealType k,int FunctionIndex)
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
typename LRCoulombAB<BreakupBasis>::RealType
LRCoulombAB<BreakupBasis>::evalXk(RealType k,int FunctionIndex)
{
  RealType FatK;
  FatK = -4.0*M_PI/(Basis.get_CellVolume()*k*k)*
         std::cos(k*Basis.get_rc());
  return (FatK);
}


template<class BreakupBasis>
typename LRCoulombAB<BreakupBasis>::RealType
LRCoulombAB<BreakupBasis>::evalLR()
{
  RealType LR=0.0,temp;
  //Set up Fk for the k vectors used in the structure factor.
  //This should only be done once unless Lattice changes...
  //GNU C++ needs this-> to access template base-class functions.
  //Species loop. Species loops are over distinct particle-sets.
  if(kcdifferent)
  {
    LOGMSG("LRCoulombAB: Different cell-sizes not yet coded");
    exit(0);
  }
  else
  {
    for(int specion=0; specion<NumSpeciesIons; specion++)
    {
      RealType Zion = Zspec[specion];
      for(int speceln=0; speceln<NumSpeciesElns; speceln++)
      {
        RealType Zeln = -1.0; //All electrons have same charge.
        temp = 0.0;
        //For each k, find -k and add term to LR.
        for(int ki=0; ki<IonsRhoK.KLists.kpts_cart.size(); ki++)
        {
          int kj=IonsRhoK.KLists.minusk[ki];
          temp += (IonsRhoK.rhok(specion,ki)*ElnsRhoK.rhok(speceln,kj)).real()*Fk[0][ki];
        } //ki
        LR += Zion*Zeln*temp;
      } //spec2
    }//spec1
  }
  return LR;
}


template<class BreakupBasis>
typename LRCoulombAB<BreakupBasis>::RealType
LRCoulombAB<BreakupBasis>::evalSR()
{
  RealType SR=0.0;
  //Loop over distinct eln-ion pairs
  for(int iat=0; iat<NIons; iat++)
  {
    RealType esum = 0.0;
    for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++)
    {
      RealType sep = d_table->r(nn);
      //If r>r_c then skip this pair.
      if(sep >= Basis.get_rc())
        continue;
      if(sep < 1.e-12)
      {
        LOGMSG("WARNING, excluding AB pair with distance < 1.e-12");
        continue;
      }
      //Within cutoff.
      //Now add the short-range pair-wise part.
      RealType vspair = 0.0;
      //First set equal to the bare-Coulombic interaction with q1=q2=1
      vspair = d_table->rinv(nn);
      //The subtract off the long-range term (with q1=q2=1):
      //for(int n=0; n<coefs.size(); n++)
      for(int n=0; n<Basis.NumBasisElem(); n++)
        vspair -= coefs[0][n]*Basis.h(n,sep);
      //Now multiply the species charge for atom j
      //Particle is an eln => -1 charge.
      esum += -1.0*vspair;
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[iat]*esum;
  }
  return SR;
}


template<class BreakupBasis>
typename LRCoulombAB<BreakupBasis>::RealType
LRCoulombAB<BreakupBasis>::evalConsts()
{
  RealType Consts=0.0, V0=0.0;
  V0 = Basis.get_rc()*Basis.get_rc()*0.5;
  for(int n=0; n<Basis.NumBasisElem(); n++)
    V0 -= coefs[0][n]*Basis.hintr2(n);
  V0 *= 2.0*TWOPI/Basis.get_CellVolume(); //For charge q1=q2=1
  //Can simplify this if we know a way to get number of particles with each
  //groupID.
  for(int specion=0; specion<NumSpeciesIons; specion++)
  {
    RealType Zion = Zspec[specion];
    int Nion = NofSpeciesIons[specion];
    for(int speceln=0; speceln<NumSpeciesElns; speceln++)
    {
      RealType Zeln = -1.0; //All electrons have same charge
      int Neln = NofSpeciesElns[speceln];
      Consts += -V0*Zion*Zeln*Nion*Neln;
    }
  }
  return Consts;
}
}
#endif
