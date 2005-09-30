#ifndef OHMMS_QMC_LRCOULOMBAA_H
#define OHMMS_QMC_LRCOULOMBAA_H

#include "LongRange/LRHandler.h"
#include "Particle/ParticleSet.h"
#include "LongRange/StructFact.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace ohmmsqmc {

  template<class BreakupBasis>
    class LRCoulombAA: public LRHandler<BreakupBasis> {
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
    ParticleSet& PtclRef;
    StructFact& PtclRhoK;
    DistanceTableData* d_aa;
    SpeciesSet& tspecies;
    int NumSpecies;
    int ChargeAttribIndx;
    int MemberAttribIndx;
    int NParticles;

  public:
    RealType evalLR();
    RealType evalSR();
    RealType evalConsts();

    void InitBreakup();

    //Constructor
    LRCoulombAA(ParticleSet& ref): 
      //non-default base-class construct. We use 1 function for each species.
      LRHandler<BreakupBasis>(ref.Lattice), 
      PtclRhoK(*ref.SK),
      PtclRef(ref),
      tspecies(ref.getSpeciesSet()) { 

      //Initialise the breakup. Can be recalled later if lattice changes.
      //Otherwise all subsequent potential evaluations can reuse existing data.
      InitBreakup();
      
      //Things that don't change with lattice are done here instead of InitBreakup()
      ChargeAttribIndx = tspecies.addAttribute("charge");
      MemberAttribIndx = tspecies.addAttribute("membersize");
      NParticles = PtclRef.getTotalNum();
    }
      
  private:
    RealType evalFk(RealType k,int FunctionIndex);
    RealType evalXk(RealType k,int FunctionIndex);

  };

template<class BreakupBasis> 
void 
LRCoulombAA<BreakupBasis>::InitBreakup() {

  //In this case, all functionality of the base-class can be reused.
  NumSpecies = tspecies.TotalNum;

  //We only breakup 1 function: The bare interaction for q1=q2=1. 
  //We put the charges in manually when the potentials are evaluated.
  LRHandler<BreakupBasis>::InitBreakup(PtclRef.Lattice,1); 

  //Additional work:
  //Fill Fk with the FT of V_l(r). This is used in evalLR.
  //This fills Fk for all species.
  LRHandler<BreakupBasis>::fillFk(PtclRhoK.KLists);
  
  //Set up the internal distance-table for the particle pair.
  //This is used in evalSR.
  d_aa = DistanceTable::getTable(DistanceTable::add(PtclRef));
  //TODO: this is only needed because ee hasn't been added to
  //DistanceTable:: set before create is called.
  d_aa->create(1);
  // DistanceTable::update(PtclRef);
}


//Evaluate F(k) using basis+coefs.
//This is for the Bare Coulomb Interaction with q1=q2=1
template<class BreakupBasis> 
typename LRCoulombAA<BreakupBasis>::RealType 
LRCoulombAA<BreakupBasis>::evalFk(RealType k,int FunctionIndex) {
  RealType FatK;
  FatK = 4.0*M_PI/(PtclRef.Lattice.Volume*k*k)*
    cos(k*Basis.get_rc());
  for(int n=0; n<Basis.NumBasisElem(); n++)
    FatK += coefs[0][n]*Basis.c(n,k);
  return(FatK);
} 

//Evaluate x(k) using analytic FT of -V(r) from rc->infinity.
//This is for the Bare Coulomb Interaction with q1=q2=1
template<class BreakupBasis> 
typename LRCoulombAA<BreakupBasis>::RealType 
LRCoulombAA<BreakupBasis>::evalXk(RealType k,int FunctionIndex) {
  RealType FatK;
  FatK = -4.0*M_PI/(PtclRef.Lattice.Volume*k*k)*
    cos(k*Basis.get_rc());
  return (FatK);
}


template<class BreakupBasis>
typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalLR() {
  RealType LR=0.0,temp;

  //Set up Fk for the k vectors used in the structure factor.
  //This should only be done once!!
  //Unless Lattice changes...
  //GNU C++ needs this-> to access template base-class functions.
  this->fillFk(PtclRhoK.KLists);
  //Evaluate LR
  //Species loop
  for(int spec1=0; spec1<NumSpecies; spec1++) {
    RealType Z1 = tspecies(ChargeAttribIndx,spec1);
    for(int spec2=0; spec2<NumSpecies; spec2++) {
      RealType Z2 = tspecies(ChargeAttribIndx,spec2);
      
      temp = 0.0;
      //For each k, find -k and add term to LR.
      for(int ki=0; ki<PtclRhoK.KLists.kpts_cart.size(); ki++) {
	int kj=PtclRhoK.KLists.minusk[ki];
	temp += (PtclRhoK.rhok(ki,spec1)*PtclRhoK.rhok(kj,spec2)).real()*Fk[0][ki];
      } //ki
      LR += Z1*Z2*temp;

    } //spec2
  }//spec1

  LR*=0.5;
  return LR;
}


template<class BreakupBasis>
typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalSR() {
  RealType SR=0.0;

  //Update the distance-table
  d_aa->evaluate(PtclRef);

  for(int ipart=0; ipart<NParticles; ipart++){

    
    RealType esum = 0.0;
    for(int nn=d_aa->M[ipart], jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++) {
      //If r>r_c then skip this pair.
      if(d_aa->r(nn) >= Basis.get_rc())continue;
      if(abs(d_aa->r(nn)) < 1.e-12){
	LOGMSG("WARNING, excluding AB pair with distance < 1.e-12");
	continue;
      }
      
      //Within cutoff.
      //Now add the short-range pair-wise part.
      double vspair = 0.0; 
      //First set equal to the bare-Coulombic interaction with q1=q2=1
      vspair = d_aa->rinv(nn);
      //The subtract off the long-range term (with q1=q2=1):
      for(int n=0; n<coefs.size(); n++)
	vspair -= coefs[0][n]*Basis.h(n,d_aa->r(nn));
      
      //Now multiply the species charge for atom j
      esum += tspecies(ChargeAttribIndx,PtclRef.GroupID[jpart])*vspair;	      
    }
    //Accumulate pair sums...species charge for atom i.
    SR += tspecies(ChargeAttribIndx,PtclRef.GroupID[ipart])*esum;
  }

  return SR;
}


template<class BreakupBasis>
typename LRCoulombAA<BreakupBasis>::RealType
LRCoulombAA<BreakupBasis>::evalConsts() {
  RealType Consts=0.0, V0=0.0;
  
  for(int n=0; n<Basis.NumBasisElem(); n++)
    V0 += coefs[0][n]*Basis.h(n,0.0); //For charge q1=q2=1

  for(int spec=0; spec<NumSpecies; spec++) {
    //  for(int i=0; i<NParticles; i++) {
    RealType Z = tspecies(ChargeAttribIndx,spec);
    RealType NofSpecies = tspecies(MemberAttribIndx,spec);
    Consts += -V0*0.5*Z*Z*NofSpecies;
  }

  V0 = Basis.get_rc()*Basis.get_rc()*0.5;
  for(int n=0; n<Basis.NumBasisElem(); n++)
    V0 -= coefs[0][n]*Basis.hintr2(n);
  V0 *= 2.0*TWOPI/PtclRef.Lattice.Volume; //For charge q1=q2=1
  
  //Can simplify this if we know a way to get number of particles with each
  //groupID.
  
  for(int spec=0; spec<NumSpecies; spec++){
    RealType Z = tspecies(ChargeAttribIndx,spec);
    RealType NofSpecies = tspecies(MemberAttribIndx,spec);
    Consts += -V0*Z*Z*0.5*NofSpecies*NofSpecies;
  }
  
  return Consts;

}
}
#endif
