#ifndef OHMMS_QMC_LRCOULOMBAB_H
#define OHMMS_QMC_LRCOULOMBAB_H

#include "LongRange/LRHandler.h"
#include "Particle/ParticleSet.h"
#include "LongRange/StructFact.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace ohmmsqmc {

  template<class BreakupBasis>
    class LRCoulombAB: public LRHandler<BreakupBasis> {
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
    ParticleSet& PtclRefA;
    ParticleSet& PtclRefB;
    StructFact& PtclRhoKA;
    StructFact& PtclRhoKB;
    DistanceTableData* d_ab;
    SpeciesSet& tspeciesA;
    SpeciesSet& tspeciesB;
    int NumSpeciesA;
    int NumSpeciesB;
    int ChargeAttribIndxA;
    int ChargeAttribIndxB;
    int MemberAttribIndxA;
    int MemberAttribIndxB;
    int NParticlesA;
    int NParticlesB;
    //This is set to true if the K_c of structure-factors are different
    bool kcdifferent; 
    RealType minkc;

  public:
    RealType evalLR();
    RealType evalSR();
    RealType evalConsts();

    void InitBreakup();

    //Constructor
    LRCoulombAB(ParticleSet& refa,
		ParticleSet& refb): 
      //non-default base-class construct. We use 1 function for each species.
      LRHandler<BreakupBasis>(refa.Lattice), 
      PtclRhoKA(*refa.SK),
      PtclRefA(refa),
      PtclRhoKA(*refb.SK),
      PtclRefA(refb),
      tspeciesA(refa.getSpeciesSet()),
      tspeciesB(refb.getSpeciesSet()) { 

      //Safety check: ensure SK was created.
      if(!refa.SK || !refb.SK){
	LOGMSG("Structure factor hasn't been created for LRCoulombAB");
	OHMMS::Controller->abort();
      }

      //Initialise the breakup. Can be re-called later if lattice changes.
      //Otherwise all subsequent potential evaluations can reuse existing data.
      InitBreakup();

      //Things that don't change with lattice are done here instead of InitBreakup()
      ChargeAttribIndxA = tspeciesA.addAttribute("charge");
      MemberAttribIndxA = tspeciesA.addAttribute("membersize");
      NParticlesA = PtclRefA.getTotalNum();
      ChargeAttribIndxB = tspeciesB.addAttribute("charge");
      MemberAttribIndxB = tspeciesB.addAttribute("membersize");
      NParticlesB = PtclRefB.getTotalNum();

      kcdifferent = false;
      if(fabs(PtclRhoKA.KLists.kcutoff - PtclRhoKB.KLists.kcutoff) < 1.e-6){
	kcdifferent = true;
	minkc = std::min(PtclRhoKA.KLists.kcutoff,PtclRhoKB.KLists.kcutoff);
      }
    }
      
  private:
    RealType evalFk(RealType k,int FunctionIndex);
    RealType evalXk(RealType k,int FunctionIndex);

  };

  template<class BreakupBasis> 
    void 
    LRCoulombAB<BreakupBasis>::InitBreakup() {

      //In this case, all functionality of the base-class can be reused.
      NumSpeciesA = tspeciesA.TotalNum;
      NumSpeciesB = tspeciesB.TotalNum;

      //We only breakup 1 function: The bare interaction for q1=q2=1. 
      //We put the charges in manually when the potentials are evaluated.
      //Breakup the potential using the lattice with smallest volume.
      //Allows extra periodicity of ions to be exploited.
      StructFact *BreakupRhoK;
      ParticleSet *BreakupPtclSet;
      if(PtclRefA.Lattice.Volume < PtclRefB.Lattice.Volume){
        BreakupPtclSet = &PtclRefA;
        BreakupRhoK = &PtclRhoKA;
      } else {
        BreakupPtclSet = &PtclRefB;
        BreakupRhoK = &PtclRhoKB;
      }

      LRHandler<BreakupBasis>::InitBreakup(BreakupPtclSet->Lattice,1); 

      //Additional work:
      //Fill Fk with the FT of V_l(r). This is used in evalLR.
      //This fills Fk for all species.
      LRHandler<BreakupBasis>::fillFk(BreakupRhoK->KLists);

      //Set up the internal distance-table for the particle pair.
      //This is used in evalSR.
      d_ab = DistanceTable::getTable(DistanceTable::add(PtclRefA,PtclRefB));
    }


  //Evaluate F(k) using basis+coefs.
  //This is for the Bare Coulomb Interaction with q1=q2=1
  template<class BreakupBasis> 
    typename LRCoulombAB<BreakupBasis>::RealType 
    LRCoulombAB<BreakupBasis>::evalFk(RealType k,int FunctionIndex) {
      RealType FatK;
      FatK = 4.0*M_PI/(Basis.Lattice.Volume*k*k)*
        cos(k*Basis.get_rc());
      for(int n=0; n<Basis.NumBasisElem(); n++)
        FatK += coefs[0][n]*Basis.c(n,k);
      return(FatK);
    } 

  //Evaluate x(k) using analytic FT of -V(r) from rc->infinity.
  //This is for the Bare Coulomb Interaction with q1=q2=1
  template<class BreakupBasis> 
    typename LRCoulombAB<BreakupBasis>::RealType 
    LRCoulombAB<BreakupBasis>::evalXk(RealType k,int FunctionIndex) {
      RealType FatK;
      FatK = -4.0*M_PI/(Basis.Lattice.Volume*k*k)*
        cos(k*Basis.get_rc());
      return (FatK);
    }


  template<class BreakupBasis>
    typename LRCoulombAB<BreakupBasis>::RealType
    LRCoulombAB<BreakupBasis>::evalLR() {
      RealType LR=0.0,temp;

      //Set up Fk for the k vectors used in the structure factor.
      //This should only be done once unless Lattice changes...
      //GNU C++ needs this-> to access template base-class functions.
      LOGMSG("TODO: Which RhoK to use here?");
      exit(0);
      this->fillFk(PtclRhoKA.KLists);
      //Evaluate LR

      LOGMSG("TODO: Deal with case that StructFacts have different Kc or rc");

      //Species loop. Species loops are over distinct particle-sets.
      if(kcdifferent){
        LOGMSG("Not Coded");
      }
      else {
        for(int speca=0; speca<NumSpeciesA; speca++) {
          RealType Za = tspeciesA(ChargeAttribIndxA,speca);
          for(int specb=0; specb<NumSpeciesB; specb++) {
            RealType Zb = tspeciesB(ChargeAttribIndxB,specb);

            temp = 0.0;
            //For each k, find -k and add term to LR.
            for(int ki=0; ki<PtclRhoKA.KLists.kpts_cart.size(); ki++) {
              int kj=PtclRhoKA.KLists.minusk[ki];
              temp += (PtclRhoKA.rhok(ki,speca)*PtclRhoKB.rhok(kj,specb)).real()*Fk[0][ki];
            } //ki
            LR += Za*Zb*temp;
          } //spec2
        }//spec1
      }

      LR*=0.5;
      return LR;
    }


  template<class BreakupBasis>
    typename LRCoulombAB<BreakupBasis>::RealType
    LRCoulombAB<BreakupBasis>::evalSR() {
      RealType SR=0.0;

      //Update the distance-table
      d_ab->evaluate(PtclRefA,PtclRefB);


      for(int ipart=0; ipart<NParticlesA; ipart++){
        RealType esum = 0.0;
        for(int nn=d_ab->M[ipart], jpart=ipart+1; nn<d_ab->M[ipart+1]; nn++,jpart++) {
          //If r>r_c then skip this pair.
          if(d_ab->r(nn) >= Basis.get_rc())continue;
          if(abs(d_ab->r(nn)) < 1.e-12){
            LOGMSG("WARNING, excluding AB pair with distance < 1.e-12");
            continue;
          }

          //Within cutoff.
          //Now add the short-range pair-wise part.
          double vspair = 0.0; 
          //First set equal to the bare-Coulombic interaction with q1=q2=1
          vspair = d_ab->rinv(nn);
          //The subtract off the long-range term (with q1=q2=1):
          for(int n=0; n<coefs.size(); n++)
            vspair -= coefs[0][n]*Basis.h(n,d_ab->r(nn));

          //Now multiply the species charge for atom j
          esum += tspeciesB(ChargeAttribIndxB,PtclRefB.GroupID[jpart])*vspair;	      
        }
        //Accumulate pair sums...species charge for atom i.
        SR += tspeciesA(ChargeAttribIndxA,PtclRefA.GroupID[ipart])*esum;
      }

      return SR;
    }


  template<class BreakupBasis>
    typename LRCoulombAB<BreakupBasis>::RealType
    LRCoulombAB<BreakupBasis>::evalConsts() {
      RealType Consts=0.0, V0=0.0;


      V0 = Basis.get_rc()*Basis.get_rc()*0.5;
      for(int n=0; n<Basis.NumBasisElem(); n++)
        V0 -= coefs[0][n]*Basis.hintr2(n);
      V0 *= 2.0*TWOPI/PtclRefA.Lattice.Volume; //For charge q1=q2=1


      //Can simplify this if we know a way to get number of particles with each
      //groupID.
      for(int speca=0; speca<NumSpeciesA; speca++) {
        RealType Za = tspeciesA(ChargeAttribIndxA,speca);
        int Na = tspeciesA(MemberAttribIndxA,speca);
        for(int specb=0; specb<NumSpeciesB; specb++) {
          RealType Zb = tspeciesB(ChargeAttribIndxB,specb);
          int Nb = tspeciesB(MemberAttribIndxB,specb);
          Consts += -V0*Za*Zb*Na*Nb;
        }
      }

      return Consts;

    }
}
#endif
