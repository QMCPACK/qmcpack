//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/CoulombPBCABTemp.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"

namespace qmcplusplus {

  CoulombPBCABTemp::CoulombPBCABTemp(ParticleSet& ions, ParticleSet& elns): 
    PtclA(&ions), PtclB(&elns), FirstTime(true), myConst(0.0){
      LOGMSG("  Performing long-range breakup for CoulombAATemp potential");
      //Use singleton pattern 
      //AB = new LRHandlerType(ions);
      d_ab = DistanceTable::add(ions,elns);
      initBreakup();
      LOGMSG("    Done\n");
    }

  /// copy constructor
  CoulombPBCABTemp::CoulombPBCABTemp(const CoulombPBCABTemp& c): 
    PtclA(c.PtclA),PtclB(c.PtclB),d_ab(c.d_ab), FirstTime(true), myConst(0.0){
      initBreakup();
    }
    
  CoulombPBCABTemp:: ~CoulombPBCABTemp() { 
  }

  void CoulombPBCABTemp::resetTargetParticleSet(ParticleSet& P) {
    //Update the internal particleref
    PtclB = &P;
    d_ab = DistanceTable::add(*PtclA,P);
    AB->resetTargetParticleSet(P);
  }

  CoulombPBCABTemp::Return_t 
    CoulombPBCABTemp::evaluate(ParticleSet& P) {
      return Value = evalLR()+evalSR()+myConst;
    }


  void CoulombPBCABTemp::initBreakup() {
    SpeciesSet& tspeciesA(PtclA->getSpeciesSet());
    SpeciesSet& tspeciesB(PtclB->getSpeciesSet());

    ChargeAttribIndxA = tspeciesA.addAttribute("charge");
    MemberAttribIndxA = tspeciesA.addAttribute("membersize");
    ChargeAttribIndxB = tspeciesB.addAttribute("charge");
    MemberAttribIndxB = tspeciesB.addAttribute("membersize");

    NptclA = PtclA->getTotalNum();
    NptclB = PtclB->getTotalNum();

    if(NptclA != NptclB) {
      LOGMSG("PBCs not yet finished for non-neutral cells");
      OHMMS::Controller->abort();
    }

    NumSpeciesA = tspeciesA.TotalNum;
    NumSpeciesB = tspeciesB.TotalNum;

    //Store information about charges and number of each species
    Zat.resize(NptclA); Zspec.resize(NumSpeciesA);
    Qat.resize(NptclB); Qspec.resize(NumSpeciesB);

    NofSpeciesA.resize(NumSpeciesA);
    NofSpeciesB.resize(NumSpeciesB);

    for(int spec=0; spec<NumSpeciesA; spec++) { 
      Zspec[spec] = tspeciesA(ChargeAttribIndxA,spec);
      NofSpeciesA[spec] = static_cast<int>(tspeciesA(MemberAttribIndxA,spec));
    }
    for(int spec=0; spec<NumSpeciesB; spec++) {
      Qspec[spec] = tspeciesB(ChargeAttribIndxB,spec);
      NofSpeciesB[spec] = static_cast<int>(tspeciesB(MemberAttribIndxB,spec));
    }
    for(int iat=0; iat<NptclA; iat++)
      Zat[iat] = Zspec[PtclA->GroupID[iat]];
    for(int iat=0; iat<NptclB; iat++)
      Qat[iat] = Qspec[PtclB->GroupID[iat]];

    ////Test if the box sizes are same (=> kcut same for fixed dimcut)
    kcdifferent = (std::abs(PtclA->Lattice.LR_kc - PtclB->Lattice.LR_kc) > 1e-6);
    minkc = std::min(PtclA->Lattice.LR_kc,PtclB->Lattice.LR_kc);

    //AB->initBreakup(*PtclB);
    //initBreakup is called only once
    AB = LRCoulombSingleton::getHandler(*PtclB);
    myConst=evalConsts();
    myRcut=AB->Basis.get_rc();
  }

  CoulombPBCABTemp::Return_t
    CoulombPBCABTemp::evalLR() {
      RealType LR=0.0;
      const StructFact& RhoKA(*(PtclA->SK));
      const StructFact& RhoKB(*(PtclB->SK));
      for(int i=0; i<NumSpeciesA; i++) {
        RealType esum=0.0;
        for(int j=0; j<NumSpeciesB; j++) {
          esum += Qspec[j]*AB->evaluate(RhoKA.KLists.minusk, RhoKA.rhok[i],RhoKB.rhok[j]);
        } //speceln
        LR += Zspec[i]*esum;
      }//specion
      return LR;
    }

  CoulombPBCABTemp::Return_t
    CoulombPBCABTemp::evalSR() {
      RealType SR=0.0;
      //Loop over distinct eln-ion pairs
      for(int iat=0; iat<NptclA; iat++){
        RealType esum = 0.0;
        for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++) {
          if(d_ab->r(nn)>=myRcut) continue;
          esum += Qat[jat]*AB->evaluate(d_ab->r(nn),d_ab->rinv(nn));
        }
        //Accumulate pair sums...species charge for atom i.
        SR += Zat[iat]*esum;
      }
      return SR;
    }

  CoulombPBCABTemp::Return_t
    CoulombPBCABTemp::evalConsts() {
      LRHandlerType::BreakupBasisType &Basis(AB->Basis);
      const Vector<RealType> &coefs(AB->coefs);
      RealType V0 = Basis.get_rc()*Basis.get_rc()*0.5;
      for(int n=0; n<coefs.size(); n++)
        V0 -= coefs[n]*Basis.hintr2(n);
      V0 *= 2.0*TWOPI/Basis.get_CellVolume(); //For charge q1=q2=1

      //Can simplify this if we know a way to get number of particles with each
      //groupID.
      RealType Consts=0.0;
      for(int i=0; i<NumSpeciesA; i++) {
        RealType q=Zspec[i]*NofSpeciesA[i];
        for(int j=0; j<NumSpeciesB; j++) {
          Consts += -V0*Qspec[j]*NofSpeciesB[j]*q;
        }
      }
      return Consts;
    }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

