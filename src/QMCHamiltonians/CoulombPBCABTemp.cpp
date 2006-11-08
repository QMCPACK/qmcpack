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
    PtclA(&ions), PtclB(&elns), FirstTime(true), myConst(0.0), myGrid(0),V0(0){
      LOGMSG("    Performing long-range breakup for CoulombABTemp potential");
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
    //remove Vspec
    //remove grid
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

    RealType totQ=0.0;
    for(int iat=0; iat<NptclA; iat++)
      totQ+=Zat[iat] = Zspec[PtclA->GroupID[iat]];
    for(int iat=0; iat<NptclB; iat++)
      totQ+=Qat[iat] = Qspec[PtclB->GroupID[iat]];

    if(totQ>numeric_limits<RealType>::epsilon()) {
      LOGMSG("PBCs not yet finished for non-neutral cells");
      OHMMS::Controller->abort();
    }

    ////Test if the box sizes are same (=> kcut same for fixed dimcut)
    kcdifferent = (std::abs(PtclA->Lattice.LR_kc - PtclB->Lattice.LR_kc) > numeric_limits<RealType>::epsilon());
    minkc = std::min(PtclA->Lattice.LR_kc,PtclB->Lattice.LR_kc);

    //AB->initBreakup(*PtclB);
    //initBreakup is called only once
    AB = LRCoulombSingleton::getHandler(*PtclB);
    myConst=evalConsts();
    myRcut=AB->Basis.get_rc();

    if(V0==0) {
      V0 = LRCoulombSingleton::createSpline4RbyVs(AB,myRcut,myGrid);
      if(Vat.size()) {
        app_log() << "  Vat is not empty. Something is wrong" << endl;
        OHMMS::Controller->abort();
      }
      Vat.resize(NptclA,V0);
      Vspec.resize(NumSpeciesA,0);
    }
  }

  void CoulombPBCABTemp::add(int groupID, RadFunctorType* ppot) {
    //add a numerical functor
    if(Vspec[groupID]==0){
      int ng=myGrid->size();
      vector<RealType> v(ng);
      for(int ig=0; ig<ng; ig++) {
        RealType r=(*myGrid)[ig];
        //need to multiply r for the LR
        v[ig]=r*AB->evaluateLR(r)+ppot->splint(r);;
      }
      v[ng-1]=0.0;
      RadFunctorType* rfunc=new RadFunctorType(myGrid,v);
      RealType deriv=(v[1]-v[0])/((*myGrid)[1]-(*myGrid)[0]);
      rfunc->spline(0,deriv,ng-1,0.0);
      Vspec[groupID]=rfunc;
      for(int iat=0; iat<NptclA; iat++) {
        if(PtclA->GroupID[iat]==groupID) Vat[iat]=rfunc;
      }
    }
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
        RadFunctorType* rVs=Vat[iat];
        for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++) {
          //if(d_ab->r(nn)>=myRcut) continue;
          esum += Qat[jat]*d_ab->rinv(nn)*rVs->splint(d_ab->r(nn));
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

      app_log() << "   Constant of PBCAB " << Consts << endl;
      return Consts;
    }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

