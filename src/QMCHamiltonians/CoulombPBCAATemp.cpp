//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/CoulombPBCAATemp.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus {

  CoulombPBCAATemp::CoulombPBCAATemp(ParticleSet& ref, bool active): 
    PtclRef(&ref), AA(0), myGrid(0), rVs(0), 
    is_active(active), FirstTime(true), myConst(0.0)
  {
    ReportEngine PRE("CoulombPBCAATemp","CoulombPBCAATemp");
    //AA = new LRHandlerType(ref);
    d_aa = DistanceTable::add(ref);
    initBreakup();
    app_log() << "  Maximum K shell " << AA->MaxKshell << endl;
    app_log() << "  Number of k vectors " << AA->Fk.size() << endl;

    if(!active) 
    {
      Value = evalLR()+evalSR()+myConst;
      app_log() << "  Constant pair-potential " << ref.getName() << " " << Value << endl;
    }
  }


  ///// copy constructor
  //CoulombPBCAATemp::CoulombPBCAATemp(const CoulombPBCAATemp& c): 
  //  PtclRef(c.PtclRef),d_aa(c.d_aa),myGrid(0),rVs(0), FirstTime(true), myConst(0.0)
  //  {
  //    //AA = new LRHandlerType(*PtclRef);
  //    initBreakup();
  //  }
    
  CoulombPBCAATemp:: ~CoulombPBCAATemp() { }

  void CoulombPBCAATemp::resetTargetParticleSet(ParticleSet& P) {
    if(is_active)
    {
      //Update the internal particleref
      PtclRef = &P;
      d_aa = DistanceTable::add(P);
      AA->resetTargetParticleSet(P);
    }
  }

  CoulombPBCAATemp::Return_t 
    CoulombPBCAATemp::evaluate(ParticleSet& P) {
      if(is_active) Value = evalLR()+evalSR()+myConst;
      return Value;
    }


  void CoulombPBCAATemp::initBreakup() {
    SpeciesSet& tspecies(PtclRef->getSpeciesSet());
    //Things that don't change with lattice are done here instead of InitBreakup()
    ChargeAttribIndx = tspecies.addAttribute("charge");
    MemberAttribIndx = tspecies.addAttribute("membersize");
    NParticles = PtclRef->getTotalNum();
    NumSpecies = tspecies.TotalNum;

    Zat.resize(NParticles);
    Zspec.resize(NumSpecies);
    NofSpecies.resize(NumSpecies);
    for(int spec=0; spec<NumSpecies; spec++) {
      Zspec[spec] = tspecies(ChargeAttribIndx,spec);
      NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx,spec));
    }
    for(int iat=0; iat<NParticles; iat++)
      Zat[iat] = Zspec[PtclRef->GroupID[iat]];

    AA = LRCoulombSingleton::getHandler(*PtclRef);
    //AA->initBreakup(*PtclRef);
    myConst=evalConsts();
    myRcut=AA->Basis.get_rc();
    if(rVs==0) {
      rVs = LRCoulombSingleton::createSpline4RbyVs(AA,myRcut,myGrid);
    }
  }

  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalLR() {
      RealType LR=0.0;
      const StructFact& PtclRhoK(*(PtclRef->SK));
      for(int spec1=0; spec1<NumSpecies; spec1++) {
        RealType Z1 = Zspec[spec1];
        for(int spec2=spec1; spec2<NumSpecies; spec2++) {
          RealType Z2 = Zspec[spec2];
          //RealType temp=AA->evaluate(PtclRhoK.KLists.minusk, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
          RealType temp=AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
          if(spec2==spec1)
            LR += 0.5*Z1*Z2*temp;    
          else
            LR += Z1*Z2*temp;
        } //spec2
      }//spec1
      //LR*=0.5;
      return LR;
    }

  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalSR() {
      RealType SR=0.0;
      for(int ipart=0; ipart<NParticles; ipart++){
        RealType esum = 0.0;
        for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++) {
          //if(d_aa->r(nn)>=myRcut) continue;
          //esum += Zat[jpart]*AA->evaluate(d_aa->r(nn),d_aa->rinv(nn));
          esum += Zat[jpart]*d_aa->rinv(nn)*rVs->splint(d_aa->r(nn));
        }
        //Accumulate pair sums...species charge for atom i.
        SR += Zat[ipart]*esum;
      }
      return SR;
    }

  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalConsts() {

      LRHandlerType::BreakupBasisType &Basis(AA->Basis);
      const Vector<RealType> &coefs(AA->coefs);
      RealType Consts=0.0, V0=0.0;

      for(int n=0; n<coefs.size(); n++)
        V0 += coefs[n]*Basis.h(n,0.0); //For charge q1=q2=1

      for(int spec=0; spec<NumSpecies; spec++) {
        RealType z = Zspec[spec];
        RealType n = NofSpecies[spec];
        Consts += -V0*0.5*z*z*n;
      }

      V0 = Basis.get_rc()*Basis.get_rc()*0.5;
      for(int n=0; n<Basis.NumBasisElem(); n++)
        V0 -= coefs[n]*Basis.hintr2(n);
      V0 *= 2.0*TWOPI/Basis.get_CellVolume(); //For charge q1=q2=1

      for(int spec=0; spec<NumSpecies; spec++){
        RealType z = Zspec[spec];
        int n = NofSpecies[spec];
        Consts += -V0*z*z*0.5*n*n;
      }

      //If we have more than one species in this particleset then there is also a 
      //single AB term that should be added to the last constant...
      //=-Na*Nb*V0*Za*Zb
      //This accounts for the partitioning of the neutralizing background...
      for(int speca=0;speca<NumSpecies;speca++) {
        RealType za = Zspec[speca];
        int na = NofSpecies[speca];
        for(int specb=speca+1;specb<NumSpecies;specb++) {
          RealType zb = Zspec[specb];
          int nb = NofSpecies[specb];
          Consts += -V0*za*zb*na*nb;
        }
      }

      app_log() << "   Constant of PBCAA " << Consts << endl;
      return Consts;
    }

    QMCHamiltonianBase* CoulombPBCAATemp::makeClone(ParticleSet& qp, TrialWaveFunction& psi) 
    {
      if(is_active)
        return new CoulombPBCAATemp(qp,is_active);
      else
        return new CoulombPBCAATemp(*this);
    }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

