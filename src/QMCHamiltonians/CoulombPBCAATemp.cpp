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
#include <numeric>

namespace qmcplusplus {

  CoulombPBCAATemp::CoulombPBCAATemp(ParticleSet& ref, bool active,
				     bool computeForces) : 
    AA(0), myGrid(0), rVs(0), 
    is_active(active), FirstTime(true), myConst(0.0),
    ComputeForces(computeForces), ForceBase(ref,ref)
  {
    ReportEngine PRE("CoulombPBCAATemp","CoulombPBCAATemp");

    //create a distance table: just to get the table name
    DistanceTableData *d_aa = DistanceTable::add(ref);
    PtclRefName=d_aa->Name;
    initBreakup(ref);
    prefix="F_AA";
    
    app_log() << "  Maximum K shell " << AA->MaxKshell << endl;
    app_log() << "  Number of k vectors " << AA->Fk.size() << endl;
    
    if(!is_active)
    {
      d_aa->evaluate(ref);
      RealType eL(0.0), eS(0.0);
      if (computeForces) {
	forces = 0.0;
	eS=evalSRwithForces(ref);
	// 1.3978248322
        eL=evalLRwithForces(ref);
	// 2.130267378 
      }
      else {
	eL=evalLR(ref);
	eS=evalSR(ref);
      }
      NewValue=Value = eL+eS+myConst;
      app_log() << "  Fixed Coulomb potential for " << ref.getName();
//       app_log() << "\n    V(short) =" << eS 
//                 << "\n    V(long)  =" << eL
      app_log() << "\n    e-e Madelung Const. =" << MC0
                << "\n    Vtot     =" << Value << endl;
    }
  }

  CoulombPBCAATemp:: ~CoulombPBCAATemp() { }

  void CoulombPBCAATemp::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    addValue(plist);
    if (ComputeForces) addObservablesF(plist); 
  }

  void CoulombPBCAATemp::resetTargetParticleSet(ParticleSet& P) {
    if(is_active)
    {
      PtclRefName=P.DistTables[0]->Name;
      AA->resetTargetParticleSet(P);
    }
  }

  CoulombPBCAATemp::Return_t 
    CoulombPBCAATemp::evaluate(ParticleSet& P) 
    {
      if(is_active) Value = evalLR(P)+evalSR(P)+myConst;
      return Value;
    }

  CoulombPBCAATemp::Return_t 
    CoulombPBCAATemp::registerData(ParticleSet& P, BufferType& buffer) 
    {
      if(is_active)
      {
        P.SK->DoUpdate=true;
        SR2.resize(NumCenters,NumCenters);
        dSR.resize(NumCenters);
        del_eikr.resize(P.SK->KLists.numk);

        Value=evaluateForPbyP(P);

        buffer.add(SR2.begin(),SR2.end());
        buffer.add(Value);
      }
      return Value;
    }

  CoulombPBCAATemp::Return_t 
    CoulombPBCAATemp::updateBuffer(ParticleSet& P, BufferType& buffer) 
    {
      if(is_active)
      {
        Value=evaluateForPbyP(P);

        buffer.put(SR2.begin(),SR2.end());
        buffer.put(Value);
      }
      return Value;
    }

  void CoulombPBCAATemp::copyFromBuffer(ParticleSet& P, BufferType& buffer) 
  {
    if(is_active)
    {
      buffer.get(SR2.begin(),SR2.end());
      buffer.get(Value);
    }
  }

  void CoulombPBCAATemp::copyToBuffer(ParticleSet& P, BufferType& buffer) 
  {
    if(is_active)
    {
      buffer.put(SR2.begin(),SR2.end());
      buffer.put(Value);
    }
  }

  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evaluateForPbyP(ParticleSet& P)
    {
      if(is_active)
      {
        SR2=0.0;
        Return_t res=myConst;
        const DistanceTableData *d_aa = P.DistTables[0];
        for(int iat=0; iat<NumCenters; iat++)
        {
          Return_t z=0.5*Zat[iat];
          for(int nn=d_aa->M[iat],jat=iat+1; nn<d_aa->M[iat+1]; ++nn,++jat) 
          {
            Return_t e=z*Zat[jat]*d_aa->rinv(nn)*rVs->splint(d_aa->r(nn));
            SR2(iat,jat)=e;
            SR2(jat,iat)=e;
            res+=e+e;
          }
        }
        return res+evalLR(P);
      }
      else
        return Value;
    }

  CoulombPBCAATemp::Return_t 
    CoulombPBCAATemp::evaluatePbyP(ParticleSet& P, int active)
    {
      if(is_active)
      {
        const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[0]->Temp);
        Return_t z=0.5*Zat[active];
        Return_t sr=0;
        const Return_t* restrict sr_ptr=SR2[active];
        for(int iat=0; iat<NumCenters; ++iat,++sr_ptr) 
        {
          if(iat==active) 
            dSR[active]=0.0;
          else
            sr+=dSR[iat]=(z*Zat[iat]*temp[iat].rinv1*rVs->splint(temp[iat].r1)- (*sr_ptr));
        }

        const StructFact& PtclRhoK(*(P.SK));
        const ComplexType* restrict eikr_new=PtclRhoK.eikr_temp.data();
        const ComplexType* restrict eikr_old=PtclRhoK.eikr[active];
        ComplexType* restrict d_ptr=del_eikr.data();
        for(int k=0; k<del_eikr.size(); ++k) *d_ptr++ = (*eikr_new++ - *eikr_old++);
        int spec2=SpeciesID[active];
        for(int spec1=0; spec1<NumSpecies; ++spec1)
        {
          Return_t zz=z*Zspec[spec1];
          sr += zz*AA->evaluate(PtclRhoK.KLists.kshell,PtclRhoK.rhok[spec1],del_eikr.data());
        }
        sr+= z*z*AA->evaluate(PtclRhoK.KLists.kshell,del_eikr.data(),del_eikr.data());
        //// const StructFact& PtclRhoK(*(PtclRef->SK));
        //const StructFact& PtclRhoK(*(P.SK));
        //const ComplexType* restrict eikr_new=PtclRhoK.eikr_temp.data();
        //const ComplexType* restrict eikr_old=PtclRhoK.eikr[active];
        //ComplexType* restrict d_ptr=del_eikr.data();
        //for(int k=0; k<del_eikr.size(); ++k) *d_ptr++ = (*eikr_new++ - *eikr_old++);

        //for(int iat=0;iat<NumCenters; ++iat)
        //{
        //  if(iat!=active)
        //    sr += z*Zat[iat]*AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.eikr[iat],del_eikr.data());
        //}
        return NewValue=Value+2.0*sr;
      }
      else
        return Value;
    }

  void CoulombPBCAATemp::acceptMove(int active)
  {
    if(is_active)
    {
      Return_t* restrict sr_ptr=SR2[active];
      Return_t* restrict pr_ptr=SR2.data()+active;
      for(int iat=0; iat<NumCenters; ++iat, ++sr_ptr,pr_ptr+=NumCenters)
        *pr_ptr = *sr_ptr += dSR[iat];
      Value=NewValue;
    }
  }

  void CoulombPBCAATemp::initBreakup(ParticleSet& P) 
  {
    //SpeciesSet& tspecies(PtclRef->getSpeciesSet());
    SpeciesSet& tspecies(P.getSpeciesSet());
    //Things that don't change with lattice are done here instead of InitBreakup()
    ChargeAttribIndx = tspecies.addAttribute("charge");
    MemberAttribIndx = tspecies.addAttribute("membersize");
    NumCenters = P.getTotalNum();
    NumSpecies = tspecies.TotalNum;

    Zat.resize(NumCenters);
    Zspec.resize(NumSpecies);
    NofSpecies.resize(NumSpecies);
    for(int spec=0; spec<NumSpecies; spec++) {
      Zspec[spec] = tspecies(ChargeAttribIndx,spec);
      NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx,spec));
    }

    SpeciesID.resize(NumCenters);
    for(int iat=0; iat<NumCenters; iat++)
    {
      SpeciesID[iat]=P.GroupID[iat];
      Zat[iat] = Zspec[P.GroupID[iat]];
    }

    AA = LRCoulombSingleton::getHandler(P);
    //AA->initBreakup(*PtclRef);
    myConst=evalConsts();
    myRcut=AA->get_rc();//Basis.get_rc();
    if(rVs==0) {
      rVs = LRCoulombSingleton::createSpline4RbyVs(AA,myRcut,myGrid);
    }
  }

  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalSR(ParticleSet& P) {
      const DistanceTableData *d_aa = P.DistTables[0];
      RealType SR=0.0;
      for(int ipart=0; ipart<NumCenters; ipart++){
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
    CoulombPBCAATemp::evalLR(ParticleSet& P) {
      RealType LR=0.0;
      const StructFact& PtclRhoK(*(P.SK));
      for(int spec1=0; spec1<NumSpecies; spec1++) {
        RealType Z1 = Zspec[spec1];
        for(int spec2=spec1; spec2<NumSpecies; spec2++) {
          RealType Z2 = Zspec[spec2];
          //RealType temp=AA->evaluate(PtclRhoK.KLists.minusk, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
          RealType temp=AA->evaluate(PtclRhoK.KLists.kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
          if(spec2==spec1)
            temp*=0.5;
          LR += Z1*Z2*temp;
        } //spec2
      }//spec1
      //LR*=0.5;
      return LR;
    }


  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalLRwithForces(ParticleSet& P) {
    RealType LR=0.0;
      const StructFact& PtclRhoK(*(P.SK));
      vector<TinyVector<RealType,DIM> > grad(P.getTotalNum());
      for(int spec2=0; spec2<NumSpecies; spec2++) {
        RealType Z2 = Zspec[spec2];
	for (int iat=0; iat<grad.size(); iat++) 
	  grad[iat] = TinyVector<RealType,DIM>(0.0);
	AA->evaluateGrad(P, P, spec2, Zat, grad);
	for (int iat=0; iat<grad.size(); iat++) 
	  forces[iat] += Z2*grad[iat];
      } //spec2
      return evalLR(P);
    }



  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalSRwithForces(ParticleSet& P) {
      const DistanceTableData *d_aa = P.DistTables[0];
      RealType SR=0.0;
      for(int ipart=0; ipart<NumCenters; ipart++){
        RealType esum = 0.0;
        for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++) {
	  RealType rV, d_rV_dr, d2_rV_dr2;
	  rV = rVs->splint(d_aa->r(nn), d_rV_dr, d2_rV_dr2);
	  RealType V = rV *d_aa->rinv(nn);
          esum += Zat[jpart]*d_aa->rinv(nn)*rV;
	  PosType grad = Zat[jpart]*Zat[ipart]*
	    (d_rV_dr - V)*d_aa->rinv(nn)*d_aa->rinv(nn)*d_aa->dr(nn);
	  forces[ipart] += grad;
	  forces[jpart] -= grad;
        }
        //Accumulate pair sums...species charge for atom i.
        SR += Zat[ipart]*esum;
      }
      return SR;
    }


  /** evaluate the constant term that does not depend on the position
   *
   * \htmlonly
   * <ul>
   * <li> self-energy: \f$ -\frac{1}{2}\sum_{i} v_l(r=0) q_i^2 = -\frac{1}{2}v_l(r=0) \sum_{alpha} N^{\alpha} q^{\alpha}^2\f$
   * <li> background term \f$ V_{bg} = -\frac{1}{2}\sum_{\alpha}\sum_{\beta} N^{\alpha}q^{\alpha}N^{\beta}q^{\beta} v_s(k=0)\f$
   * </ul>
   * \endhtmlonly
   * CoulombPBCABTemp contributes additional background term which completes the background term
   */
  CoulombPBCAATemp::Return_t
    CoulombPBCAATemp::evalConsts() {

      //LRHandlerType::BreakupBasisType &Basis(AA->Basis);
      //const Vector<RealType> &coefs(AA->coefs);
      RealType Consts=0.0; // constant term

      //v_l(r=0)
      RealType vl_r0 = AA->evaluateLR_r0();
      for(int spec=0; spec<NumSpecies; spec++) 
      {
        RealType z = Zspec[spec];
        RealType n = NofSpecies[spec];
        Consts -= 0.5*vl_r0*z*z*n;
      }

      //Compute Madelung constant: this is not correct for general cases
      MC0 = 0.0;
      for(int i=0;i<AA->Fk.size();i++) MC0 += AA->Fk[i];
      MC0 = 0.5*(MC0 - vl_r0);
      
      //Neutraling background term
      RealType vs_k0=AA->evaluateSR_k0(); //v_s(k=0)
      for(int speca=0;speca<NumSpecies;speca++) 
      {
        RealType za = Zspec[speca];
        RealType na = NofSpecies[speca];
        Consts -= 0.5*vs_k0*za*na*za*na;
        for(int specb=speca+1;specb<NumSpecies;specb++) {
          RealType zb = Zspec[specb];
          int nb = NofSpecies[specb];
          Consts -= vs_k0*za*zb*na*nb;
        }
      }

      //app_log() << "   Constant of PBCAA " << Consts << endl;
      //app_log() << "   MC0 of PBCAA " << MC0 << endl;
      return Consts;
    }

    QMCHamiltonianBase* CoulombPBCAATemp::makeClone(ParticleSet& qp, TrialWaveFunction& psi) 
    {
      if(is_active)
        return new CoulombPBCAATemp(qp,is_active,ComputeForces);
      else
        return new CoulombPBCAATemp(*this);//nothing needs to be re-evaluated
    }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

