//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
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
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus {

  LRTwoBodyJastrow::LRTwoBodyJastrow(ParticleSet& p, HandlerType* inHandler):
  NumPtcls(0), NumSpecies(0), skRef(0) {
    Optimizable=false;
    handler=inHandler;
    NumSpecies=p.groups();
    skRef=p.SK;
    handler = LRJastrowSingleton::getHandler(p);
    if(skRef) {
      //Rs  = Omega
      OneOverCellVolume = 1.0 / p.Lattice.Volume;
      Omega = std::pow(3.0/(4.0*M_PI)*p.Lattice.Volume/static_cast<RealType>(p.getTotalNum()),1.0/3.0);
      //Rs=std::pow(3.0/4.0/M_PI*p.Lattice.Volume/static_cast<RealType>(p.getTotalNum()),1.0/3.0);
      //Omega=std::sqrt(4.0*M_PI*static_cast<RealType>(p.getTotalNum())/p.Lattice.Volume);
      OneOverOmega=1.0/Omega;
      FourPiOmega=4.0*M_PI*Omega;
      NumPtcls=p.getTotalNum();
      NumKpts=skRef->KLists.numk;
      NormConstant=FourPiOmega*NumPtcls*(NumPtcls-1)*0.5;
      resize();
      resetInternals();
    }
  }
  
  void LRTwoBodyJastrow::resize() {
    Rhok.resize(NumKpts);
    rokbyF.resize(NumPtcls,NumKpts);
    U.resize(NumPtcls);
    dU.resize(NumPtcls);
    d2U.resize(NumPtcls);
    FirstAddressOfdU=&(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU+NumPtcls*DIM;
    
    offU.resize(NumPtcls);
    offdU.resize(NumPtcls);
    offd2U.resize(NumPtcls);
  }
  
  
  /** update Fk using the handler
   */
  void LRTwoBodyJastrow::resetInternals() 
  {
    Fk_0.resize(handler->Fk.size());
    Fk_0 = -1.0 * handler->Fk;
    Fk.resize(Fk_0.size());
    Fk=Fk_0;
  }
  
  void LRTwoBodyJastrow::resetParameters(OptimizableSetType& optVariables) 
  {
    ///DO NOTHING FOR NOW
  }

  void LRTwoBodyJastrow::resetTargetParticleSet(ParticleSet& P) {
    // update handler as well, should there also be a reset?
    skRef=P.SK;
    handler->initBreakup(P);
    resetInternals();
  }
  
  LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::evaluateLog(ParticleSet& P, 
				      ParticleSet::ParticleGradient_t& G, 
				      ParticleSet::ParticleLaplacian_t& L) {
      
      Rhok=0.0;
      for(int spec1=0; spec1<NumSpecies; spec1++) {
        const ComplexType* restrict rhok(P.SK->rhok[spec1]);
        for(int ki=0; ki<MaxK; ki++) {
          Rhok[ki] += rhok[ki];
        }
      }
      
      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
      const KContainer::SContainer_t& ksq(P.SK->KLists.ksq);
      
      ValueType sum(0.0);
      for(int iat=0; iat<NumPtcls; iat++) {
        ValueType res(0.0),l(0.0);
        GradType g;
        const ComplexType* restrict eikr(P.SK->eikr[iat]);
        for(int ki=0; ki<MaxK; ki++) {
          ComplexType skp((Fk[ki]*conj(eikr[ki])*Rhok[ki]));
#if defined(QMC_COMPLEX)
          res +=  skp;
          l += ksq[ki]*(Fk[ki]-skp);
          g += ComplexType(skp.imag(),-skp.real())*kpts[ki];
#else
          res +=  skp.real();
          g += kpts[ki]*skp.imag();
          l += ksq[ki]*(Fk[ki]-skp.real());
#endif
        }
        sum+=(U[iat]=res);
        G[iat]+=(dU[iat]=g);
        L[iat]+=(d2U[iat]=l);
      }
      
      return sum*0.5;
    }
  
  
  LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::ratio(ParticleSet& P, int iat) {
      //restore, if called should do nothing
      NeedToRestore=false;
      curVal=0.0;

//      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
//      const Vector<ComplexType>& eikr1(P.SK->eikr_new);
//      const Vector<ComplexType>& del_eikr(P.SK->delta_eikr);
//      //Rhok += del_eikr;
//      for(int ki=0; ki<MaxK; ki++) {
//        //ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
//        ComplexType skp((Fk[ki]*conj(eikr1[ki])*(Rhok[ki]+del_eikr[ki])));
//#if defined(QMC_COMPLEX)
//        curVal +=  skp;
//#else
//        curVal +=  skp.real();
//#endif
//      }
      return std::exp(curVal-U[iat]);
    }
  
  
  LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::logRatio(ParticleSet& P, int iat,
				   ParticleSet::ParticleGradient_t& dG,
				   ParticleSet::ParticleLaplacian_t& dL) {
      
      NeedToRestore=true;
      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
      const KContainer::SContainer_t& ksq(P.SK->KLists.ksq);
      const Vector<ComplexType>& eikr1(P.SK->eikr_new);
      const Vector<ComplexType>& del_eikr(P.SK->delta_eikr);
      
      //add the difference
      Rhok += del_eikr;
      
      curVal=0.0; curLap=0.0; curGrad=0.0;
      for(int ki=0; ki<MaxK; ki++) {
#if defined(QMC_COMPLEX)
        ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
        curVal +=  skp;
        curGrad += ComplexType(skp.imag(),-skp.real())*kpts[ki];
        curLap += ksq[ki]*(Fk[ki]-skp);
#else
        RealType skp_r=Fk[ki]*(eikr1[ki].real()*Rhok[ki].real()+eikr1[ki].imag()*Rhok[ki].imag());
        RealType skp_i=Fk[ki]*(eikr1[ki].real()*Rhok[ki].imag()-eikr1[ki].imag()*Rhok[ki].real());
        curVal +=  skp_r;
        curLap += ksq[ki]*(Fk[ki]-skp_r);
        curGrad += kpts[ki]*skp_i;
#endif
      }

      for(int jat=0;jat<NumPtcls; jat++) 
      {
        if(jat == iat) continue;
        const ComplexType* restrict eikrj(P.SK->eikr[jat]);
        GradType g;
        ValueType l(0.0), v(0.0);
        for(int ki=0; ki<MaxK; ki++) {
#if defined(QMC_COMPLEX)
          ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
          GradType dg(skp.imag()*kpts[ki]);
          ValueType dl(skp.real()*ksq[ki]);
          v += skp.real();
          g +=dg;
          l -= dl;
#else
          ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
          //GradType dg(skp.imag()*kpts[ki]);
          //ValueType dl(skp.real()*ksq[ki]);
          v += skp.real();
          l -= skp.real()*ksq[ki];
          g += skp.imag()*kpts[ki];
#endif
          //dG[jat] += Fk[ki]*skp.imag()*kpts[ki];
          //dL[jat] -= Fk[ki]*skp.real()*ksq[ki];
        }
        offU[jat]=v;
        offdU[jat]=g;
        offd2U[jat]=l;
        dG[jat] += g;
        dL[jat] += l;
      }
//      for(int jat=0;jat<NumPtcls; jat++) {
//        if(iat==jat) {
//          for(int ki=0; ki<MaxK; ki++) {
//            //ComplexType rhok_new(Rhok[ki]+del_eikr[ki]);
//            //ComplexType skp((Fk[ki]*conj(eikr1[ki])*rhok_new));
//#if defined(QMC_COMPLEX)
//            ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
//            curVal +=  skp;
//            curGrad += ComplexType(skp.imag(),-skp.real())*kpts[ki];
//            curLap += ksq[ki]*(Fk[ki]-skp);
//#else
//            RealType skp_r=Fk[ki]*(eikr1[ki].real()*Rhok[ki].real()+eikr1[ki].imag()*Rhok[ki].imag());
//            RealType skp_i=Fk[ki]*(eikr1[ki].real()*Rhok[ki].imag()-eikr1[ki].imag()*Rhok[ki].real());
//            curVal +=  skp_r;
//            curLap += ksq[ki]*(Fk[ki]-skp_r);
//            curGrad += kpts[ki]*skp_i;
//            //curVal +=  skp.real();
//            //curLap += ksq[ki]*(Fk[ki]-skp.real());
//            //curGrad += skp.imag()*kpts[ki];
//#endif
//          }
//        } else {
//          const ComplexType* restrict eikrj(P.SK->eikr[jat]);
//          GradType g;
//          ValueType l(0.0), v(0.0);
//          for(int ki=0; ki<MaxK; ki++) {
//#if defined(QMC_COMPLEX)
//            ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
//            GradType dg(skp.imag()*kpts[ki]);
//            ValueType dl(skp.real()*ksq[ki]);
//            v += skp.real();
//            g +=dg;
//            l -= dl;
//#else
//            ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
//            //GradType dg(skp.imag()*kpts[ki]);
//            //ValueType dl(skp.real()*ksq[ki]);
//            v += skp.real();
//            l -= skp.real()*ksq[ki];
//            g += skp.imag()*kpts[ki];
//#endif
//            //dG[jat] += Fk[ki]*skp.imag()*kpts[ki];
//            //dL[jat] -= Fk[ki]*skp.real()*ksq[ki];
//          }
//          offU[jat]=v;
//          offdU[jat]=g;
//          offd2U[jat]=l;
//          dG[jat] += g;
//          dL[jat] += l;
//        }
//      }
//      
      dG[iat] += offdU[iat] = curGrad-dU[iat];
      dL[iat] += offd2U[iat] = curLap-d2U[iat];
      return offU[iat] = curVal-U[iat];
    }
  
  void LRTwoBodyJastrow::restore(int iat) {
    //substract the addition in logRatio
    if(NeedToRestore) Rhok -= skRef->delta_eikr;
  }

  void LRTwoBodyJastrow::acceptMove(ParticleSet& P, int iat) {
    U += offU;
    dU += offdU;
    d2U += offd2U;
  }

  void LRTwoBodyJastrow::update(ParticleSet& P, 
				    ParticleSet::ParticleGradient_t& dG, 
				    ParticleSet::ParticleLaplacian_t& dL,
				    int iat) {
    app_error() << "LRTwoBodyJastrow::update is INCOMPLETE " << endl;
  }
  
  
  LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::registerData(ParticleSet& P, PooledData<RealType>& buf) {
      LogValue=evaluateLog(P,P.G,P.L); 
      buf.add(Rhok.first_address(), Rhok.last_address());
      buf.add(U.first_address(), U.last_address());
      buf.add(d2U.first_address(), d2U.last_address());
      buf.add(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }

  LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      LogValue=evaluateLog(P,P.G,P.L); 
      buf.put(Rhok.first_address(), Rhok.last_address());
      buf.put(U.first_address(), U.last_address());
      buf.put(d2U.first_address(), d2U.last_address());
      buf.put(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }

  void LRTwoBodyJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    buf.get(Rhok.first_address(), Rhok.last_address());
    buf.get(U.first_address(), U.last_address());
    buf.get(d2U.first_address(), d2U.last_address());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
  }
  
  LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
      buf.put(Rhok.first_address(), Rhok.last_address());
      buf.put(U.first_address(), U.last_address());
      buf.put(d2U.first_address(), d2U.last_address());
      buf.put(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }
  
  
  bool
    LRTwoBodyJastrow::put(xmlNodePtr cur, VarRegistry<RealType>& vlist) {
      
      if(skRef == 0) {
        app_error() << "  LRTowBodyJastrow should not be used for non periodic systems." << endl;
        return false;
      }
      
      int nkshells=skRef->KLists.kshell.size()-1;
      Fk_symm.resize(nkshells);
      
      bool foundCoeff=false;
      if(cur != NULL)
      {
        xmlNodePtr tcur=cur->children;
        while(tcur != NULL) {
          string cname((const char*)(tcur->name));
          if(cname == "parameter") {
            const xmlChar* kptr=xmlGetProp(tcur,(const xmlChar *)"name");
            const xmlChar* idptr=xmlGetProp(tcur,(const xmlChar *)"id");
            if(idptr!= NULL && kptr != NULL) {
              int ik=atoi((const char*)kptr);
              if(ik<Fk_symm.size()) { // only accept valid ik 
                RealType x;
                putContent(x,tcur);
                Fk_symm[ik]=x;
                vlist[(const char*)idptr]=x;
                //vlist.add((const char*)idptr,Fk_symm.data()+ik);
              }
              foundCoeff=true;
            }
          }
          tcur=tcur->next;
        }
      }
      Fk.resize(NumKpts);
      if(foundCoeff) {
        resetInternals();
      } else {
        int ki=0; 
        char coeffname[128];
        MaxK=0;
        int ish=0;
        //RealType kc2=(skRef->KLists.Lattice.LR_kc)*(skRef->KLists.Lattice.LR_kc);
        RealType kc2=2.0;
        while(ish<nkshells && ki<NumKpts)
        {
          Fk_symm[ish]=-1.0*handler->Fk[ki];
          sprintf(coeffname,"rpa_k%d",ish);
          vlist[coeffname]=Fk_symm[ish];
          //vlist.add(coeffname,Fk_symm.data()+ik);
          std::ostringstream kname,val;
          kname << ish;
          val<<Fk_symm[ish];
          xmlNodePtr p_ptr = xmlNewTextChild(cur,NULL,(const xmlChar*)"parameter",
              (const xmlChar*)val.str().c_str());
          xmlNewProp(p_ptr,(const xmlChar*)"id",(const xmlChar*)coeffname);
          xmlNewProp(p_ptr,(const xmlChar*)"name",(const xmlChar*)kname.str().c_str());

          if(skRef->KLists.ksq[ki]<kc2) MaxK=skRef->KLists.kshell[ish+1];;
          for(; ki<skRef->KLists.kshell[ish+1]; ki++) Fk[ki]=Fk_symm[ish];

          if(abs(Fk_symm[ish])>1e-3) MaxK=ki;
          ++ish;
        }
      }
    
      app_log() << "  Long-range Two-Body Jastrow coefficients " << endl;
      app_log() << "  MaxK = " << MaxK << " NumKpts = " << NumKpts << endl;
      for(int ikpt=0; ikpt<NumKpts; ikpt++) {
        app_log() <<  skRef->KLists.ksq[ikpt] << " " << Fk[ikpt] << endl;
      }
      return true;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
