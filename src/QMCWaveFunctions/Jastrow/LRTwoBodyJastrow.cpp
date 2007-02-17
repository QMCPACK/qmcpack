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
    Optimizable=true;
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
  
  
  /**
   * update Fk using the handler
   */
  void LRTwoBodyJastrow::resetInternals() 
  {
    Fk.resize(handler->Fk.size());
    Fk = -1.0 * handler->Fk;
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
        for(int ki=0; ki<NumKpts; ki++) {
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
        for(int ki=0; ki<NumKpts; ki++) {
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
      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
      const Vector<ComplexType>& eikr1(P.SK->eikr_new);
      const Vector<ComplexType>& del_eikr(P.SK->delta_eikr);
      //Rhok += del_eikr;
      for(int ki=0; ki<NumKpts; ki++) {
        //ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
        ComplexType skp((Fk[ki]*conj(eikr1[ki])*(Rhok[ki]+del_eikr[ki])));
#if defined(QMC_COMPLEX)
        curVal +=  skp;
#else
        curVal +=  skp.real();
#endif
      }
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
      for(int jat=0;jat<NumPtcls; jat++) {
        if(iat==jat) {
          for(int ki=0; ki<NumKpts; ki++) {
            //ComplexType rhok_new(Rhok[ki]+del_eikr[ki]);
            //ComplexType skp((Fk[ki]*conj(eikr1[ki])*rhok_new));
            ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
#if defined(QMC_COMPLEX)
            curVal +=  skp;
            curGrad += ComplexType(skp.imag(),-skp.real())*kpts[ki];
            curLap += ksq[ki]*(Fk[ki]-skp);
#else
            curVal +=  skp.real();
            curGrad += kpts[ki]*skp.imag();
            curLap += ksq[ki]*(Fk[ki]-skp.real());
#endif
          }
        } else {
          const ComplexType* restrict eikrj(P.SK->eikr[jat]);
          GradType g;
          ValueType l(0.0), v(0.0);
          for(int ki=0; ki<NumKpts; ki++) {
            ComplexType skp(Fk[ki]*del_eikr[ki]*conj(eikrj[ki]));
            GradType dg(skp.imag()*kpts[ki]);
            ValueType dl(skp.real()*ksq[ki]);
            v += skp.real();
            g +=dg;
            l -= dl;
            //dG[jat] += Fk[ki]*skp.imag()*kpts[ki];
            //dL[jat] -= Fk[ki]*skp.real()*ksq[ki];
          }
          offU[jat]=v;
          offdU[jat]=g;
          offd2U[jat]=l;
          dG[jat] += g;
          dL[jat] += l;
        }
      }
      
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
      
      std::map<int,std::vector<int>*>& kpts_sorted(skRef->KLists.kpts_sorted);
      Fk_symm.resize(kpts_sorted.size());
      
      bool foundCoeff=false;
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
      
      Fk.resize(NumKpts);
      if(foundCoeff) {
        resetInternals();
      } else {
        std::map<int,std::vector<int>*>::iterator it(kpts_sorted.begin());
        int uniqueK=0;
        while(it != kpts_sorted.end()) {
          std::vector<int>::iterator vit((*it).second->begin());
          int ik=(*vit);
          Fk_symm[uniqueK]=Fk[ik]=-1.0*handler->Fk[ik];
          ++vit;
          while(vit != (*it).second->end()) {
            int ik=(*vit);
            Fk[ik]=-1.0*handler->Fk[ik];
            ++vit;
          }
          ++it;++uniqueK;
        }
        char coeffname[128];
        for(int ik=0; ik<Fk_symm.size(); ik++) {
          sprintf(coeffname,"rpa_k%d",ik);
          vlist[coeffname]=Fk_symm[ik];
          //vlist.add(coeffname,Fk_symm.data()+ik);
	  
          std::ostringstream kname,val;
          kname << ik;
          val<<Fk_symm[ik];
          xmlNodePtr p_ptr = xmlNewTextChild(cur,NULL,(const xmlChar*)"parameter",
					     (const xmlChar*)val.str().c_str());
          xmlNewProp(p_ptr,(const xmlChar*)"id",(const xmlChar*)coeffname);
          xmlNewProp(p_ptr,(const xmlChar*)"name",(const xmlChar*)kname.str().c_str());
        }
      }
      
      app_log() << "  Long-range Two-Body Jastrow coefficients " << endl;
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
