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
#include "QMCWaveFunctions/LRTwoBodyJastrow.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus {

    LRTwoBodyJastrow::LRTwoBodyJastrow(ParticleSet& p):
      NumPtcls(0), NumSpecies(0), skRef(0) {
      NumSpecies=p.groups();
      skRef=p.SK;
      if(skRef) {
        Rs=std::pow(3.0/4.0/M_PI*p.Lattice.Volume/static_cast<RealType>(p.getTotalNum()),1.0/3.0);
        Omega=std::sqrt(4.0*M_PI*static_cast<RealType>(p.getTotalNum())/p.Lattice.Volume);
        OneOverOmega=1.0/Omega;
        FourPiOmega=4.0*M_PI*Omega;
        NumPtcls=p.getTotalNum();
        NumKpts=skRef->KLists.numk;
        resize();
      }
    }

    void LRTwoBodyJastrow::resize() {
      Rhok.resize(NumKpts);
      rokbyF.resize(NumSpecies,NumKpts);
      U.resize(NumPtcls);
      dU.resize(NumPtcls);
      d2U.resize(NumPtcls);
    }

    void LRTwoBodyJastrow::reset() {
    }

    void LRTwoBodyJastrow::resetTargetParticleSet(ParticleSet& P) {
      skRef=P.SK;
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
          l += ksq[ki]*(Fk[ki]*Rhok[ki]-skp);
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
      cout<< " Sum " << sum << endl;
      return sum;//use trait
    }


    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::ratio(ParticleSet& P, int iat) {
      curVal=0.0;
      const ComplexType* restrict eikr(P.SK->eikr[iat]);
      const ComplexType* restrict rtemp(rokbyF[iat]);
      for(int ki=0; ki<NumKpts; ki++) {
        ComplexType skp(rtemp[ki]*conj(eikr[ki]));
#if defined(QMC_COMPLEX)
        curVal += skp;
#else
        curVal += skp.real();
#endif
      }
      return std::exp(curVal-U[iat]);
    }


    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL) {
      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
      const KContainer::SContainer_t& ksq(P.SK->KLists.ksq);

      const ComplexType* restrict eikr(P.SK->eikr[iat]);
      curVal=0.0; curLap=0.0; curGrad=0.0;
      const ComplexType* restrict rtemp(rokbyF[iat]);
      for(int ki=0; ki<NumKpts; ki++) {
        ComplexType skp(rtemp[ki]*conj(eikr[ki]));
#if defined(QMC_COMPLEX)
        curVal += skp;
        curGrad+=skp*kpts[ki];
        curLap+=skp*ksq[ki];
#else
        curVal += skp.real();
        curGrad+=kpts[ki]*skp.real();
        curLap +=ksq[ki]*skp.real();
#endif
      }

      dG[iat] += curGrad-dU[iat];
      dL[iat] += curLap-d2U[iat];
      return curVal-U[iat];
    }

    void LRTwoBodyJastrow::restore(int iat) {
      //do nothing
    }

    void LRTwoBodyJastrow::acceptMove(ParticleSet& P, int iat) {
      U[iat]=curVal;
      dU[iat]=curGrad;
      d2U[iat]=curLap;
      //update rhokbyF
      const ComplexType* restrict eikrOld(P.SKOld->eikr[iat]);
      const ComplexType* restrict eikrNew(P.SK->eikr[iat]);
      ComplexType* restrict t(rokbyF[iat]);
      for(int ki=0; ki<NumKpts; ki++) {
        t[ki]=(Rhok[ki]-eikrNew[ki]+eikrOld[ki])*Fk[ki];
      }
    }

    void LRTwoBodyJastrow::update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
    }


    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::registerData(ParticleSet& P, PooledData<RealType>& buf) {
      LogValue=evaluateLog(P,P.G,P.L); 

      //buf.add(U.begin(), U.end());
      //buf.add(d2U.begin(), d2U.end());
      //buf.add(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }

    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      return 1.0;
    }

    void 
    LRTwoBodyJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    }

    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
      LogValue=evaluateLog(P,P.G,P.L); 
      //buf.put(U.begin(), U.end());
      //buf.put(d2U.begin(), d2U.end());
      //buf.put(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }


    bool
    LRTwoBodyJastrow::put(xmlNodePtr cur, VarRegistry<RealType>& vlist) {

      if(skRef == 0) {
        app_error() << "  LRTowBodyJastrow should not be used for non periodic systems." << endl;
        return false;
      }

      //not optimization yet
      Fk.resize(NumKpts);
      app_log() << "  LRTwoBodyJastrow: kx ky kz U[k] " << endl;
      for(int ikpt=0; ikpt<NumKpts; ikpt++) {
        Fk[ikpt]= getRPACoeff(skRef->KLists.ksq[ikpt]);
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
