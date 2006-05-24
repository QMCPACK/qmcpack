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

namespace qmcplusplus {

    LRTwoBodyJastrow::LRTwoBodyJastrow(ParticleSet& p):NeverInitialized(true) {
      NumSpecies=p.groups();
    }

    void LRTwoBodyJastrow::reset() {
    }

    void LRTwoBodyJastrow::resetTargetParticleSet(ParticleSet& P) {
    }

    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L) {

      if(NeverInitialized) {
        Rhok.resize(P.SK->rhok.cols());
        rokbyF.resize(P.SK->rhok.rows(),P.SK->rhok.cols());
        U.resize(P.getTotalNum());
        dU.resize(P.getTotalNum());
        d2U.resize(P.getTotalNum());
        NeverInitialized=false;
      }

      Rhok=0.0;
      for(int spec1=0; spec1<NumSpecies; spec1++) {
        const ComplexType* restrict rhok(P.SK->rhok[spec1]);
        for(int ki=0; ki<nk; ki++) {
          Rhok[ki] += rhok[ki];
        }
      }

      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
      const KContainer::SContainer_t& ksq(P.SK->KLists.ksq);

      ValueType sum(0.0);
      for(int iat=0; iat<NumPtcls; iat++) {
        ValueType res(0.0);
        GradType g;
        ValueType l(0.0);
        const ComplexType* restrict eikr(P.SK->eikr[iat]);
        for(int ki=0; ki<nk; ki++) {
          //ComplexType skp((Rhok[ki]-eikr[ki])*Fk[ki]*conj(eikr[ki]));
          ComplexType skp(rokbyF(iat,ki)=(Rhok[ki]-eikr[ki])*Fk[ki]);
          skp*=conj(eikr[ki]);
#if defined(QMC_COMPLEX)
          res += skp;
          g+=skp*kpts[ki];
          l+=skp*ksq[ki];
#else
          res += skp.real();
          g+=kpts[ki]*skp.real();
          l+=ksq[ki]*skp.real();
#endif
        }
        sum+=(U[iat]=res);
        G[iat]+=(dU[iat]=g);
        L[iat]+=(d2U[iat]=lap);
      }
      return sum;//use trait
    }


    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::ratio(ParticleSet& P, int iat) {
      curVal=0.0;
      const ComplexType* restrict eikr(P.SK->eikr[iat]);
      for(int ki=0; ki<nk; ki++) {
        ComplexType skp(rhokbyF[ki]*conj(eikr[ki]));
#if defined(QMC_COMPLEX)
        curVal += skp;
#else
        curVal += skp.real();
#endif
      }
      return exp(curVal-U[iat]);
    }


    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL) {
      const ComplexType* restrict eikr(P.SK->eikr[iat]);
      curVal=0.0; curLap=0.0; curGrad=0.0;
      for(int ki=0; ki<nk; ki++) {
        ComplexType skp(rhokbyF[ki]*conj(eikr[ki]));
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
      for(int ki=0; ki<nk; ki++) {
        t[ki]=(Rhok[ki]-eikrNew[ki]+eikrOld[ki])*Fk[ki];
      }
    }

    void LRTwoBodyJastrow::update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
      return LogValue;
    }


    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::registerData(ParticleSet& P, PooledData<RealType>& buf) {
      LogValue=evaluateLog(P,P.G,P.L); 
      buf.add(U.begin(), U.end());
      buf.add(d2U.begin(), d2U.end());
      buf.add(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }

    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    }

    void 
    LRTwoBodyJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    }

    LRTwoBodyJastrow::ValueType 
    LRTwoBodyJastrow::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
      LogValue=evaluateLog(P,P.G,P.L); 
      buf.put(U.begin(), U.end());
      buf.put(d2U.begin(), d2U.end());
      buf.put(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
