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
        NormConstant=FourPiOmega*NumPtcls*(NumPtcls-1)*0.5;
        resize();
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
      curVal=0.0;
      const KContainer::VContainer_t& kpts(P.SK->KLists.kpts_cart);
      const Vector<ComplexType>& eikr1(P.SK->eikr_new);
      const Vector<ComplexType>& del_eikr(P.SK->delta_eikr);
      Rhok += del_eikr;
      for(int ki=0; ki<NumKpts; ki++) {
        ComplexType skp((Fk[ki]*conj(eikr1[ki])*Rhok[ki]));
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
      Rhok -= skRef->delta_eikr;
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

    void 
    LRTwoBodyJastrow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
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

      //not optimization yet
      Fk.resize(NumKpts);
      app_log() << "  LRTwoBodyJastrow: kx ky kz U[k] " << endl;
      for(int ikpt=0; ikpt<NumKpts; ikpt++) {
        Fk[ikpt]= -0.001*getRPACoeff(skRef->KLists.ksq[ikpt]);
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
