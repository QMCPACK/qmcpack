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

    LRTwoBodyJastrow::LRTwoBodyJastrow(ParticleSet& p) {
      PtclRef=&p;
      NumSpecies=p.groups();
    }

    void LRTwoBodyJastrow::reset() {
    }

    void LRTwoBodyJastrow::resetTargetParticleSet(ParticleSet& P) {
    }

    LRTwoBodyJastrow::ValueType LRTwoBodyJastrow::evaluateLog(ParticleSet& P,
		         ParticleSet::ParticleGradient_t& G, 
		         ParticleSet::ParticleLaplacian_t& L) {

      Sk=0.0;
      for(int ki=0; ki<PtclRhoK.KLists.kpts_cart.size(); ki++) {
        int kj=PtclRhoK.KLists.minusk[ki];
        for(int spec1=0; spec1<NumSpecies; spec1++) {
          Sk[ki] += PtclRhoK.rhok(ki,spec1);
        }
      }

      ComplexType sum=0.0;
      for(int iat=0; iat<NumPtcls; iat++) {
        for(int ki=0; ki<PtclRhoK.KLists.kpts_cart.size(); ki++) {
          ComplexType skp=(Sk[ki]-SkP(iat,ki))*Fk[ki]*conj(SkP(iat,ki));
          sum += skp;
          g+=K[ki]*skp;
          l+=L[ki]*skp;
        }
        G[iat]+=g;
        L[iat]+=lap;
      }

      return ValueType(sum);//use trait
    }

    inline ValueType evaluate(ParticleSet& P,
			      ParticleSet::ParticleGradient_t& G, 
			      ParticleSet::ParticleLaplacian_t& L) {
      return exp(evaluateLog(P,G,L));
    }

    ValueType ratio(ParticleSet& P, int iat);

    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL) ;

    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL);

    void restore(int iat);
    void acceptMove(ParticleSet& P, int iat);
    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat);


    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);
    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);
  };
}
#endif
