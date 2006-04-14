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
#ifndef QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/DummyBasisSet.h"

namespace qmcplusplus {

  /** OrbitalBuilder for Slater determinants of electron-gas 
  */
  class ElectronGasComplexOrbitalBuilder: public OrbitalBuilderBase {
  public:

    ///constructor
    ElectronGasComplexOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);

    ///implement vritual function
    bool put(xmlNodePtr cur);

    struct EGOSet {

      typedef DummyBasisSet BasisSet_t;

      int KptMax;
      RealType kdotr;
      vector<PosType> K;
      vector<RealType> mK2;

      EGOSet(const vector<PosType>& k, const vector<RealType>& k2);

      inline void reset() { }
      inline void resetTargetParticleSet(ParticleSet& P) { }
      inline void resizeByWalkers(int nwalkers) {}

      inline ValueType
      evaluate(const ParticleSet& P, int iat, int jorb) {
        cout << "EGOSet::this should not be used" << endl;
        kdotr=dot(K[jorb],P.R[iat]);
        return ValueType(std::cost(kdotr),std::sin(kdotr));
      }

      template<class VV>
        inline void 
        evaluate(const ParticleSet& P, int iat, VV& psi) {
          for(int ik=0; ik<KptMax; ik++) {
            kdotr=dot(K[ik],P.R[iat]);
            psi[ik]=ValueType(std::cos(kdotr),std::sin(kdotr));
          }
        }

      template<class VV, class GV>
        inline void 
        evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi) {
          for(int ik=0; ik<KptMax; ik++) {
            kdotr=dot(K[ik],P.R[iat]);
            RealType coskr=std::cos(kdotr);
            RealType sinkr=std::sin(kdotr);
            psi[ik]=ValueType(coskr,sinkr);
            dpsi[ik]=ValueType(-sinkr,coskr)*K[ik];
            d2psi[ik]=ValueType(mK2[ik]*coskr,mK2[ik]*sinkr);
          }
        }

      template<class VM, class GM>
        inline void 
        evaluate(const ParticleSet& P, int first, int last,
            VM& logdet, GM& dlogdet, VM& d2logdet) {
          for(int i=0,iat=first; iat<last; i++,iat++) {
          for(int ik=0; ik<KptMax; ik++) {
            kdotr=dot(K[ik],P.R[iat]);
            RealType coskr=std::cos(kdotr);
            RealType sinkr=std::sin(kdotr);
            logdet(ik,i)=ValueType(coskr,sinkr);
            dlogdet(i,ik)=ValueType(-sinkr,coskr)*K[ik];
            d2logdet(i,ik)=ValueType(mK2[ik]*coskr,mK2[ik]*sinkr);
          }
        }
    };

  private:

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
