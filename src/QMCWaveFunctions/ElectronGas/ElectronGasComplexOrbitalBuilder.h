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
    typedef VarRegistry<RealType> OptimizableSetType;

    ///implement vritual function
    bool put(xmlNodePtr cur);

    struct EGOSet {

      typedef DummyBasisSet BasisSet_t;

      int KptMax;
      vector<PosType> K;
      vector<RealType> mK2;

      EGOSet(const vector<PosType>& k, const vector<RealType>& k2);

      EGOSet* makeClone() const
      {
        return new EGOSet(*this);
      }

      inline void resetParameters(OptimizableSetType& vlist) { }
      inline void resetTargetParticleSet(ParticleSet& P) { }

      inline ValueType
      evaluate(const ParticleSet& P, int iat, int jorb) {
        cout << "EGOSet::this should not be used" << endl;
        RealType kdotr=dot(K[jorb],P.R[iat]);
        return ValueType(std::cos(kdotr),std::sin(kdotr));
      }

      template<class VV>
        inline void 
        evaluate(const ParticleSet& P, int iat, VV& psi) {
          RealType sinkr,coskr;
          for(int ik=0; ik<KptMax; ik++) {
            //RealType kdotr=dot(K[ik],P.R[iat]);
            //psi[ik]=ValueType(std::cos(kdotr),std::sin(kdotr));
            sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
            psi[ik]=ValueType(coskr,sinkr);
          }
        }

      template<class VV, class GV>
        inline void 
        evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi) {
          RealType sinkr,coskr;
          for(int ik=0; ik<KptMax; ik++) {
            //RealType kdotr=dot(K[ik],P.R[iat]);
            //RealType coskr=std::cos(kdotr);
            //RealType sinkr=std::sin(kdotr);
            sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
            psi[ik]=ValueType(coskr,sinkr);
            dpsi[ik]=ValueType(-sinkr,coskr)*K[ik];
            d2psi[ik]=ValueType(mK2[ik]*coskr,mK2[ik]*sinkr);
          }
        }

      template<class VM, class GM>
        inline void 
        evaluate(const ParticleSet& P, int first, int last,
            VM& logdet, GM& dlogdet, VM& d2logdet) {
          RealType sinkr,coskr;
          for(int i=0,iat=first; iat<last; i++,iat++) {
            for(int ik=0; ik<KptMax; ik++) {
              //RealType kdotr=dot(K[ik],P.R[iat]);
              //RealType coskr=std::cos(kdotr);
              //RealType sinkr=std::sin(kdotr);
              sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
              logdet(ik,i)=ValueType(coskr,sinkr);
              dlogdet(i,ik)=ValueType(-sinkr,coskr)*K[ik];
              d2logdet(i,ik)=ValueType(mK2[ik]*coskr,mK2[ik]*sinkr);
            }
          }
        }
    };

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
