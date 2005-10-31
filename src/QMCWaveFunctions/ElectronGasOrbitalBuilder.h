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
#ifndef QMCPLUSPLUS_ELECTRONGAS_ORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_ORBITALS_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/DummyBasisSet.h"
namespace qmcplusplus {

  /** OrbitalBuilder for Slater determinants of electron-gas 
  */
  class ElectronGasOrbitalBuilder: public OrbitalBuilderBase {
  public:

    ///constructor
    ElectronGasOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);

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
        if(jorb) {
          if(jorb&1) {
            kdotr=dot(K[jorb/2],P.R[iat]);
            return cos(kdotr);
          }
          else {
            return sin(kdotr);
          }
        } else {
          return 1.0;
        }
      }

      template<class VV>
        inline void 
        evaluate(const ParticleSet& P, int iat, VV& psi) {
          psi[0]=1.0;
          for(int ik=0, j=1; ik<KptMax; ik++) {
            kdotr=dot(K[ik],P.R[iat]);
            psi[j++]=cos(kdotr);
            psi[j++]=sin(kdotr);
          }
        }

      template<class VV, class GV>
        inline void 
        evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi) {
          psi[0]=1.0;
          dpsi[0]=0.0;
          d2psi[0]=0.0;
          for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2) {
            int j2=j1+1;
            kdotr=dot(K[ik],P.R[iat]);
            RealType coskr=cos(kdotr);
            RealType sinkr=sin(kdotr);
            psi[j1]=coskr;
            psi[j2]=sinkr;
            dpsi[j1]=-sinkr*K[ik];
            dpsi[j2]= coskr*K[ik];
            d2psi[j1]=mK2[ik]*coskr;
            d2psi[j2]=mK2[ik]*sinkr;
          }
        }

      template<class VM, class GM>
        inline void 
        evaluate(const ParticleSet& P, int first, int last,
            VM& logdet, GM& dlogdet, VM& d2logdet) {
          for(int i=0,iat=first; iat<last; i++,iat++) {
            logdet(0,i)=1.0;
            dlogdet(i,0)=0.0;
            d2logdet(i,0)=0.0;
            for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2) {
              kdotr=dot(K[ik],P.R[iat]);
              RealType coskr=cos(kdotr);
              RealType sinkr=sin(kdotr);
              int j2=j1+1;
              logdet(j1,i)=coskr;
              logdet(j2,i)=sinkr;
              dlogdet(i,j1)=-sinkr*K[ik];
              dlogdet(i,j2)= coskr*K[ik];
              d2logdet(i,j1)=mK2[ik]*coskr;
              d2logdet(i,j2)=mK2[ik]*sinkr;
            }
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
