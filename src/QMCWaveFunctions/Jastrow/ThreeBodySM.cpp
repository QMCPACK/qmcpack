//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and John Gergely
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
/**@file ThreeBodySM.cpp
 * @brief Definition of ThreeBodySM.
 */
#include "QMCWaveFunctions/Jastrow/ThreeBodySM.h"

namespace qmcplusplus {

  ThreeBodySM::ThreeBodySM(ParticleSet& ions, ParticleSet& els): CenterRef(ions) {
     dist_ee = DistanceTable::add(els);
     dist_ie = DistanceTable::add(ions,els);
  }

  void ThreeBodySM::resetTargetParticleSet(ParticleSet& P) {
     dist_ee = DistanceTable::add(P);
     dist_ie = DistanceTable::add(CenterRef,P);
  }
  
  ThreeBodySM::ValueType ThreeBodySM::evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, 
      ParticleSet::ParticleLaplacian_t& L) {

    LogValue=0.0;
    RealType dudr, d2udr2;

    int nc(CenterRef.getTotalNum()), nptcl(P.getTotalNum());

    //first fill the matrix AA(i,j) where j is a composite index
    for(int I=0; I<nc; I++) {
      BasisType& a(*ieBasis[CenterRef.GroupID[I]]);
      int offset(0);
      for(int nn=dist_ie->M[I]; nn<dist_ie->M[I+1]; nn++) {
        RealType sep(dist_ie->r(nn));
        RealType rinv(dist_ie->rinv(nn));
        int i(dist_ie->J[nn]);
        int offset(ieBasisOffset[I]);
        for(int k=0; k<a.size(); k++,offset++) {
          AA(i,offset)=a[k]->evaluate(sep,dudr,d2udr2);
          dudr *= rinv;
          dAA(i,offset)=dudr*dist_ie->dr(nn);
          d2AA(i,offset)=d2udr2+2.0*dudr;
        }
      }
    }

    for(int i=0; i<nptcl; i++) {
      for(int nn=dist_ee->M[i]; nn<dist_ee->M[i]; nn++) {
        int j(dist_ee->J[nn]);
        RealType sep(dist_ee->r(nn));
        RealType rinv(dist_ee->rinv(nn));
        for(int m=0; m<eeBasis.size(); m++) {
          RealType psum=0,lapmi=0,lapmj=0;
          PosType grmi,grmj;
          for(int I=0; I<nc; I++) {
            const Matrix<RealType>& cblock(*C(m,CenterRef.GroupID[I]));
            int offsetI(ieBasisOffSet[I]);
            for(int k=0; k< ieBasisSize[I],kb=offsetI; k++,kb++) {
              RealType vall=0,valk=AA(i,kb);
              for(int l=0; l<ieBasisSize[I],lb=offsetI; l++,lb++) {
                vall += cblock(k,l)*AA(j,lb);
                grmj += valk*cblock(k,l)*dAA(j,lb);
                lapmj += valk*cblock(k,l)*d2AA(j,lb);
              }//l
              psum += valk*vall;
              grmi += dAA(i,kb)*vall;
              lampi += d2AA(i,kb)*vall;
            }//k
          }//I

          RealType bm =eeBasis[m]->evaluate(sep,dudr,d2udr2);
          dudr *= rinv;
          PosType dbm=dudr*dist_ee->dr(nn);
          RealType d2bm=d2udr2+2.0*dudr;

          LogValue += bm*psum;

          G[i] += bm*grmi-dbm*psum;
          G[j] += bm*grmj+dbm*psum;
          L[i] += b2bm*psum+bm*lapi;
          L[j] += b2bm*psum+bm*lapj;

        }
      }
    }
    return LogValue;

  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
