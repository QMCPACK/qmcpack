//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_ELECTRONGAS_ORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_ORBITALS_H

#include <QMCWaveFunctions/OrbitalBuilderBase.h>
#include <QMCWaveFunctions/SPOSetBase.h>
#include <config/stdlib/math.h>

namespace qmcplusplus {

  struct RealEGOSet: public SPOSetBase
  {

    int KptMax;
    RealType kdotr;
    vector<PosType> K;
    vector<RealType> mK2;

    RealEGOSet(const vector<PosType>& k, const vector<RealType>& k2);

    void resetParameters(const opt_variables_type& optVariables){}
    inline void resetTargetParticleSet(ParticleSet& P) { }
    void setOrbitalSetSize(int norbs) { }

    SPOSetBase* makeClone() const
    {
      return new RealEGOSet(*this);
    }

    inline ValueType f(const PosType& pos,int i)
    {
      if(i>0)
      {
        int ik=(i-1)/2;
        int even=(i-1)%2;
        kdotr=dot(K[ik],pos);
        if(even)
          return std::cos(kdotr);
        else
          return std::sin(kdotr);
      } 
      else
        return 1.0;
    }

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi) 
    {
      RealType sinkr,coskr;
      psi[0]=1.0;
      for(int ik=0, j=1; ik<KptMax; ik++) 
      {
        sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
        psi[j++]=coskr;
        psi[j++]=sinkr;
      }
    }

    /** generic inline function to handle a row
     * @param r position of the particle
     * @param psi value row
     * @param dpsi gradient row
     * @param d2psi laplacian row
     */
    inline void evaluate_p(const PosType& r, ValueType* restrict psi, GradType* restrict dpsi, ValueType* restrict d2psi)
    {
       psi[0]=1.0;
       dpsi[0]=0.0;
       d2psi[0]=0.0;
       RealType coskr, sinkr;
       for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2)
       {
         int j2=j1+1;
         sincos(dot(K[ik],r),&sinkr,&coskr);
         psi[j1]=coskr;
         psi[j2]=sinkr;
         dpsi[j1]=-sinkr*K[ik];
         dpsi[j2]= coskr*K[ik];
         d2psi[j1]=mK2[ik]*coskr;
         d2psi[j2]=mK2[ik]*sinkr;
       }
    }

    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      evaluate_p(P.R[iat],psi.data(),dpsi.data(),d2psi.data());
    }

    /*
    void evaluate(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      RealType coskr, sinkr;
      for(int i=0,iat=first; iat<last; i++,iat++) {
        logdet(0,i)=1.0;
        dlogdet(i,0)=0.0;
        d2logdet(i,0)=0.0;
        for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2) {
          sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
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
    */

    void evaluate_notranspose(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      for(int i=0,iat=first; iat<last; i++,iat++)
        evaluate_p(P.R[iat],logdet[i],dlogdet[i],d2logdet[i]);
    }

  };

  /** OrbitalBuilder for Slater determinants of electron-gas 
  */
  class ElectronGasOrbitalBuilder: public OrbitalBuilderBase {
  public:

    ///constructor
    ElectronGasOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);

    ///implement vritual function
    bool put(xmlNodePtr cur);

  private:

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
