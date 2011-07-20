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
#ifndef QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/BasisSetBase.h"
//#include "QMCWaveFunctions/DummyBasisSet.h"

namespace qmcplusplus 
{

  struct EGOSet: public SPOSetBase
  {

    int KptMax;
    vector<PosType> K;
    vector<RealType> mK2;

    EGOSet(const vector<PosType>& k, const vector<RealType>& k2);

    SPOSetBase* makeClone() const
    {
      return new EGOSet(*this);
    }

    void resetParameters(const opt_variables_type& optVariables){}
    inline void resetTargetParticleSet(ParticleSet& P) { }
    void setOrbitalSetSize(int norbs) { }

    inline void 
      evaluate(const ParticleSet& P, int iat, ValueVector_t& psi) 
      {
        RealType sinkr,coskr;
        for(int ik=0; ik<KptMax; ik++) {
          sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
          psi[ik]=ValueType(coskr,sinkr);
        }
      }

    /** generic inline function to handle a row
     * @param r position of the particle
     * @param psi value row
     * @param dpsi gradient row
     * @param d2psi laplacian row
     */
    void evaluate_p(const PosType& r, ValueType* restrict psi, GradType* restrict dpsi
        , ValueType* restrict d2psi)
    {
      RealType sinkr,coskr;
      for(int ik=0; ik<KptMax; ik++) 
      {
        sincos(dot(K[ik],r),&sinkr,&coskr);
        psi[ik]  =ValueType(coskr,sinkr);
        dpsi[ik] =ValueType(-sinkr,coskr)*K[ik];
        d2psi[ik]=ValueType(mK2[ik]*coskr,mK2[ik]*sinkr);
      }
    }

    inline void 
      evaluate(const ParticleSet& P, int iat, 
          ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
      {
        evaluate_p(P.R[iat],psi.data(),dpsi.data(),d2psi.data());
      }

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi
        , HessVector_t& grad_grad_psi)
    {
      APP_ABORT("Incomplete implementation EGOSet::evaluate(P,iat,psi,dpsi,grad_grad_psi)");
    }
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      for(int i=0,iat=first; iat<last; i++,iat++) 
        evaluate_p(P.R[iat],logdet[i],dlogdet[i],d2logdet[i]);
    }
    void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
    {
      APP_ABORT("Incomplete implementation EGOSet::evaluate(P,first,last,lodget,dlodet,grad_grad_logdet,grad_grad_grad_logdet");
    }
  };


  /** OrbitalBuilder for Slater determinants of electron-gas 
  */
  class ElectronGasComplexOrbitalBuilder: public OrbitalBuilderBase {
  public:

    ///constructor
    ElectronGasComplexOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);
    //typedef VarRegistry<RealType> OptimizableSetType;

    ///implement vritual function
    bool put(xmlNodePtr cur);

  };

  /** OrbitalBuilder for Slater determinants of electron-gas 
  */
  class ElectronGasBasisBuilder: public BasisSetBuilder {
    protected:
      typedef map<string,ParticleSet*> PtclPoolType;
      typedef map<string,SPOSetBase*>  SPOPoolType;
      ParticleSet *targetPtcl;
      EGOSet* myBasis;
    public:
      ///constructor
      ElectronGasBasisBuilder(ParticleSet& p, xmlNodePtr cur);
      
      ///implement vritual function
      bool put(xmlNodePtr cur);
      /** initialize the Antisymmetric wave function for electrons
      *@param cur the current xml node
      */
      SPOSetBase* createSPOSet(xmlNodePtr cur) {return myBasis;}; 
      
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
