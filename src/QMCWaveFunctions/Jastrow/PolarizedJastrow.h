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
#ifndef QMCPLUSPLUS_POLARIZED_ONEBODYJASTROW_H
#define QMCPLUSPLUS_POLARIZED_ONEBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "OhmmsData/ParameterSet.h"
//#include "OhmmsData/libxmldefs.h"
//#include <libxml++/libxml++.h>

namespace qmcplusplus {

  //  class OneBodyJastrow: public OrbitalBase {
  class PolarizedJastrow: public OrbitalBase {

  public:

    ParameterSet m_param;
    RealType alpha;

    ///constructor
    PolarizedJastrow():alpha(0.0){ 
      m_param.add(alpha,"alpha","none");
    }

    ~PolarizedJastrow(){
      DEBUGMSG("PolarizedJastrow::~PolarizedJastrow")
      //for(int i=0; i<F.size(); i++) delete F[i];
    }

    void resetParameters(OptimizableSetType& optVariables) 
    {
    }

    void resetTargetParticleSet(ParticleSet& P) {}

    void put(xmlNodePtr cur, VarRegistry<RealType>& vlist){
      m_param.put(cur);
      vlist.add("C_alpha",&alpha,1);
    }

    ValueType evaluateLog(ParticleSet& P,
		         ParticleSet::ParticleGradient_t& G, 
		         ParticleSet::ParticleLaplacian_t& L) {
      LogValue=0.0;
      for(int i=0; i<P.getTotalNum(); i++) {
	RealType p = 1.0+alpha*P.R[i][2];
	LogValue+=log(p);
	p=alpha/p;
	G[i][2] += p;
	L[i] -= p*p;
      }
      return LogValue;
    }
    ValueType evaluate(ParticleSet& P,
		       ParticleSet::ParticleGradient_t& G, 
		       ParticleSet::ParticleLaplacian_t& L) {
      return exp(evaluateLog(P,G,L));
    }

#ifdef USE_FASTWALKER
    void evaluate(WalkerSetRef& W, 
		  ValueVectorType& psi,
		  WalkerSetRef::WalkerGradient_t& G,
		  WalkerSetRef::WalkerLaplacian_t& L) {
    }
#else
    void evaluate(WalkerSetRef& W,
                  ValueVectorType& psi,
                  WalkerSetRef::WalkerGradient_t& G,
                  WalkerSetRef::WalkerLaplacian_t& L) {
    }
#endif


    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {
      std::cerr << "PolarizedJastrow::ratio for particle-by-particle is empty " << std::endl;
      return 1.0;
    }

    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {
      std::cerr << "PolarizedJastrow::logRatio for particle-by-particle is empty " << std::endl;
      return 0.0;
    }

    inline void restore(int iat) { }

    //@todo implement the virutal functions for particle-by-particle move
    void acceptMove(ParticleSet& P, int iat) {
      std::cerr << "PolarizedJastrow::update for particle-by-particle is empty " << std::endl;
    }

    ValueType ratio(ParticleSet& P, int iat) {
      std::cerr << "PolarizedJastrow::ratio for particle-by-particle is empty " << std::endl;
      return 1.0;
    }

    //@todo implement the virutal functions for particle-by-particle move
    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
      std::cerr << "PolarizedJastrow::update for particle-by-particle is empty " << std::endl;
    }


    ValueType evaluate(ParticleSet& P,PooledData<RealType>& buf) {
      std::cerr << "PolarizedJastrow::evaluate for particle-by-particle is empty " << std::endl;
      return 1.0;
    }

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf){
      std::cerr << "PolarizedJastrow::registerData for particle-by-particle is empty " << std::endl;
      return 0.0;
    }

    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf){
      std::cerr << "PolarizedJastrow::updateBuffer for particle-by-particle is empty " << std::endl;
      return 0.0;
    }

    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      std::cerr << "PolarizedJastrow::copyFromBuffer for particle-by-particle is empty " << std::endl;
    }

  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

