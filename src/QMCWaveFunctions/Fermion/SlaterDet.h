//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/DummyBasisSet.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"

namespace qmcplusplus {

  class SlaterDet: public OrbitalBase {

  public:

    typedef DiracDeterminantBase Determinant_t;

    /// constructor
    SlaterDet() {M.resize(3,0);Optimizable=false;}

    ///destructor
    ~SlaterDet();

    ///add a new DiracDeterminant to the list of determinants
    void add(Determinant_t* det);

    ///reset all the Dirac determinants, Optimizable is true
    void reset();

    void resetTargetParticleSet(ParticleSet& P);

    ValueType 
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L);

    ValueType 
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L);

    void resizeByWalkers(int nwalkers) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->resizeByWalkers(nwalkers);
    }

    ///return the total number of Dirac determinants
    inline int size() const { return Dets.size();}

    ///return the dimension of the i-th Dirac determinant
    inline int size(int i) const { return Dets[i]->cols();}

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);
    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    inline ValueType ratio(ParticleSet& P, int iat,
			   ParticleSet::ParticleGradient_t& dG, 
			   ParticleSet::ParticleLaplacian_t& dL) { 
      return Dets[DetID[iat]]->ratio(P,iat,dG,dL);
    }

    inline ValueType logRatio(ParticleSet& P, int iat,
			   ParticleSet::ParticleGradient_t& dG, 
			   ParticleSet::ParticleLaplacian_t& dL) { 
      ValueType r = Dets[DetID[iat]]->ratio(P,iat,dG,dL);
      SignValue = (r<0.0)?-1.0:1.0;
      return log(abs(r));
    }
    
    inline void restore(int iat) {
      return Dets[DetID[iat]]->restore(iat);
    }

    inline void update(ParticleSet& P, int iat) {
      Dets[DetID[iat]]->update(P,iat);
    }

    ValueType
    ratio(ParticleSet& P, int iat) {
      return Dets[DetID[iat]]->ratio(P,iat);
    } 	  

    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
      return Dets[DetID[iat]]->update(P,dG,dL,iat);
    }

  private:
    vector<int> M;
    vector<int> DetID;
    ///container for the DiracDeterminants
    vector<Determinant_t*>  Dets;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
