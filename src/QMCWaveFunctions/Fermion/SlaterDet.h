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
#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"

namespace qmcplusplus {

  class SlaterDet: public OrbitalBase {

  public:

    typedef DiracDeterminantBase Determinant_t;

    /// constructor
    SlaterDet();

    ///destructor
    ~SlaterDet();

    ///add a new DiracDeterminant to the list of determinants
    void add(Determinant_t* det);

    void checkInVariables(opt_variables_type& active);

    void checkOutVariables(const opt_variables_type& active);

    ///reset all the Dirac determinants, Optimizable is true
    void resetParameters(const opt_variables_type& optVariables);

    void reportStatus(ostream& os);

    void resetTargetParticleSet(ParticleSet& P);

    ValueType 
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L);

    RealType 
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L);

    ///return the total number of Dirac determinants
    inline int size() const { return Dets.size();}

    ///return the column dimension of the i-th Dirac determinant
    inline int size(int i) const { return Dets[i]->cols();}

    RealType registerData(ParticleSet& P, PooledData<RealType>& buf);
    RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf);
    void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

    inline ValueType ratio(ParticleSet& P, int iat,
			   ParticleSet::ParticleGradient_t& dG, 
			   ParticleSet::ParticleLaplacian_t& dL) { 
      return Dets[DetID[iat]]->ratio(P,iat,dG,dL);
    }

    inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      return Dets[DetID[iat]]->ratioGrad(P,iat,grad_iat);
    }

    GradType evalGrad(ParticleSet& P, int iat)
    {
      return Dets[DetID[iat]]->evalGrad(P,iat);
    }

    GradType evalGradSource(ParticleSet& P, ParticleSet &src, int iat)
    {
      GradType G = GradType();
      for (int iz=0; iz < size(); iz++) 
	G += Dets[iz]->evalGradSource (P, src, iat);
      return G;
    }



    inline ValueType logRatio(ParticleSet& P, int iat,
			   ParticleSet::ParticleGradient_t& dG, 
			   ParticleSet::ParticleLaplacian_t& dL) { 
      ValueType r = Dets[DetID[iat]]->ratio(P,iat,dG,dL);
      return evaluateLogAndPhase(r,PhaseValue);
    }
    
    inline void restore(int iat) 
    {
      return Dets[DetID[iat]]->restore(iat);
    }

    inline void acceptMove(ParticleSet& P, int iat) 
    {
      Dets[DetID[iat]]->acceptMove(P,iat);
    }

    inline ValueType
    ratio(ParticleSet& P, int iat) 
    {
      return Dets[DetID[iat]]->ratio(P,iat);
    } 	  

    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) 
    {
      return Dets[DetID[iat]]->update(P,dG,dL,iat);
    }

    OrbitalBasePtr makeClone(ParticleSet& tqp) const;

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
