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
#ifndef QMCPLUSPLUS_SINGLEPARTICLEORBITALSET_H
#define QMCPLUSPLUS_SINGLEPARTICLEORBITALSET_H
#include <vector>
#include "QMCWaveFunctions/DummyBasisSet.h"
namespace qmcplusplus {

/**a set of single-particle orbitals. 
 *
 * This class provides necessary interfaces for SlaterDeterminant<SPOSet>
 * and can be used by any orbital representation that has an one-to-one
 * mapping between the function evaluations and the column indeices of a
 * Dirac determinant.
 * The template parameter is any function with a member function
 * which can provde the value, gradient and laplacian at a point.
 * value_type OT::evaluate(const point_type&, gradient_type&, value_type& ) 
 * Example classes can be found in Numerics/CosineFunction.h
 */
template<class OT>
struct SingleParticleOrbitalSet {

  ///the type of single-particle orbtials 
  typedef OT                      SPOrbital_t;
  typedef vector<OT*>             SPOContainer_t;
  typedef typename OT::value_type value_type;
  typedef DummyBasisSet           BasisSet_t;

  SPOContainer_t Phi;

  ///constructor
  SingleParticleOrbitalSet(){ }

  /**add a single-particle orbital */
  int add(SPOrbital_t* afunction ) {
    Phi.push_back(afunction);
    return Phi.size()-1;
  }

  inline int size() const { return Phi.size();}

  void reset() { for(int i=0; i<Phi.size(); i++) Phi[i]->reset(); }
  void resetTargetParticleSet(ParticleSet& P) { }

  inline value_type
  evaluate(const ParticleSet& P, int iat, int jorb) {
    return Phi[jorb]->evaluate(P.R[iat]);
  }

  template<class VV>
  inline void 
  evaluate(const ParticleSet& P, int iat, VV& psi) {
    vector<SPOrbital_t*>::iterator it(Phi.begin()),it_end(Phi.end());
    int j(0);
    while(it != it_end) {
      psi[j]=(*it)->evaluate(P.R[iat]);++it;j++;
    }
  }

  template<class VV, class GV>
  inline void
  evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi) {
    vector<SPOrbital_t*>::iterator it(Phi.begin()),it_end(Phi.end());
    int j(0);
    while(it != it_end) {
      psi[j]=(*it)->evaluate(P.R[iat],dpsi[j],d2psi[j]);++it;j++;
    }
  }

  template<class VM, class GM>
  inline void evaluate(const ParticleSet& P, int first, int last,
		       VM& logdet, GM& dlogdet, VM& d2logdet) {
    int n = last-first;
    int iat = first;
    vector<SPOrbital_t*>::iterator it_end(Phi.end());
    for(int i=0; i<n; i++,iat++) {
      vector<SPOrbital_t*>::iterator it(Phi.begin());
      int j(0);
      while(it != it_end) {
	logdet(j,i)= (*it)->evaluate(P.R[iat], dlogdet(i,j),d2logdet(i,j));
        ++it;++j;
      }
    }
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
